#include "m2424/seal_adapter.hpp"

#include <seal/seal.h>
#include <seal/util/iterator.h>
#include <seal/util/ntt.h>
#include <seal/util/rns.h>
#include <seal/util/uintarith.h>
#include <seal/util/uintarithmod.h>
#include <seal/util/uintarithsmallmod.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace m2424 {

// ---- Plain ----
struct Plain::Impl {
    seal::Plaintext pt;
};

Plain::Plain() : pimpl_(new Impl{}) {}
Plain::~Plain() = default;
Plain::Plain(const Plain& other) : pimpl_(other.pimpl_ ? std::make_unique<Impl>(*other.pimpl_) : nullptr) {}
Plain& Plain::operator=(const Plain& other) {
    if (this != &other) {
        pimpl_ = other.pimpl_ ? std::make_unique<Impl>(*other.pimpl_) : nullptr;
    }
    return *this;
}
Plain::Plain(Plain&&) noexcept = default;
Plain& Plain::operator=(Plain&&) noexcept = default;

// ---- Cipher ----
struct Cipher::Impl {
    seal::Ciphertext ct;
};

Cipher::Cipher() : pimpl_(new Impl{}) {}
Cipher::~Cipher() = default;
Cipher::Cipher(const Cipher& other) : pimpl_(other.pimpl_ ? std::make_unique<Impl>(*other.pimpl_) : nullptr) {}
Cipher& Cipher::operator=(const Cipher& other) {
    if (this != &other) {
        pimpl_ = other.pimpl_ ? std::make_unique<Impl>(*other.pimpl_) : nullptr;
    }
    return *this;
}
Cipher::Cipher(Cipher&&) noexcept = default;
Cipher& Cipher::operator=(Cipher&&) noexcept = default;

// ---- SealAdapter ----
struct SealAdapter::Impl {
    CkksProfile profile{};

    std::shared_ptr<seal::SEALContext> context;
    std::unique_ptr<seal::CKKSEncoder> encoder;
    std::unique_ptr<seal::Evaluator> evaluator;
    std::unique_ptr<seal::Encryptor> encryptor;
    std::unique_ptr<seal::Decryptor> decryptor;

    seal::SecretKey sk;
    seal::PublicKey pk;
    seal::RelinKeys rlk;
    seal::GaloisKeys gk;

    double scale{0.0};
    std::size_t slotCount{0};
    bool has_public{false};
    bool has_secret{false};
    bool has_relin{false};
    bool has_galois{false};
};

SealAdapter::SealAdapter() : pimpl_(new Impl{}) {}
SealAdapter::~SealAdapter() = default;
// non-copyable
SealAdapter::SealAdapter(SealAdapter&&) noexcept = default;
SealAdapter& SealAdapter::operator=(SealAdapter&&) noexcept = default;

static seal::EncryptionParameters make_ckks_parms(const CkksProfile& prof) {
    using namespace seal;
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(prof.polyModulusDegree);
    parms.set_coeff_modulus(CoeffModulus::Create(prof.polyModulusDegree, prof.coeffModulusBits));
    return parms;
}

static bool is_power_of_two(std::size_t value) {
    return value != 0 && (value & (value - 1)) == 0;
}

static void validate_profile_shape(const CkksProfile& profile) {
    if (!is_power_of_two(profile.polyModulusDegree)) {
        throw std::invalid_argument("polyModulusDegree must be a non-zero power of two");
    }
    if (profile.coeffModulusBits.empty()) {
        throw std::invalid_argument("coeffModulusBits must not be empty");
    }
    if (!std::isfinite(profile.scale) || profile.scale <= 0.0) {
        throw std::invalid_argument("scale must be a positive finite value");
    }
    for (int bits : profile.coeffModulusBits) {
        if (bits <= 0) {
            throw std::invalid_argument("coeffModulusBits entries must be positive");
        }
    }
}

static void validate_real_values(const std::vector<double>& vals, std::size_t slotCount, std::size_t profile_slots) {
    if (vals.empty()) {
        throw std::invalid_argument("input vector must not be empty");
    }
    if (vals.size() > slotCount || (profile_slots && vals.size() > profile_slots)) {
        throw std::invalid_argument("input size exceeds configured slots");
    }
    for (double value : vals) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("input values must be finite");
        }
    }
}

static void validate_complex_values(const std::vector<std::complex<double>>& vals,
                                    std::size_t slotCount, std::size_t profile_slots) {
    if (vals.empty()) {
        throw std::invalid_argument("input vector must not be empty");
    }
    if (vals.size() > slotCount || (profile_slots && vals.size() > profile_slots)) {
        throw std::invalid_argument("input size exceeds configured slots");
    }
    for (const auto& value : vals) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            throw std::invalid_argument("input values must be finite");
        }
    }
}

template <class T>
static std::size_t serialized_size_of(const T& value) {
    std::ostringstream out(std::ios::binary);
    value.save(out);
    return static_cast<std::size_t>(out.tellp());
}

template <class T>
static SerializedBuffer serialize_to_buffer(const T& value) {
    std::ostringstream out(std::ios::binary);
    value.save(out);
    const auto bytes = out.str();
    return SerializedBuffer(bytes.begin(), bytes.end());
}

static void require_non_empty_buffer(const SerializedBuffer& buffer, const char* name) {
    if (buffer.empty()) {
        throw std::invalid_argument(std::string(name) + " buffer must not be empty");
    }
}

template <class T>
static void load_from_buffer(T& value, const seal::SEALContext& context,
                             const SerializedBuffer& buffer, const char* name) {
    require_non_empty_buffer(buffer, name);
    value.load(context, reinterpret_cast<const seal::seal_byte*>(buffer.data()), buffer.size());
}

SealAdapter SealAdapter::create(const CkksProfile& profile) {
    validate_profile_shape(profile);

    SealAdapter a;
    a.pimpl_->profile = profile;

    auto parms = make_ckks_parms(profile);
    a.pimpl_->context = std::make_shared<seal::SEALContext>(parms, /*expand_mod_chain*/ true);
    if (!a.pimpl_->context->parameters_set()) {
        throw std::invalid_argument("invalid CKKS encryption parameters");
    }
    a.pimpl_->encoder = std::make_unique<seal::CKKSEncoder>(*a.pimpl_->context);
    a.pimpl_->evaluator = std::make_unique<seal::Evaluator>(*a.pimpl_->context);
    a.pimpl_->scale = profile.scale;
    a.pimpl_->slotCount = a.pimpl_->encoder->slot_count();
    if (profile.slots > a.pimpl_->slotCount) {
        throw std::invalid_argument("configured slots exceed CKKS slotCount");
    }
    return a;
}

void SealAdapter::generateKeys(bool needRelin, bool needGalois) {
    if (!pimpl_->context) throw std::runtime_error("SealAdapter not initialized");
    seal::KeyGenerator generateKeys(*pimpl_->context);
    pimpl_->sk = generateKeys.secret_key();
    generateKeys.create_public_key(pimpl_->pk);
    pimpl_->has_relin = false;
    pimpl_->has_galois = false;
    pimpl_->has_public = true;
    pimpl_->has_secret = true;
    if (needRelin) {
        generateKeys.create_relin_keys(pimpl_->rlk);
        pimpl_->has_relin = true;
    }
    if (needGalois) {
        generateKeys.create_galois_keys(pimpl_->gk);
        pimpl_->has_galois = true;
    }

    pimpl_->encryptor = std::make_unique<seal::Encryptor>(*pimpl_->context, pimpl_->pk);
    pimpl_->decryptor = std::make_unique<seal::Decryptor>(*pimpl_->context, pimpl_->sk);
}

void SealAdapter::generateKeys(const std::vector<int>& rotationSteps, bool needRelin) {
    if (!pimpl_->context) throw std::runtime_error("SealAdapter not initialized");
    seal::KeyGenerator generateKeys(*pimpl_->context);
    pimpl_->sk = generateKeys.secret_key();
    generateKeys.create_public_key(pimpl_->pk);
    pimpl_->has_relin = false;
    pimpl_->has_galois = false;
    pimpl_->has_public = true;
    pimpl_->has_secret = true;
    if (needRelin) {
        generateKeys.create_relin_keys(pimpl_->rlk);
        pimpl_->has_relin = true;
    }
    if (!rotationSteps.empty()) {
        generateKeys.create_galois_keys(rotationSteps, pimpl_->gk);
        pimpl_->has_galois = true;
    }

    pimpl_->encryptor = std::make_unique<seal::Encryptor>(*pimpl_->context, pimpl_->pk);
    pimpl_->decryptor = std::make_unique<seal::Decryptor>(*pimpl_->context, pimpl_->sk);
}

std::size_t SealAdapter::slotCount() const {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    return pimpl_->slotCount;
}

Plain SealAdapter::encode(const std::vector<double>& vals) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_real_values(vals, pimpl_->slotCount, pimpl_->profile.slots);
    Plain out;
    pimpl_->encoder->encode(vals, pimpl_->scale, out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encodeComplex(const std::vector<std::complex<double>>& vals) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_complex_values(vals, pimpl_->slotCount, pimpl_->profile.slots);
    Plain out;
    pimpl_->encoder->encode(vals, pimpl_->scale, out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encodeFor(const std::vector<double>& vals, const Cipher& cipher) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_real_values(vals, pimpl_->slotCount, pimpl_->profile.slots);
    Plain out;
    pimpl_->encoder->encode(vals, cipher.pimpl_->ct.parms_id(), cipher.pimpl_->ct.scale(), out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encodeComplexFor(const std::vector<std::complex<double>>& vals, const Cipher& cipher) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_complex_values(vals, pimpl_->slotCount, pimpl_->profile.slots);
    Plain out;
    pimpl_->encoder->encode(vals, cipher.pimpl_->ct.parms_id(), cipher.pimpl_->ct.scale(), out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encodeComplexAtScaleFor(const std::vector<std::complex<double>>& vals,
                                                double scale,
                                                const Cipher& cipher) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_complex_values(vals, pimpl_->slotCount, pimpl_->profile.slots);
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::invalid_argument("scale must be a positive finite value");
    }
    Plain out;
    pimpl_->encoder->encode(vals, cipher.pimpl_->ct.parms_id(), scale, out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encodeScalarFor(double value, const Cipher& cipher) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    if (!std::isfinite(value)) throw std::invalid_argument("scalar must be finite");
    Plain out;
    pimpl_->encoder->encode(value, cipher.pimpl_->ct.parms_id(), cipher.pimpl_->ct.scale(), out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encodeScalarAtScaleFor(double value, double scale, const Cipher& cipher) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    if (!std::isfinite(value)) throw std::invalid_argument("scalar must be finite");
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw std::invalid_argument("scale must be a positive finite value");
    }
    Plain out;
    pimpl_->encoder->encode(value, cipher.pimpl_->ct.parms_id(), scale, out.pimpl_->pt);
    return out;
}

Cipher SealAdapter::encrypt(const Plain& plain) {
    if (!pimpl_->has_public || !pimpl_->encryptor) throw std::runtime_error("public key not loaded");
    Cipher out;
    pimpl_->encryptor->encrypt(plain.pimpl_->pt, out.pimpl_->ct);
    return out;
}

Plain SealAdapter::decrypt(const Cipher& cipher) {
    if (!pimpl_->has_secret || !pimpl_->decryptor) throw std::runtime_error("secret key not loaded");
    Plain out;
    pimpl_->decryptor->decrypt(cipher.pimpl_->ct, out.pimpl_->pt);
    return out;
}

std::vector<double> SealAdapter::decode(const Plain& plain) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    std::vector<double> result;
    pimpl_->encoder->decode(plain.pimpl_->pt, result);
    return result;
}

std::vector<std::complex<double>> SealAdapter::decodeComplex(const Plain& plain) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    std::vector<std::complex<double>> result;
    pimpl_->encoder->decode(plain.pimpl_->pt, result);
    return result;
}

Cipher SealAdapter::add(const Cipher& a, const Cipher& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->add(a.pimpl_->ct, b.pimpl_->ct, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::sub(const Cipher& a, const Cipher& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->sub(a.pimpl_->ct, b.pimpl_->ct, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::addPlain(const Cipher& a, const Plain& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->add_plain(a.pimpl_->ct, b.pimpl_->pt, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::subPlain(const Cipher& a, const Plain& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->sub_plain(a.pimpl_->ct, b.pimpl_->pt, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::multiplyPlain(const Cipher& a, const Plain& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->multiply_plain(a.pimpl_->ct, b.pimpl_->pt, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::multiply(const Cipher& a, const Cipher& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->multiply(a.pimpl_->ct, b.pimpl_->ct, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::relinearize(const Cipher& cipher) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    Cipher out = cipher;
    pimpl_->evaluator->relinearize_inplace(out.pimpl_->ct, pimpl_->rlk);
    return out;
}

Cipher SealAdapter::rescaleToNext(const Cipher& cipher) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out = cipher;
    pimpl_->evaluator->rescale_to_next_inplace(out.pimpl_->ct);
    return out;
}


Cipher SealAdapter::modSwitchTo(const Cipher& cipher, const Cipher& target) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out = cipher;
    pimpl_->evaluator->mod_switch_to_inplace(out.pimpl_->ct, target.pimpl_->ct.parms_id());
    return out;
}

Cipher SealAdapter::alignForAddition(const Cipher& cipher, const Cipher& target) {
    Cipher out = modSwitchTo(cipher, target);
    const double source_scale = out.pimpl_->ct.scale();
    const double target_scale = target.pimpl_->ct.scale();
    if (!std::isfinite(source_scale) || !std::isfinite(target_scale) || source_scale <= 0.0 || target_scale <= 0.0) {
        throw std::runtime_error("cannot align ciphertexts for addition: scale must be positive and finite");
    }
    const double relative_error = std::fabs(source_scale - target_scale) / std::max(source_scale, target_scale);
    if (relative_error > 1e-2) {
        throw std::runtime_error("cannot align ciphertexts for addition: relative scale mismatch is too large");
    }
    out.pimpl_->ct.scale() = target_scale;
    return out;
}

Cipher SealAdapter::rotate(const Cipher& c, int steps) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    Cipher out;
    pimpl_->evaluator->rotate_vector(c.pimpl_->ct, steps, pimpl_->gk, out.pimpl_->ct);
    return out;
}

std::size_t SealAdapter::serializedSize(const Cipher& cipher) const {
    return serialized_size_of(cipher.pimpl_->ct);
}

CipherInfo SealAdapter::info(const Cipher& cipher) const {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    const auto context_data = pimpl_->context->get_context_data(cipher.pimpl_->ct.parms_id());
    if (!context_data) throw std::runtime_error("ciphertext parameters are not valid for this context");
    double coeffModulusLog2 = 0.0;
    for (const auto& modulus : context_data->parms().coeff_modulus()) {
        coeffModulusLog2 += std::log2(static_cast<double>(modulus.value()));
    }
    return CipherInfo{
        cipher.pimpl_->ct.scale(),
        context_data->chain_index(),
        context_data->parms().coeff_modulus().size(),
        cipher.pimpl_->ct.size(),
        coeffModulusLog2
    };
}

double SealAdapter::scale(const Cipher& cipher) const {
    return info(cipher).scale;
}

double SealAdapter::coeffModulusLog2(const Cipher& cipher) const {
    return info(cipher).coeffModulusLog2;
}

std::vector<int> SealAdapter::coeffModulusBits() const {
    return pimpl_->profile.coeffModulusBits;
}

std::size_t SealAdapter::chainIndex(const Cipher& cipher) const {
    return info(cipher).chainIndex;
}

std::size_t SealAdapter::coeffModulusSize(const Cipher& cipher) const {
    return info(cipher).coeffModulusSize;
}

std::size_t SealAdapter::publicKeySize() const {
    if (!pimpl_->has_public) throw std::runtime_error("public key not loaded");
    return serialized_size_of(pimpl_->pk);
}

std::size_t SealAdapter::relinKeysSize() const {
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    return serialized_size_of(pimpl_->rlk);
}

std::size_t SealAdapter::galoisKeysSize() const {
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    return serialized_size_of(pimpl_->gk);
}

bool SealAdapter::hasRelinKeys() const noexcept {
    return pimpl_ && pimpl_->has_relin;
}

bool SealAdapter::hasRotationKeys(const std::vector<int>& rotationSteps) const {
    if (rotationSteps.empty()) {
        return true;
    }
    if (!pimpl_ || !pimpl_->context || !pimpl_->has_galois) {
        return false;
    }
    const auto key_context = pimpl_->context->key_context_data();
    if (!key_context || !key_context->galois_tool()) {
        return false;
    }
    for (int step : rotationSteps) {
        if (step == 0 || !pimpl_->gk.has_key(key_context->galois_tool()->get_elt_from_step(step))) {
            return false;
        }
    }
    return true;
}

SerializedBuffer SealAdapter::savePublicKey() const {
    if (!pimpl_->has_public) throw std::runtime_error("public key not loaded");
    return serialize_to_buffer(pimpl_->pk);
}

SerializedBuffer SealAdapter::saveSecretKey() const {
    if (!pimpl_->has_secret) throw std::runtime_error("secret key not loaded");
    return serialize_to_buffer(pimpl_->sk);
}

SerializedBuffer SealAdapter::saveRelinKeys() const {
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    return serialize_to_buffer(pimpl_->rlk);
}

SerializedBuffer SealAdapter::saveGaloisKeys() const {
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    return serialize_to_buffer(pimpl_->gk);
}

SerializedBuffer SealAdapter::saveCipher(const Cipher& cipher) const {
    return serialize_to_buffer(cipher.pimpl_->ct);
}

void SealAdapter::loadPublicKey(const SerializedBuffer& buffer) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    load_from_buffer(pimpl_->pk, *pimpl_->context, buffer, "public key");
    pimpl_->encryptor = std::make_unique<seal::Encryptor>(*pimpl_->context, pimpl_->pk);
    pimpl_->has_public = true;
}

void SealAdapter::loadSecretKey(const SerializedBuffer& buffer) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    load_from_buffer(pimpl_->sk, *pimpl_->context, buffer, "secret key");
    pimpl_->decryptor = std::make_unique<seal::Decryptor>(*pimpl_->context, pimpl_->sk);
    pimpl_->has_secret = true;
}

void SealAdapter::loadRelinKeys(const SerializedBuffer& buffer) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    load_from_buffer(pimpl_->rlk, *pimpl_->context, buffer, "relin keys");
    pimpl_->has_relin = true;
}

void SealAdapter::loadGaloisKeys(const SerializedBuffer& buffer) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    load_from_buffer(pimpl_->gk, *pimpl_->context, buffer, "galois keys");
    pimpl_->has_galois = true;
}

Cipher SealAdapter::loadCipher(const SerializedBuffer& buffer) const {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    Cipher out;
    load_from_buffer(out.pimpl_->ct, *pimpl_->context, buffer, "ciphertext");
    return out;
}

} // namespace m2424
