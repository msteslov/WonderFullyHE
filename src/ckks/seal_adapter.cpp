#include "m2424/seal_adapter.hpp"

#include <seal/seal.h>
#include <seal/util/iterator.h>
#include <seal/util/ntt.h>
#include <seal/util/rns.h>
#include <seal/util/uintarith.h>
#include <seal/util/uintarithmod.h>
#include <seal/util/uintarithsmallmod.h>
#include <seal/util/uintcore.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <limits>

namespace m2424 {
namespace {

std::size_t significantBitCount(const std::uint64_t* limbs, std::size_t limbCount) {
    while (limbCount > 0 && limbs[limbCount - 1] == 0) {
        --limbCount;
    }
    if (limbCount == 0) {
        return 0;
    }
    return (limbCount - 1) * 64
        + static_cast<std::size_t>(64 - __builtin_clzll(limbs[limbCount - 1]));
}

bool bitAt(const std::uint64_t* limbs, std::size_t bit) {
    return ((limbs[bit / 64] >> (bit % 64)) & 1U) != 0;
}

double roundedUnsignedIntegerToDouble(const std::uint64_t* limbs, std::size_t limbCount) {
    constexpr std::size_t precision = 53;
    const std::size_t bitCount = significantBitCount(limbs, limbCount);
    if (bitCount == 0) return 0.0;
    const std::size_t retained = std::min(bitCount, precision);
    const std::size_t discarded = bitCount - retained;
    std::uint64_t significand = 0;
    for (std::size_t offset = 0; offset < retained; ++offset) {
        significand = (significand << 1)
            | static_cast<std::uint64_t>(bitAt(limbs, bitCount - 1 - offset));
    }
    if (discarded != 0) {
        const bool roundBit = bitAt(limbs, discarded - 1);
        bool sticky = false;
        for (std::size_t bit = 0; bit + 1 < discarded && !sticky; ++bit) {
            sticky = bitAt(limbs, bit);
        }
        if (roundBit && (sticky || (significand & 1U))) {
            ++significand;
            if (significand == (std::uint64_t{1} << precision)) {
                significand >>= 1;
                return std::ldexp(static_cast<double>(significand),
                                  static_cast<int>(discarded + 1));
            }
        }
    }
    return std::ldexp(static_cast<double>(significand), static_cast<int>(discarded));
}

double scaledIntegerToDouble(const std::uint64_t* limbs, std::size_t limbCount,
                             bool negative, double scale) {
    constexpr double maxOracleConversionError = 1e-12;
    const std::size_t bitCount = significantBitCount(limbs, limbCount);
    const double magnitude = roundedUnsignedIntegerToDouble(limbs, limbCount);
    const double result = (negative ? -magnitude : magnitude) / scale;
    const std::size_t discarded = bitCount > 53 ? bitCount - 53 : 0;
    const double integerRoundingBound =
        discarded == 0 ? 0.0 : std::ldexp(0.5, static_cast<int>(discarded)) / scale;
    const double divisionRoundingBound =
        std::abs(result) * std::numeric_limits<double>::epsilon();
    if (integerRoundingBound + divisionRoundingBound > maxOracleConversionError) {
        throw std::overflow_error("raw coefficient conversion exceeds oracle precision budget");
    }
    return result;
}

} // namespace

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

RaisedCipher::RaisedCipher(Cipher&& cipher, std::size_t sourceCoeffModulusSize)
    : cipher_(std::move(cipher)), sourceCoeffModulusSize_(sourceCoeffModulusSize) {}
RaisedCipher::RaisedCipher(RaisedCipher&&) noexcept = default;
RaisedCipher& RaisedCipher::operator=(RaisedCipher&&) noexcept = default;

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

std::array<std::uint64_t, 4> SealAdapter::contextFingerprint() const {
    if (!pimpl_->context || !pimpl_->context->key_context_data()) {
        throw std::runtime_error("SEALContext not initialized");
    }
    return pimpl_->context->key_context_data()->parms_id();
}

std::array<std::uint64_t, 4>
SealAdapter::parmsFingerprint(const RaisedCipher& cipher) const {
    if (!cipher.cipher_.pimpl_ || !pimpl_->context
        || !pimpl_->context->get_context_data(cipher.cipher_.pimpl_->ct.parms_id())) {
        throw std::invalid_argument("RaisedCipher does not belong to this SEAL context");
    }
    return cipher.cipher_.pimpl_->ct.parms_id();
}

Plain SealAdapter::encodeComplexAtKeyScale(
        const std::vector<std::complex<double>>& values, double scale) {
    if (!pimpl_->encoder || !pimpl_->context) {
        throw std::runtime_error("CKKS encoder not initialized");
    }
    if (values.empty() || values.size() > pimpl_->slotCount
        || (pimpl_->profile.slots && values.size() > pimpl_->profile.slots)) {
        throw std::invalid_argument("invalid CKKS preparation vector size");
    }
    for (const auto& value : values) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            throw std::invalid_argument("CKKS preparation values must be finite");
        }
    }
    Plain result;
    pimpl_->encoder->encode(
        values, pimpl_->context->key_parms_id(), scale, result.pimpl_->pt);
    return result;
}

double SealAdapter::rescalePlaintextScaleAtChainIndex(std::size_t chainIndex) const {
    if (!pimpl_->context) {
        throw std::runtime_error("SEALContext not initialized");
    }
    auto contextData = pimpl_->context->first_context_data();
    while (contextData && contextData->chain_index() != chainIndex) {
        contextData = contextData->next_context_data();
    }
    if (!contextData || contextData->parms().coeff_modulus().empty()) {
        throw std::invalid_argument("unknown CKKS chain index");
    }
    return std::exp2(static_cast<double>(
        contextData->parms().coeff_modulus().back().bit_count()));
}

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

Plain SealAdapter::encodeScalarRnsAtScaleFor(const std::vector<std::uint64_t>& residues,
                                             double scale, const Cipher& reference,
                                             std::size_t levelsConsumed) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    if (!std::isfinite(scale) || scale <= 0.0)
        throw std::invalid_argument("scale must be a positive finite value");
    auto data = pimpl_->context->get_context_data(reference.pimpl_->ct.parms_id());
    if (!data) throw std::runtime_error("reference ciphertext parameters are not valid for this context");
    for (std::size_t level = 0; level < levelsConsumed; ++level) {
        data = data->next_context_data();
        if (!data) throw std::invalid_argument("prepared scalar level exceeds modulus chain");
    }
    const auto& moduli = data->parms().coeff_modulus();
    if (residues.size() != moduli.size())
        throw std::invalid_argument("prepared scalar RNS residue count mismatch");
    if (static_cast<int>(std::log2(scale)) >= data->total_coeff_modulus_bit_count())
        throw std::invalid_argument("prepared scalar scale out of bounds");
    const std::size_t coefficientCount = data->parms().poly_modulus_degree();
    Plain out;
    out.pimpl_->pt.parms_id() = seal::parms_id_zero;
    out.pimpl_->pt.resize(coefficientCount * moduli.size());
    for (std::size_t limb = 0; limb < moduli.size(); ++limb) {
        if (residues[limb] >= moduli[limb].value())
            throw std::invalid_argument("prepared scalar RNS residue is not reduced");
        std::fill_n(out.pimpl_->pt.data() + limb * coefficientCount,
                    coefficientCount, residues[limb]);
    }
    out.pimpl_->pt.parms_id() = data->parms_id();
    out.pimpl_->pt.scale() = scale;
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

std::vector<double> SealAdapter::decryptRaisedCoefficientsAtRaisedModulus(const RaisedCipher& cipher) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    const auto contextData = pimpl_->context->get_context_data(cipher.cipher_.pimpl_->ct.parms_id());
    if (!contextData) throw std::runtime_error("RaisedCipher parameters are not valid for this context");
    return decryptRaisedCoefficients(cipher, contextData->parms().coeff_modulus().size());
}

std::vector<double> SealAdapter::decryptRaisedCoefficientsAtSourceModulus(const RaisedCipher& cipher) {
    return decryptRaisedCoefficients(cipher, cipher.sourceCoeffModulusSize_);
}

std::vector<double> SealAdapter::decryptRaisedCoefficients(
    const RaisedCipher& cipher,
    std::size_t reconstructionModulusCount) {
    if (!pimpl_->has_secret || !pimpl_->decryptor) throw std::runtime_error("secret key not loaded");
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");

    seal::Plaintext plaintext;
    pimpl_->decryptor->decrypt(cipher.cipher_.pimpl_->ct, plaintext);
    const auto contextData = pimpl_->context->get_context_data(plaintext.parms_id());
    if (!contextData || !plaintext.is_ntt_form()) {
        throw std::runtime_error("raised CKKS plaintext must be in NTT form");
    }
    const auto& targetModuli = contextData->parms().coeff_modulus();
    const std::size_t sourceSize = reconstructionModulusCount;
    const std::size_t degree = contextData->parms().poly_modulus_degree();
    if (sourceSize == 0 || sourceSize > targetModuli.size()) {
        throw std::runtime_error("RaisedCipher has invalid reconstruction RNS base");
    }

    std::vector<std::uint64_t> coefficientRns(plaintext.data(),
                                               plaintext.data() + degree * targetModuli.size());
    seal::util::inverse_ntt_negacyclic_harvey(
        seal::util::RNSIter(coefficientRns.data(), degree),
        targetModuli.size(),
        contextData->small_ntt_tables());
    coefficientRns.resize(degree * sourceSize);

    std::vector<seal::Modulus> sourceModuli(targetModuli.begin(),
                                            targetModuli.begin() + static_cast<std::ptrdiff_t>(sourceSize));
    const auto pool = seal::MemoryManager::GetPool();
    seal::util::RNSBase sourceBase(sourceModuli, pool);
    sourceBase.compose_array(coefficientRns.data(), degree, pool);

    std::vector<std::uint64_t> halfModulus(sourceSize);
    seal::util::half_round_up_uint(sourceBase.base_prod(), sourceSize, halfModulus.data());
    std::vector<double> result(degree);
    std::vector<std::uint64_t> negativeMagnitude(sourceSize);
    for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
        const auto* value = coefficientRns.data() + coefficient * sourceSize;
        const bool negative = seal::util::is_greater_than_or_equal_uint(
            value, halfModulus.data(), sourceSize);
        const auto* magnitudeWords = value;
        if (negative) {
            seal::util::sub_uint(
                sourceBase.base_prod(), value, sourceSize, negativeMagnitude.data());
            magnitudeWords = negativeMagnitude.data();
        }
        const double rounded = scaledIntegerToDouble(
            magnitudeWords, sourceSize, negative, cipher.cipher_.pimpl_->ct.scale());
        if (!std::isfinite(rounded)) {
            throw std::overflow_error("raw raised coefficient does not fit finite double");
        }
        result[coefficient] = rounded;
    }
    return result;
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

RaisedCipher SealAdapter::modRaiseToTop(const Cipher& cipher) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    const auto sourceData = pimpl_->context->get_context_data(cipher.pimpl_->ct.parms_id());
    const auto targetData = pimpl_->context->first_context_data();
    if (!sourceData || !targetData) throw std::runtime_error("ciphertext parameters are not valid for this context");
    const auto& sourceModuli = sourceData->parms().coeff_modulus();
    const auto& targetModuli = targetData->parms().coeff_modulus();
    if (sourceModuli.size() >= targetModuli.size()) {
        throw std::invalid_argument("modRaiseToTop requires a ciphertext below the top level");
    }
    const std::size_t degree = sourceData->parms().poly_modulus_degree();
    for (std::size_t index = 0; index < sourceModuli.size(); ++index) {
        if (sourceModuli[index] != targetModuli[index]) {
            throw std::invalid_argument("source modulus is not a prefix of the top modulus");
        }
    }

    const auto pool = seal::MemoryManager::GetPool();
    seal::util::RNSBase sourceBase(sourceModuli, pool);
    std::vector<seal::Modulus> extraModuli(targetModuli.begin() + static_cast<std::ptrdiff_t>(sourceModuli.size()),
                                           targetModuli.end());
    std::vector<std::uint64_t> halfModulus(sourceBase.size());
    seal::util::half_round_up_uint(sourceBase.base_prod(), sourceBase.size(), halfModulus.data());

    Cipher result;
    result.pimpl_->ct.resize(*pimpl_->context, targetData->parms_id(), cipher.pimpl_->ct.size());
    result.pimpl_->ct.is_ntt_form() = true;
    result.pimpl_->ct.scale() = cipher.pimpl_->ct.scale();
    for (std::size_t component = 0; component < cipher.pimpl_->ct.size(); ++component) {
        std::vector<std::uint64_t> sourceCoefficients(
            cipher.pimpl_->ct.data(component),
            cipher.pimpl_->ct.data(component) + degree * sourceModuli.size());
        seal::util::inverse_ntt_negacyclic_harvey(
            seal::util::RNSIter(sourceCoefficients.data(), degree), sourceModuli.size(), sourceData->small_ntt_tables());
        std::vector<std::uint64_t> centered = sourceCoefficients;
        sourceBase.compose_array(centered.data(), degree, pool);
        std::vector<std::uint64_t> extraCoefficients(degree * extraModuli.size());
        for (std::size_t extraIndex = 0; extraIndex < extraModuli.size(); ++extraIndex) {
            const auto qModuloExtra = seal::util::modulo_uint(sourceBase.base_prod(), sourceBase.size(), extraModuli[extraIndex]);
            for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
                // centered хранит точный CRT-представитель x из [0, q). Новая
                // residue строится непосредственно из него, без приближённого
                // fast base conversion.
                auto value = seal::util::modulo_uint(
                    centered.data() + coefficient * sourceBase.size(), sourceBase.size(), extraModuli[extraIndex]);
                if (seal::util::is_greater_than_or_equal_uint(centered.data() + coefficient * sourceBase.size(),
                                                               halfModulus.data(), sourceBase.size())) {
                    value = seal::util::sub_uint_mod(value, qModuloExtra, extraModuli[extraIndex]);
                }
                extraCoefficients[extraIndex * degree + coefficient] = value;
            }
        }
        auto* output = result.pimpl_->ct.data(component);
        std::copy(sourceCoefficients.begin(), sourceCoefficients.end(), output);
        std::copy(extraCoefficients.begin(), extraCoefficients.end(), output + sourceCoefficients.size());
        seal::util::ntt_negacyclic_harvey(
            seal::util::RNSIter(output, degree), targetModuli.size(), targetData->small_ntt_tables());
#ifdef M2424_ENABLE_MOD_RAISE_CHECKS
        // В тестовой сборке проверяем и сохранённую исходную базу, и каждую
        // новую residue centered lift. Иначе ошибка в знаке новых limbs могла
        // бы пройти проверку prefix-инварианта.
        std::vector<std::uint64_t> verification(output, output + degree * targetModuli.size());
        seal::util::inverse_ntt_negacyclic_harvey(
            seal::util::RNSIter(verification.data(), degree), targetModuli.size(), targetData->small_ntt_tables());
        if (!std::equal(sourceCoefficients.begin(), sourceCoefficients.end(), verification.begin())) {
            throw std::logic_error("modRaise residue invariant failed");
        }
        if (!std::equal(extraCoefficients.begin(), extraCoefficients.end(),
                        verification.begin() + static_cast<std::ptrdiff_t>(sourceCoefficients.size()))) {
            throw std::logic_error("modRaise centered extension invariant failed");
        }
#endif
    }
    return RaisedCipher(std::move(result), sourceModuli.size());
}

Cipher SealAdapter::rotate(const Cipher& c, int steps) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    Cipher out;
    pimpl_->evaluator->rotate_vector(c.pimpl_->ct, steps, pimpl_->gk, out.pimpl_->ct);
    return out;
}

std::vector<Cipher> SealAdapter::rotateManyHoisted(
        const Cipher& cipher, const std::vector<int>& steps) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    std::vector<seal::Ciphertext> nativeOutputs;
    pimpl_->evaluator->rotate_vector_many_hoisted(
        cipher.pimpl_->ct, steps, pimpl_->gk, nativeOutputs);
    std::vector<Cipher> outputs(nativeOutputs.size());
    for (std::size_t index = 0; index < nativeOutputs.size(); ++index) {
        outputs[index].pimpl_->ct = std::move(nativeOutputs[index]);
    }
    return outputs;
}

Cipher SealAdapter::applyBsgsDoubleHoisted(
        const Cipher& cipher, const std::vector<HoistedBsgsGroup>& groups) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    std::vector<seal::HoistedBSGSGroup> nativeGroups;
    nativeGroups.reserve(groups.size());
    for (const auto& group : groups) {
        seal::HoistedBSGSGroup nativeGroup;
        nativeGroup.giant_step = group.giantRotation;
        nativeGroup.terms.reserve(group.terms.size());
        for (const auto& term : group.terms) {
            if (!term.extendedDiagonal) {
                throw std::invalid_argument("null extended BSGS diagonal");
            }
            nativeGroup.terms.push_back({
                term.babyRotation, &term.extendedDiagonal->pimpl_->pt});
        }
        nativeGroups.push_back(std::move(nativeGroup));
    }
    Cipher output;
    pimpl_->evaluator->apply_bsgs_double_hoisted(
        cipher.pimpl_->ct, nativeGroups, pimpl_->gk, output.pimpl_->ct);
    return output;
}

Cipher SealAdapter::conjugate(const Cipher& cipher) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!hasConjugationKey()) throw std::runtime_error("conjugation key not generated");
    Cipher out;
    pimpl_->evaluator->complex_conjugate(cipher.pimpl_->ct, pimpl_->gk, out.pimpl_->ct);
    return out;
}

std::size_t SealAdapter::serializedSize(const Cipher& cipher) const {
    return serialized_size_of(cipher.pimpl_->ct);
}

std::size_t SealAdapter::serializedSize(const Plain& plain) const {
    return serialized_size_of(plain.pimpl_->pt);
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

CipherInfo SealAdapter::info(const RaisedCipher& cipher) const {
    return info(cipher.cipher_);
}

double SealAdapter::scale(const Cipher& cipher) const {
    return info(cipher).scale;
}

Cipher SealAdapter::normalizeScale(const Cipher& cipher, double targetScale) const {
    const double current = scale(cipher);
    if (!std::isfinite(targetScale) || targetScale <= 0.0 || !std::isfinite(current)
        || std::abs(std::log2(current / targetScale)) > 0.1) {
        throw std::invalid_argument("unsafe ciphertext scale normalization");
    }
    Cipher result = cipher;
    result.pimpl_->ct.scale() = targetScale;
    return result;
}

double SealAdapter::coeffModulusLog2(const Cipher& cipher) const {
    return info(cipher).coeffModulusLog2;
}

std::vector<int> SealAdapter::coeffModulusBits() const {
    return pimpl_->profile.coeffModulusBits;
}

std::vector<std::uint64_t> SealAdapter::dataModulusValues() const {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    const auto data = pimpl_->context->first_context_data();
    if (!data) throw std::runtime_error("CKKS data modulus chain is not initialized");
    std::vector<std::uint64_t> values;
    values.reserve(data->parms().coeff_modulus().size());
    for (const auto& modulus : data->parms().coeff_modulus()) values.push_back(modulus.value());
    return values;
}

std::uint64_t SealAdapter::specialKeyModulusValue() const {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    const auto key = pimpl_->context->key_context_data();
    const auto data = pimpl_->context->first_context_data();
    if (!key || !data) throw std::runtime_error("CKKS modulus chain is not initialized");
    const auto& keyModuli = key->parms().coeff_modulus();
    const auto& dataModuli = data->parms().coeff_modulus();
    if (keyModuli.size() != dataModuli.size() + 1)
        throw std::runtime_error("CKKS key modulus does not contain one special prime");
    for (std::size_t index = 0; index < dataModuli.size(); ++index) {
        if (keyModuli[index].value() != dataModuli[index].value())
            throw std::runtime_error("CKKS data modulus is not a prefix of key modulus");
    }
    return keyModuli.back().value();
}

std::vector<std::uint64_t> SealAdapter::coeffModulusValues(const Cipher& cipher) const {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    const auto data = pimpl_->context->get_context_data(cipher.pimpl_->ct.parms_id());
    if (!data) throw std::runtime_error("ciphertext parameters are not valid for this context");
    std::vector<std::uint64_t> values;
    values.reserve(data->parms().coeff_modulus().size());
    for (const auto& modulus : data->parms().coeff_modulus()) values.push_back(modulus.value());
    return values;
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

bool SealAdapter::hasConjugationKey() const {
    if (!pimpl_ || !pimpl_->context || !pimpl_->has_galois) {
        return false;
    }
    const auto keyContext = pimpl_->context->key_context_data();
    return keyContext && keyContext->galois_tool()
        && pimpl_->gk.has_key(keyContext->galois_tool()->get_elt_from_step(0));
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
