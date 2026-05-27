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
    std::size_t slot_count{0};
    bool has_keys{false};
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
    parms.set_poly_modulus_degree(prof.poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(prof.poly_modulus_degree, prof.coeff_modulus_bits));
    return parms;
}

static bool is_power_of_two(std::size_t value) {
    return value != 0 && (value & (value - 1)) == 0;
}

static void validate_profile_shape(const CkksProfile& profile) {
    if (!is_power_of_two(profile.poly_modulus_degree)) {
        throw std::invalid_argument("poly_modulus_degree must be a non-zero power of two");
    }
    if (profile.coeff_modulus_bits.empty()) {
        throw std::invalid_argument("coeff_modulus_bits must not be empty");
    }
    if (!std::isfinite(profile.scale) || profile.scale <= 0.0) {
        throw std::invalid_argument("scale must be a positive finite value");
    }
    for (int bits : profile.coeff_modulus_bits) {
        if (bits <= 0) {
            throw std::invalid_argument("coeff_modulus_bits entries must be positive");
        }
    }
}

static void validate_real_values(const std::vector<double>& vals, std::size_t slot_count, std::size_t profile_slots) {
    if (vals.empty()) {
        throw std::invalid_argument("input vector must not be empty");
    }
    if (vals.size() > slot_count || (profile_slots && vals.size() > profile_slots)) {
        throw std::invalid_argument("input size exceeds configured slots");
    }
    for (double value : vals) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("input values must be finite");
        }
    }
}

static void validate_complex_values(const std::vector<std::complex<double>>& vals,
                                    std::size_t slot_count, std::size_t profile_slots) {
    if (vals.empty()) {
        throw std::invalid_argument("input vector must not be empty");
    }
    if (vals.size() > slot_count || (profile_slots && vals.size() > profile_slots)) {
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
    a.pimpl_->slot_count = a.pimpl_->encoder->slot_count();
    if (profile.slots > a.pimpl_->slot_count) {
        throw std::invalid_argument("configured slots exceed CKKS slot_count");
    }
    return a;
}

void SealAdapter::keygen(bool need_relin, bool need_galois) {
    if (!pimpl_->context) throw std::runtime_error("SealAdapter not initialized");
    seal::KeyGenerator keygen(*pimpl_->context);
    pimpl_->sk = keygen.secret_key();
    keygen.create_public_key(pimpl_->pk);
    pimpl_->has_relin = false;
    pimpl_->has_galois = false;
    pimpl_->has_public = true;
    pimpl_->has_secret = true;
    if (need_relin) {
        keygen.create_relin_keys(pimpl_->rlk);
        pimpl_->has_relin = true;
    }
    if (need_galois) {
        keygen.create_galois_keys(pimpl_->gk);
        pimpl_->has_galois = true;
    }

    pimpl_->encryptor = std::make_unique<seal::Encryptor>(*pimpl_->context, pimpl_->pk);
    pimpl_->decryptor = std::make_unique<seal::Decryptor>(*pimpl_->context, pimpl_->sk);
    pimpl_->has_keys = true;
}

void SealAdapter::keygen(const std::vector<int>& rotation_steps, bool need_relin) {
    if (!pimpl_->context) throw std::runtime_error("SealAdapter not initialized");
    seal::KeyGenerator keygen(*pimpl_->context);
    pimpl_->sk = keygen.secret_key();
    keygen.create_public_key(pimpl_->pk);
    pimpl_->has_relin = false;
    pimpl_->has_galois = false;
    pimpl_->has_public = true;
    pimpl_->has_secret = true;
    if (need_relin) {
        keygen.create_relin_keys(pimpl_->rlk);
        pimpl_->has_relin = true;
    }
    if (!rotation_steps.empty()) {
        keygen.create_galois_keys(rotation_steps, pimpl_->gk);
        pimpl_->has_galois = true;
    }

    pimpl_->encryptor = std::make_unique<seal::Encryptor>(*pimpl_->context, pimpl_->pk);
    pimpl_->decryptor = std::make_unique<seal::Decryptor>(*pimpl_->context, pimpl_->sk);
    pimpl_->has_keys = true;
}

std::size_t SealAdapter::slot_count() const {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    return pimpl_->slot_count;
}

Plain SealAdapter::encode(const std::vector<double>& vals) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_real_values(vals, pimpl_->slot_count, pimpl_->profile.slots);
    Plain out;
    pimpl_->encoder->encode(vals, pimpl_->scale, out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encode_complex(const std::vector<std::complex<double>>& vals) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_complex_values(vals, pimpl_->slot_count, pimpl_->profile.slots);
    Plain out;
    pimpl_->encoder->encode(vals, pimpl_->scale, out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encode_like(const std::vector<double>& vals, const Cipher& cipher) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_real_values(vals, pimpl_->slot_count, pimpl_->profile.slots);
    Plain out;
    pimpl_->encoder->encode(vals, cipher.pimpl_->ct.parms_id(), cipher.pimpl_->ct.scale(), out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encode_complex_like(const std::vector<std::complex<double>>& vals, const Cipher& cipher) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    validate_complex_values(vals, pimpl_->slot_count, pimpl_->profile.slots);
    Plain out;
    pimpl_->encoder->encode(vals, cipher.pimpl_->ct.parms_id(), cipher.pimpl_->ct.scale(), out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encode_scalar_like(double value, const Cipher& cipher) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    if (!std::isfinite(value)) throw std::invalid_argument("scalar must be finite");
    Plain out;
    pimpl_->encoder->encode(value, cipher.pimpl_->ct.parms_id(), cipher.pimpl_->ct.scale(), out.pimpl_->pt);
    return out;
}

Plain SealAdapter::encode_scalar_at_scale_like(double value, double scale, const Cipher& cipher) {
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

std::vector<std::complex<double>> SealAdapter::decode_complex(const Plain& plain) {
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

Cipher SealAdapter::add_plain(const Cipher& a, const Plain& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->add_plain(a.pimpl_->ct, b.pimpl_->pt, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::sub_plain(const Cipher& a, const Plain& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->sub_plain(a.pimpl_->ct, b.pimpl_->pt, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::mul_plain(const Cipher& a, const Plain& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->multiply_plain(a.pimpl_->ct, b.pimpl_->pt, out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::mul_plain_rescale(const Cipher& a, const Plain& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out;
    pimpl_->evaluator->multiply_plain(a.pimpl_->ct, b.pimpl_->pt, out.pimpl_->ct);
    pimpl_->evaluator->rescale_to_next_inplace(out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::mul_relin_rescale(const Cipher& a, const Cipher& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    Cipher out;
    pimpl_->evaluator->multiply(a.pimpl_->ct, b.pimpl_->ct, out.pimpl_->ct);
    pimpl_->evaluator->relinearize_inplace(out.pimpl_->ct, pimpl_->rlk);
    pimpl_->evaluator->rescale_to_next_inplace(out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::mod_raise_to_first(const Cipher& cipher) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    const auto source_context_data = pimpl_->context->get_context_data(cipher.pimpl_->ct.parms_id());
    if (!source_context_data) {
        throw std::runtime_error("ciphertext parameters are not valid for this context");
    }
    const auto target_context_data = pimpl_->context->first_context_data();
    if (!target_context_data) {
        throw std::runtime_error("target context data is not available");
    }
    if (cipher.pimpl_->ct.parms_id() == target_context_data->parms_id()) {
        return cipher;
    }
    if (!cipher.pimpl_->ct.is_ntt_form()) {
        throw std::invalid_argument("CKKS ciphertext must be in NTT form for modulus raising");
    }

    const auto& source_parms = source_context_data->parms();
    const auto& target_parms = target_context_data->parms();
    if (source_parms.scheme() != seal::scheme_type::ckks || target_parms.scheme() != seal::scheme_type::ckks) {
        throw std::invalid_argument("mod_raise_to_first supports CKKS only");
    }

    const std::size_t coeff_count = source_parms.poly_modulus_degree();
    const auto& source_moduli = source_parms.coeff_modulus();
    const auto& target_moduli = target_parms.coeff_modulus();
    const std::size_t source_modulus_size = source_moduli.size();
    const std::size_t target_modulus_size = target_moduli.size();
    if (source_modulus_size >= target_modulus_size) {
        throw std::invalid_argument("ciphertext is not below the first modulus level");
    }
    if (coeff_count != target_parms.poly_modulus_degree()) {
        throw std::invalid_argument("source and target polynomial degrees differ");
    }
    for (std::size_t i = 0; i < source_modulus_size; ++i) {
        if (source_moduli[i] != target_moduli[i]) {
            throw std::invalid_argument("source modulus chain is not a prefix of target chain");
        }
    }

    auto pool = seal::MemoryManager::GetPool();
    seal::util::RNSBase source_base(source_moduli, pool);
    std::vector<seal::Modulus> extension_moduli(
        target_moduli.begin() + static_cast<std::ptrdiff_t>(source_modulus_size), target_moduli.end());
    seal::util::RNSBase extension_base(extension_moduli, pool);
    seal::util::BaseConverter source_to_extension(source_base, extension_base, pool);
    std::vector<std::uint64_t> half_source_modulus(source_modulus_size);
    seal::util::right_shift_uint(source_base.base_prod(), 1, source_modulus_size, half_source_modulus.data());
    std::vector<std::uint64_t> source_modulus_mod_extension(extension_moduli.size());
    for (std::size_t i = 0; i < extension_moduli.size(); ++i) {
        source_modulus_mod_extension[i] =
            seal::util::modulo_uint(source_base.base_prod(), source_modulus_size, extension_moduli[i]);
    }

    Cipher out;
    out.pimpl_->ct.resize(*pimpl_->context, target_context_data->parms_id(), cipher.pimpl_->ct.size());
    out.pimpl_->ct.is_ntt_form() = true;
    out.pimpl_->ct.scale() = cipher.pimpl_->ct.scale();
    out.pimpl_->ct.correction_factor() = cipher.pimpl_->ct.correction_factor();

    std::vector<std::uint64_t> source_normal(source_modulus_size * coeff_count);
    std::vector<std::uint64_t> extension_normal(extension_moduli.size() * coeff_count);
    std::vector<std::uint64_t> composed(source_modulus_size * coeff_count);

    for (std::size_t poly_index = 0; poly_index < cipher.pimpl_->ct.size(); ++poly_index) {
        const auto* source_poly = cipher.pimpl_->ct.data(poly_index);
        auto* target_poly = out.pimpl_->ct.data(poly_index);

        std::copy_n(source_poly, source_normal.size(), source_normal.begin());
        seal::util::RNSIter source_normal_iter(source_normal.data(), coeff_count);
        seal::util::inverse_ntt_negacyclic_harvey(
            source_normal_iter, source_modulus_size, source_context_data->small_ntt_tables());

        seal::util::RNSIter extension_normal_iter(extension_normal.data(), coeff_count);
        source_to_extension.fast_convert_array(source_normal_iter, extension_normal_iter, pool);

        composed = source_normal;
        source_base.compose_array(composed.data(), coeff_count, pool);
        for (std::size_t coeff_index = 0; coeff_index < coeff_count; ++coeff_index) {
            const auto* coeff_value = composed.data() + coeff_index * source_modulus_size;
            if (!seal::util::is_greater_than_uint(coeff_value,
                                                  half_source_modulus.data(),
                                                  source_modulus_size)) {
                continue;
            }
            for (std::size_t extension_index = 0; extension_index < extension_moduli.size(); ++extension_index) {
                auto& value = extension_normal[extension_index * coeff_count + coeff_index];
                value = seal::util::sub_uint_mod(value,
                                                 source_modulus_mod_extension[extension_index],
                                                 extension_moduli[extension_index]);
            }
        }

        seal::util::RNSIter target_iter(target_poly, coeff_count);
        for (std::size_t modulus_index = 0; modulus_index < source_modulus_size; ++modulus_index) {
            std::copy_n(source_normal.data() + modulus_index * coeff_count,
                        coeff_count,
                        target_iter[modulus_index]);
        }
        for (std::size_t extension_index = 0; extension_index < extension_moduli.size(); ++extension_index) {
            std::copy_n(extension_normal.data() + extension_index * coeff_count,
                        coeff_count,
                        target_iter[source_modulus_size + extension_index]);
        }
        seal::util::ntt_negacyclic_harvey(
            target_iter, target_modulus_size, target_context_data->small_ntt_tables());
    }

    return out;
}

Cipher SealAdapter::multiply_decoded_value(const Cipher& cipher, double multiplier) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    if (!std::isfinite(multiplier) || multiplier <= 0.0) {
        throw std::invalid_argument("decoded value multiplier must be positive and finite");
    }
    const auto context_data = pimpl_->context->get_context_data(cipher.pimpl_->ct.parms_id());
    if (!context_data) {
        throw std::runtime_error("ciphertext parameters are not valid for this context");
    }

    Cipher out = cipher;
    const double old_scale = out.pimpl_->ct.scale();
    const double new_scale = old_scale / multiplier;
    if (!std::isfinite(old_scale) || old_scale <= 0.0 || !std::isfinite(new_scale) || new_scale <= 0.0) {
        throw std::runtime_error("cannot adjust decoded value: invalid ciphertext scale");
    }

    const double log2_scale = std::log2(new_scale);
    const int total_bits = context_data->total_coeff_modulus_bit_count();
    if (log2_scale < 0.0 || log2_scale > static_cast<double>(total_bits - 2)) {
        throw std::runtime_error("cannot adjust decoded value: resulting scale is outside modulus capacity");
    }

    out.pimpl_->ct.scale() = new_scale;
    return out;
}

Cipher SealAdapter::mod_switch_to(const Cipher& cipher, const Cipher& target) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    Cipher out = cipher;
    pimpl_->evaluator->mod_switch_to_inplace(out.pimpl_->ct, target.pimpl_->ct.parms_id());
    return out;
}

Cipher SealAdapter::match_level_and_scale(const Cipher& cipher, const Cipher& target) {
    Cipher out = mod_switch_to(cipher, target);
    const double source_scale = out.pimpl_->ct.scale();
    const double target_scale = target.pimpl_->ct.scale();
    if (!std::isfinite(source_scale) || !std::isfinite(target_scale) || source_scale <= 0.0 || target_scale <= 0.0) {
        throw std::runtime_error("cannot match ciphertext scales: scale must be positive and finite");
    }
    const double relative_error = std::fabs(source_scale - target_scale) / std::max(source_scale, target_scale);
    if (relative_error > 1e-3) {
        throw std::runtime_error("cannot match ciphertext scales: relative mismatch is too large");
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

std::size_t SealAdapter::serialized_size(const Cipher& cipher) const {
    return serialized_size_of(cipher.pimpl_->ct);
}

CipherInfo SealAdapter::info(const Cipher& cipher) const {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    const auto context_data = pimpl_->context->get_context_data(cipher.pimpl_->ct.parms_id());
    if (!context_data) throw std::runtime_error("ciphertext parameters are not valid for this context");
    double coeff_modulus_log2 = 0.0;
    for (const auto& modulus : context_data->parms().coeff_modulus()) {
        coeff_modulus_log2 += std::log2(static_cast<double>(modulus.value()));
    }
    return CipherInfo{
        cipher.pimpl_->ct.scale(),
        context_data->chain_index(),
        context_data->parms().coeff_modulus().size(),
        cipher.pimpl_->ct.size(),
        coeff_modulus_log2
    };
}

double SealAdapter::scale(const Cipher& cipher) const {
    return info(cipher).scale;
}

double SealAdapter::coeff_modulus_log2(const Cipher& cipher) const {
    return info(cipher).coeff_modulus_log2;
}

double SealAdapter::bootstrap_period_log2(const Cipher& cipher) const {
    const auto cipher_info = info(cipher);
    if (!std::isfinite(cipher_info.scale) || cipher_info.scale <= 0.0) {
        throw std::runtime_error("ciphertext scale must be positive and finite");
    }
    return cipher_info.coeff_modulus_log2 - std::log2(cipher_info.scale);
}

double SealAdapter::bootstrap_period(const Cipher& cipher) const {
    const double period_log2 = bootstrap_period_log2(cipher);
    if (period_log2 > 1023.0) {
        throw std::runtime_error("bootstrap period exceeds double range");
    }
    return std::exp2(period_log2);
}

std::vector<int> SealAdapter::coeff_modulus_bits() const {
    return pimpl_->profile.coeff_modulus_bits;
}

std::size_t SealAdapter::chain_index(const Cipher& cipher) const {
    return info(cipher).chain_index;
}

std::size_t SealAdapter::coeff_modulus_size(const Cipher& cipher) const {
    return info(cipher).coeff_modulus_size;
}

std::size_t SealAdapter::public_key_size() const {
    if (!pimpl_->has_public) throw std::runtime_error("public key not loaded");
    return serialized_size_of(pimpl_->pk);
}

std::size_t SealAdapter::relin_keys_size() const {
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    return serialized_size_of(pimpl_->rlk);
}

std::size_t SealAdapter::galois_keys_size() const {
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    return serialized_size_of(pimpl_->gk);
}

SerializedBuffer SealAdapter::save_public_key() const {
    if (!pimpl_->has_public) throw std::runtime_error("public key not loaded");
    return serialize_to_buffer(pimpl_->pk);
}

SerializedBuffer SealAdapter::save_secret_key() const {
    if (!pimpl_->has_secret) throw std::runtime_error("secret key not loaded");
    return serialize_to_buffer(pimpl_->sk);
}

SerializedBuffer SealAdapter::save_relin_keys() const {
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    return serialize_to_buffer(pimpl_->rlk);
}

SerializedBuffer SealAdapter::save_galois_keys() const {
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    return serialize_to_buffer(pimpl_->gk);
}

SerializedBuffer SealAdapter::save_cipher(const Cipher& cipher) const {
    return serialize_to_buffer(cipher.pimpl_->ct);
}

void SealAdapter::load_public_key(const SerializedBuffer& buffer) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    load_from_buffer(pimpl_->pk, *pimpl_->context, buffer, "public key");
    pimpl_->encryptor = std::make_unique<seal::Encryptor>(*pimpl_->context, pimpl_->pk);
    pimpl_->has_public = true;
    pimpl_->has_keys = pimpl_->has_secret;
}

void SealAdapter::load_secret_key(const SerializedBuffer& buffer) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    load_from_buffer(pimpl_->sk, *pimpl_->context, buffer, "secret key");
    pimpl_->decryptor = std::make_unique<seal::Decryptor>(*pimpl_->context, pimpl_->sk);
    pimpl_->has_secret = true;
    pimpl_->has_keys = pimpl_->has_public;
}

void SealAdapter::load_relin_keys(const SerializedBuffer& buffer) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    load_from_buffer(pimpl_->rlk, *pimpl_->context, buffer, "relin keys");
    pimpl_->has_relin = true;
}

void SealAdapter::load_galois_keys(const SerializedBuffer& buffer) {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    load_from_buffer(pimpl_->gk, *pimpl_->context, buffer, "galois keys");
    pimpl_->has_galois = true;
}

Cipher SealAdapter::load_cipher(const SerializedBuffer& buffer) const {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    Cipher out;
    load_from_buffer(out.pimpl_->ct, *pimpl_->context, buffer, "ciphertext");
    return out;
}

} // namespace m2424
