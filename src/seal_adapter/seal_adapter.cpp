#include "m2424/seal_adapter.hpp"

#include <seal/seal.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

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

SealAdapter SealAdapter::create(const CkksProfile& profile) {
    SealAdapter a;
    a.pimpl_->profile = profile;

    auto parms = make_ckks_parms(profile);
    a.pimpl_->context = std::make_shared<seal::SEALContext>(parms, /*expand_mod_chain*/ true);
    a.pimpl_->encoder = std::make_unique<seal::CKKSEncoder>(*a.pimpl_->context);
    a.pimpl_->evaluator = std::make_unique<seal::Evaluator>(*a.pimpl_->context);
    a.pimpl_->scale = profile.scale;
    a.pimpl_->slot_count = a.pimpl_->encoder->slot_count();
    return a;
}

void SealAdapter::keygen(bool need_relin, bool need_galois) {
    if (!pimpl_->context) throw std::runtime_error("SealAdapter not initialized");
    seal::KeyGenerator keygen(*pimpl_->context);
    pimpl_->sk = keygen.secret_key();
    keygen.create_public_key(pimpl_->pk);
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

Plain SealAdapter::encode(const std::vector<double>& vals) {
    if (!pimpl_->encoder) throw std::runtime_error("CKKSEncoder not initialized");
    if (vals.size() > pimpl_->slot_count || (pimpl_->profile.slots && vals.size() > pimpl_->profile.slots)) {
        throw std::invalid_argument("input size exceeds configured slots");
    }
    Plain out;
    pimpl_->encoder->encode(vals, pimpl_->scale, out.pimpl_->pt);
    return out;
}

Cipher SealAdapter::encrypt(const Plain& plain) {
    if (!pimpl_->has_keys || !pimpl_->encryptor) throw std::runtime_error("keys not generated");
    Cipher out;
    pimpl_->encryptor->encrypt(plain.pimpl_->pt, out.pimpl_->ct);
    return out;
}

Plain SealAdapter::decrypt(const Cipher& cipher) {
    if (!pimpl_->has_keys || !pimpl_->decryptor) throw std::runtime_error("keys not generated");
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

Cipher SealAdapter::mul_relin_rescale(const Cipher& a, const Cipher& b) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    Cipher out;
    pimpl_->evaluator->multiply(a.pimpl_->ct, b.pimpl_->ct, out.pimpl_->ct);
    pimpl_->evaluator->relinearize_inplace(out.pimpl_->ct, pimpl_->rlk);
    pimpl_->evaluator->rescale_to_next_inplace(out.pimpl_->ct);
    return out;
}

Cipher SealAdapter::rotate(const Cipher& c, int steps) {
    if (!pimpl_->evaluator) throw std::runtime_error("Evaluator not initialized");
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    Cipher out;
    pimpl_->evaluator->rotate_vector(c.pimpl_->ct, steps, pimpl_->gk, out.pimpl_->ct);
    return out;
}

} // namespace m2424
