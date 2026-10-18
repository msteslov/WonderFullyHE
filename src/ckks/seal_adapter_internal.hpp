#pragma once

#include "m2424/seal_adapter.hpp"

#include <seal/seal.h>

namespace m2424 {

struct Plain::Impl {
    seal::Plaintext pt;
};

struct Cipher::Impl {
    seal::Ciphertext ct;
};

struct SealAdapter::Impl {
    CkksProfile profile{};

    std::shared_ptr<seal::SEALContext> context;
    std::unique_ptr<seal::CKKSEncoder> encoder;
    std::unique_ptr<seal::Evaluator> evaluator;
    std::unique_ptr<seal::Encryptor> encryptor;
    std::unique_ptr<seal::Decryptor> decryptor;

#ifdef M2424_ENABLE_SPARSE_DIAGNOSTICS
    std::size_t sparseRestorationOperations{};
#endif
    seal::SecretKey sparseSk;
    seal::GaloisKeys encapsulationKey, restorationKey;
    std::size_t sparseWeight{}, encryptionSamples{};
    std::uint64_t sparseGeneration{};
    bool ordinaryDistributionKnown{};
    std::vector<std::pair<std::vector<std::uint64_t>, std::size_t>> encryptionModuli;
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

} // namespace m2424
