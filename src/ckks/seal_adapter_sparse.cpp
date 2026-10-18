#include "seal_adapter_internal.hpp"

#include <seal/util/ntt.h>
#include <seal/util/rlwe.h>
#include <seal/util/uintarithsmallmod.h>
#include <seal/util/uintcore.h>

#include <algorithm>
#include <atomic>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace m2424 {
namespace {

// Identity automorphism dispatches the audited generic switch_key_inplace.
// The key container is separate from ordinary Galois keys; it encrypts the
// SOURCE secret under the DESTINATION secret, not tau(s) under itself.
seal::GaloisKeys identitySwitchKey(const seal::SEALContext& context,
                                  const seal::SecretKey& from,
                                  const seal::SecretKey& to) {
    const auto& parameters = context.key_context_data()->parms();
    const auto& moduli = parameters.coeff_modulus();
    const auto degree = parameters.poly_modulus_degree();
    const auto decompositionSize =
        context.first_context_data()->parms().coeff_modulus().size();
    seal::GaloisKeys keys;
    keys.parms_id() = context.key_parms_id();
    keys.data().resize(1);
    keys.data()[0].resize(decompositionSize);
    for (std::size_t j = 0; j < decompositionSize; ++j) {
        auto& output = keys.data()[0][j].data();
        seal::util::encrypt_zero_symmetric(
            to, context, context.key_parms_id(), true, false, output);
        const auto factor = moduli.back().value() % moduli[j].value();
        for (std::size_t i = 0; i < degree; ++i) {
            const auto value = seal::util::multiply_uint_mod(
                from.data()[j * degree + i], factor, moduli[j]);
            output.data(0)[j * degree + i] = seal::util::add_uint_mod(
                output.data(0)[j * degree + i], value, moduli[j]);
        }
    }
    return keys;
}

} // namespace

void SealAdapter::generateSparseBootstrapKeys(std::size_t h) {
    if (!pimpl_->has_secret || !pimpl_->ordinaryDistributionKnown) {
        throw std::invalid_argument(
            "Sparse keys require an adapter-generated ordinary ternary secret");
    }
    const auto& data = *pimpl_->context->key_context_data();
    const auto& parameters = data.parms();
    const auto degree = parameters.poly_modulus_degree();
    if (h < 2 || h > degree) {
        throw std::invalid_argument("Sparse weight must be in [2,N]; no K=1 shortcut");
    }
    seal::RandomToStandardAdapter rng(parameters.random_generator()->create());
    std::vector<std::size_t> positions(degree);
    std::iota(positions.begin(), positions.end(), 0);
    // Partial Fisher-Yates with unbiased bounded draws; uniform subset, iid signs.
    std::vector<int> coefficients(degree);
    for (std::size_t i = 0; i < h; ++i) {
        auto j = std::uniform_int_distribution<std::size_t>(i, degree - 1)(rng);
        std::swap(positions[i], positions[j]);
        coefficients[positions[i]] =
            std::uniform_int_distribution<int>(0, 1)(rng) ? 1 : -1;
    }
    seal::SecretKey sparse;
    const auto& moduli = parameters.coeff_modulus();
    sparse.data().resize(degree * moduli.size());
    for (std::size_t j = 0; j < moduli.size(); ++j) {
        for (std::size_t i = 0; i < degree; ++i) {
            sparse.data()[j * degree + i] = coefficients[i] < 0
                ? moduli[j].value() - 1 : coefficients[i];
        }
    }
    seal::util::ntt_negacyclic_harvey(
        seal::util::RNSIter(sparse.data().data(), degree),
        moduli.size(), data.small_ntt_tables());
    sparse.parms_id() = data.parms_id();
    auto encapsulation = identitySwitchKey(*pimpl_->context, pimpl_->sk, sparse);
    auto restoration = identitySwitchKey(*pimpl_->context, sparse, pimpl_->sk);
    static std::atomic<std::uint64_t> generation{0};
    pimpl_->sparseSk = std::move(sparse);
    pimpl_->sparseWeight = h;
    pimpl_->encapsulationKey = std::move(encapsulation);
    pimpl_->restorationKey = std::move(restoration);
    pimpl_->sparseGeneration = ++generation;
    seal::util::seal_memzero(coefficients.data(), coefficients.size() * sizeof(int));
}

bool SealAdapter::hasSparseEncapsulationKey() const noexcept {
    return pimpl_ && pimpl_->sparseGeneration && pimpl_->encapsulationKey.has_key(1);
}

bool SealAdapter::hasSparseRestorationKey() const noexcept {
    return pimpl_ && pimpl_->sparseGeneration && pimpl_->restorationKey.has_key(1);
}

SparseKeyMetadata SealAdapter::sparseKeyMetadata() const {
    SparseKeyMetadata metadata;
    metadata.degree = slotCount() * 2;
    metadata.weight = pimpl_->sparseWeight;
    metadata.context = contextFingerprint();
    metadata.generation = pimpl_->sparseGeneration;
    metadata.ordinaryDistributionKnown = pimpl_->ordinaryDistributionKnown;
    metadata.distribution =
        "uniform signed fixed-weight ternary: uniform h-subset, independent uniform {-1,+1} signs";
    auto count = [](const seal::KSwitchKeys& keys) {
        std::size_t result = 0;
        for (const auto& row : keys.data()) result += row.size();
        return result;
    };
    metadata.publicKeySamples = pimpl_->has_public ? 1 : 0;
    metadata.relinSamples = pimpl_->has_relin ? count(pimpl_->rlk) : 0;
    metadata.galoisSamples = pimpl_->has_galois ? count(pimpl_->gk) : 0;
    metadata.encryptionSamples = pimpl_->encryptionSamples;
    metadata.encryptionModuli = pimpl_->encryptionModuli;
    if (pimpl_->has_galois) {
        for (std::size_t i = 0; i < pimpl_->gk.data().size(); ++i) {
            if (!pimpl_->gk.data()[i].empty()) {
                metadata.galoisElements.push_back(static_cast<std::uint32_t>(2 * i + 1));
            }
        }
    }
    return metadata;
}

SparseCipher SealAdapter::encapsulateSparse(const Cipher& input) {
    if (!hasSparseEncapsulationKey()) {
        throw std::invalid_argument("Missing sparse encapsulation key");
    }
    if (info(input).ciphertextSize != 2) {
        throw std::invalid_argument("Encapsulation requires two components");
    }
    Cipher output = input;
    pimpl_->evaluator->apply_galois_inplace(
        output.pimpl_->ct, 1, pimpl_->encapsulationKey);
    return SparseCipher(std::move(output), pimpl_->sparseGeneration);
}

SparseRaisedCipher SealAdapter::modRaiseSparse(const SparseCipher& input) {
    if (!hasSparseEncapsulationKey() || input.generation_ != pimpl_->sparseGeneration) {
        throw std::invalid_argument("Sparse ciphertext key generation mismatch");
    }
    return SparseRaisedCipher(modRaiseToTop(input.cipher_), input.generation_);
}

RaisedCipher SealAdapter::restoreSparseForFirstFactor(const SparseRaisedCipher& input) {
    if (!hasSparseRestorationKey() || input.generation_ != pimpl_->sparseGeneration) {
        throw std::invalid_argument("Missing or mismatched sparse restoration key");
    }
    Cipher output = input.cipher_.cipher_;
    pimpl_->evaluator->apply_galois_inplace(
        output.pimpl_->ct, 1, pimpl_->restorationKey);
#ifdef M2424_ENABLE_SPARSE_DIAGNOSTICS
    ++pimpl_->sparseRestorationOperations;
#endif
    return RaisedCipher(std::move(output), input.cipher_.sourceCoeffModulusSize_);
}

#ifdef M2424_ENABLE_SPARSE_DIAGNOSTICS
RaisedCipher SealAdapter::restoreSparse(const SparseRaisedCipher& input) {
    return restoreSparseForFirstFactor(input);
}
#endif

} // namespace m2424
