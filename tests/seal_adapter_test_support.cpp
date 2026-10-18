#include "bootstrap_fixture.hpp"
#include "sparse_bootstrap_oracle.hpp"
#include "../src/ckks/seal_adapter_internal.hpp"

#include <seal/util/ntt.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m2424 {
namespace {

bool isPowerOfTwo(std::size_t value) {
    return value != 0 && (value & (value - 1)) == 0;
}

void validateFixtureProfile(const CkksProfile& profile) {
    if (!isPowerOfTwo(profile.polyModulusDegree)) {
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

seal::EncryptionParameters makeFixtureParameters(const CkksProfile& profile) {
    seal::EncryptionParameters parameters(seal::scheme_type::ckks);
    parameters.set_poly_modulus_degree(profile.polyModulusDegree);
    parameters.set_coeff_modulus(
        seal::CoeffModulus::Create(profile.polyModulusDegree, profile.coeffModulusBits));
    return parameters;
}

} // namespace

SealAdapter test::BootstrapFixture::create(const CkksProfile& profile) {
    validateFixtureProfile(profile);
    SealAdapter adapter;
    adapter.pimpl_->profile = profile;
    adapter.pimpl_->context = std::make_shared<seal::SEALContext>(
        makeFixtureParameters(profile), true, seal::sec_level_type::none);
    if (!adapter.pimpl_->context->parameters_set()) {
        throw std::invalid_argument("Invalid explicit test-only CKKS fixture");
    }
    adapter.pimpl_->encoder = std::make_unique<seal::CKKSEncoder>(*adapter.pimpl_->context);
    adapter.pimpl_->evaluator = std::make_unique<seal::Evaluator>(*adapter.pimpl_->context);
    adapter.pimpl_->scale = profile.scale;
    adapter.pimpl_->slotCount = adapter.pimpl_->encoder->slot_count();
    return adapter;
}

Cipher test::BootstrapFixture::input(SealAdapter& adapter, const Plain& plain) {
    auto data = adapter.pimpl_->context->get_context_data(plain.pimpl_->pt.parms_id());
    if (!data || !plain.pimpl_->pt.is_ntt_form()) {
        throw std::invalid_argument("Fixture requires NTT plaintext at source level");
    }
    Cipher output;
    output.pimpl_->ct.resize(*adapter.pimpl_->context, data->parms_id(), 2);
    const auto size = data->parms().poly_modulus_degree()
        * data->parms().coeff_modulus().size();
    std::copy_n(plain.pimpl_->pt.data(), size, output.pimpl_->ct.data(0));
    std::fill_n(output.pimpl_->ct.data(1), size, 1);
    output.pimpl_->ct.is_ntt_form() = true;
    output.pimpl_->ct.scale() = plain.pimpl_->pt.scale();
    return output;
}

std::vector<int> test::SparseBootstrapOracle::secret(const SealAdapter& adapter) {
    const auto& data = *adapter.pimpl_->context->key_context_data();
    const auto degree = adapter.slotCount() * 2;
    std::vector<std::uint64_t> coefficients(
        adapter.pimpl_->sparseSk.data().data(),
        adapter.pimpl_->sparseSk.data().data() + degree);
    seal::util::inverse_ntt_negacyclic_harvey(
        coefficients.data(), data.small_ntt_tables()[0]);
    std::vector<int> output;
    for (auto value : coefficients) {
        output.push_back(value == data.parms().coeff_modulus()[0].value() - 1
            ? -1 : static_cast<int>(value));
    }
    return output;
}

std::vector<double> test::SparseBootstrapOracle::original(
        SealAdapter& adapter, const Cipher& cipher) {
    RaisedCipher raised(Cipher(cipher), adapter.coeffModulusSize(cipher));
    return adapter.decryptRaisedCoefficientsAtRaisedModulus(raised);
}

std::vector<double> test::SparseBootstrapOracle::source(
        SealAdapter& adapter, const SparseCipher& cipher) {
    auto old = std::move(adapter.pimpl_->decryptor);
    adapter.pimpl_->decryptor = std::make_unique<seal::Decryptor>(
        *adapter.pimpl_->context, adapter.pimpl_->sparseSk);
    try {
        auto output = original(adapter, cipher.cipher_);
        adapter.pimpl_->decryptor = std::move(old);
        return output;
    } catch (...) {
        adapter.pimpl_->decryptor = std::move(old);
        throw;
    }
}

std::vector<double> test::SparseBootstrapOracle::raised(
        SealAdapter& adapter, const SparseRaisedCipher& cipher) {
    auto old = std::move(adapter.pimpl_->decryptor);
    adapter.pimpl_->decryptor = std::make_unique<seal::Decryptor>(
        *adapter.pimpl_->context, adapter.pimpl_->sparseSk);
    try {
        auto output = adapter.decryptRaisedCoefficientsAtRaisedModulus(cipher.cipher_);
        adapter.pimpl_->decryptor = std::move(old);
        return output;
    } catch (...) {
        adapter.pimpl_->decryptor = std::move(old);
        throw;
    }
}

std::size_t test::SparseBootstrapOracle::restorationCount(const SealAdapter& adapter) {
    return adapter.pimpl_->sparseRestorationOperations;
}

void test::SparseBootstrapOracle::removeKey(SealAdapter& adapter, bool restoration) {
    if (restoration) {
        adapter.pimpl_->restorationKey = {};
    } else {
        adapter.pimpl_->encapsulationKey = {};
    }
}

} // namespace m2424
