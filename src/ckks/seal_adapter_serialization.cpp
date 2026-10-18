#include "seal_adapter_internal.hpp"

#include <sstream>
#include <stdexcept>
#include <string>

namespace m2424 {
namespace {

template <class T>
std::size_t serializedSizeOf(const T& value) {
    std::ostringstream output(std::ios::binary);
    value.save(output);
    return static_cast<std::size_t>(output.tellp());
}

template <class T>
SerializedBuffer serializeToBuffer(const T& value) {
    std::ostringstream output(std::ios::binary);
    value.save(output);
    const auto bytes = output.str();
    return SerializedBuffer(bytes.begin(), bytes.end());
}

void requireNonEmptyBuffer(const SerializedBuffer& buffer, const char* name) {
    if (buffer.empty()) {
        throw std::invalid_argument(std::string(name) + " buffer must not be empty");
    }
}

template <class T>
void loadFromBuffer(T& value, const seal::SEALContext& context,
                    const SerializedBuffer& buffer, const char* name) {
    requireNonEmptyBuffer(buffer, name);
    value.load(context, reinterpret_cast<const seal::seal_byte*>(buffer.data()), buffer.size());
}

} // namespace

std::size_t SealAdapter::serializedSize(const Cipher& cipher) const {
    return serializedSizeOf(cipher.pimpl_->ct);
}

std::size_t SealAdapter::serializedSize(const Plain& plain) const {
    return serializedSizeOf(plain.pimpl_->pt);
}

std::size_t SealAdapter::publicKeySize() const {
    if (!pimpl_->has_public) throw std::runtime_error("public key not loaded");
    return serializedSizeOf(pimpl_->pk);
}

std::size_t SealAdapter::relinKeysSize() const {
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    return serializedSizeOf(pimpl_->rlk);
}

std::size_t SealAdapter::galoisKeysSize() const {
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    return serializedSizeOf(pimpl_->gk);
}

SerializedBuffer SealAdapter::savePublicKey() const {
    if (!pimpl_->has_public) throw std::runtime_error("public key not loaded");
    return serializeToBuffer(pimpl_->pk);
}

SerializedBuffer SealAdapter::saveSecretKey() const {
    if (!pimpl_->has_secret) throw std::runtime_error("secret key not loaded");
    return serializeToBuffer(pimpl_->sk);
}

SerializedBuffer SealAdapter::saveRelinKeys() const {
    if (!pimpl_->has_relin) throw std::runtime_error("relin keys not generated");
    return serializeToBuffer(pimpl_->rlk);
}

SerializedBuffer SealAdapter::saveGaloisKeys() const {
    if (!pimpl_->has_galois) throw std::runtime_error("galois keys not generated");
    return serializeToBuffer(pimpl_->gk);
}

SerializedBuffer SealAdapter::saveCipher(const Cipher& cipher) const {
    return serializeToBuffer(cipher.pimpl_->ct);
}

void SealAdapter::loadPublicKey(const SerializedBuffer& buffer) {
    pimpl_->ordinaryDistributionKnown = false;
    pimpl_->encapsulationKey = {};
    pimpl_->restorationKey = {};
    pimpl_->sparseGeneration = 0;
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    loadFromBuffer(pimpl_->pk, *pimpl_->context, buffer, "public key");
    pimpl_->encryptor = std::make_unique<seal::Encryptor>(*pimpl_->context, pimpl_->pk);
    pimpl_->has_public = true;
}

void SealAdapter::loadSecretKey(const SerializedBuffer& buffer) {
    pimpl_->encapsulationKey = {};
    pimpl_->restorationKey = {};
    pimpl_->sparseGeneration = 0;
    pimpl_->ordinaryDistributionKnown = false;
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    loadFromBuffer(pimpl_->sk, *pimpl_->context, buffer, "secret key");
    pimpl_->decryptor = std::make_unique<seal::Decryptor>(*pimpl_->context, pimpl_->sk);
    pimpl_->has_secret = true;
}

void SealAdapter::loadRelinKeys(const SerializedBuffer& buffer) {
    pimpl_->ordinaryDistributionKnown = false;
    pimpl_->encapsulationKey = {};
    pimpl_->restorationKey = {};
    pimpl_->sparseGeneration = 0;
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    loadFromBuffer(pimpl_->rlk, *pimpl_->context, buffer, "relin keys");
    pimpl_->has_relin = true;
}

void SealAdapter::loadGaloisKeys(const SerializedBuffer& buffer) {
    pimpl_->ordinaryDistributionKnown = false;
    pimpl_->encapsulationKey = {};
    pimpl_->restorationKey = {};
    pimpl_->sparseGeneration = 0;
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    loadFromBuffer(pimpl_->gk, *pimpl_->context, buffer, "galois keys");
    pimpl_->has_galois = true;
}

Cipher SealAdapter::loadCipher(const SerializedBuffer& buffer) const {
    if (!pimpl_->context) throw std::runtime_error("SEALContext not initialized");
    Cipher output;
    loadFromBuffer(output.pimpl_->ct, *pimpl_->context, buffer, "ciphertext");
    return output;
}

} // namespace m2424
