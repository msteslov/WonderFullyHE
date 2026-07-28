#pragma once

#include <cstddef>
#include <complex>
#include <cstdint>
#include <memory>
#include <vector>

namespace m2424 {

struct CkksProfile {
    std::size_t polyModulusDegree{};
    std::vector<int> coeffModulusBits{};
    double scale{};
    std::size_t slots{};
};

struct CipherInfo {
    double scale{};
    std::size_t chainIndex{};
    std::size_t coeffModulusSize{};
    std::size_t ciphertextSize{};
    double coeffModulusLog2{};
};

using SerializedBuffer = std::vector<std::uint8_t>;

class Plain {
public:
    Plain();
    ~Plain();
    Plain(const Plain&);
    Plain& operator=(const Plain&);
    Plain(Plain&&) noexcept;
    Plain& operator=(Plain&&) noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;

    friend class SealAdapter;
};

class Cipher {
public:
    Cipher();
    ~Cipher();
    Cipher(const Cipher&);
    Cipher& operator=(const Cipher&);
    Cipher(Cipher&&) noexcept;
    Cipher& operator=(Cipher&&) noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;

    friend class SealAdapter;
};

class SealAdapter {
public:
    /// Creates a CKKS context and validates its parameter profile.
    static SealAdapter create(const CkksProfile&);

    /// Generates public, secret, and optionally evaluation keys for all rotations.
    void generateKeys(bool needRelin = true, bool needGalois = true);
    /// Generates only the Galois keys required by the listed rotation steps.
    void generateKeys(const std::vector<int>& rotationSteps, bool needRelin = true);
    /// Returns the physical CKKS slot capacity of this context.
    std::size_t slotCount() const;

    /// Encodes real or complex values at the profile's initial level and scale.
    Plain encode(const std::vector<double>&);
    Plain encodeComplex(const std::vector<std::complex<double>>&);
    /// Encodes data at the level and scale of target; intended for a subsequent plaintext operation.
    Plain encodeFor(const std::vector<double>&, const Cipher&);
    Plain encodeComplexFor(const std::vector<std::complex<double>>&, const Cipher&);
    /// Encodes data at target's level with an explicitly chosen, validated scale.
    Plain encodeComplexAtScaleFor(const std::vector<std::complex<double>>&, double scale, const Cipher&);
    Plain encodeScalarFor(double, const Cipher&);
    Plain encodeScalarAtScaleFor(double, double scale, const Cipher&);
    /// Encrypts, decrypts, and decodes values using loaded public or secret keys.
    Cipher encrypt(const Plain&);
    Plain decrypt(const Cipher&);
    std::vector<double> decode(const Plain&);
    std::vector<std::complex<double>> decodeComplex(const Plain&);

    /// Performs the corresponding CKKS primitive without changing hidden policy.
    Cipher add(const Cipher&, const Cipher&);
    Cipher sub(const Cipher&, const Cipher&);
    Cipher addPlain(const Cipher&, const Plain&);
    Cipher subPlain(const Cipher&, const Plain&);
    Cipher multiplyPlain(const Cipher&, const Plain&);
    Cipher multiply(const Cipher&, const Cipher&);
    Cipher relinearize(const Cipher&);
    Cipher rescaleToNext(const Cipher&);
    /// Switches ciphertext to the level of target; it never switches levels upward.
    Cipher modSwitchTo(const Cipher&, const Cipher&);
    /// Aligns a ciphertext with target for addition when their CKKS scales differ by at most one percent.
    Cipher alignForAddition(const Cipher&, const Cipher&);
    /// Rotates CKKS slots using a previously generated Galois key.
    Cipher rotate(const Cipher&, int steps);

    /// Returns ciphertext metadata and serialized key/ciphertext sizes.
    std::size_t serializedSize(const Cipher&) const;
    CipherInfo info(const Cipher&) const;
    double scale(const Cipher&) const;
    double coeffModulusLog2(const Cipher&) const;
    std::vector<int> coeffModulusBits() const;
    std::size_t chainIndex(const Cipher&) const;
    std::size_t coeffModulusSize(const Cipher&) const;
    std::size_t publicKeySize() const;
    std::size_t relinKeysSize() const;
    std::size_t galoisKeysSize() const;

    SerializedBuffer savePublicKey() const;
    SerializedBuffer saveSecretKey() const;
    SerializedBuffer saveRelinKeys() const;
    SerializedBuffer saveGaloisKeys() const;
    SerializedBuffer saveCipher(const Cipher&) const;

    void loadPublicKey(const SerializedBuffer&);
    void loadSecretKey(const SerializedBuffer&);
    void loadRelinKeys(const SerializedBuffer&);
    void loadGaloisKeys(const SerializedBuffer&);
    Cipher loadCipher(const SerializedBuffer&) const;

    SealAdapter();
    ~SealAdapter();
    SealAdapter(const SealAdapter&) = delete;
    SealAdapter& operator=(const SealAdapter&) = delete;
    SealAdapter(SealAdapter&&) noexcept;
    SealAdapter& operator=(SealAdapter&&) noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

/// Explicit composition for callers that intentionally consume one CKKS level.
Cipher multiplyPlainAndRescale(SealAdapter&, const Cipher&, const Plain&);
/// Explicit composition for callers that intentionally relinearize and consume one CKKS level.
Cipher multiplyRelinearizeAndRescale(SealAdapter&, const Cipher&, const Cipher&);

} // namespace m2424
