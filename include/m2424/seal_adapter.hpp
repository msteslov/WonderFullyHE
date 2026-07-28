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
    /// Создаёт CKKS-контекст и проверяет корректность профиля параметров.
    static SealAdapter create(const CkksProfile&);

    /// Генерирует public/secret и, при необходимости, evaluation keys для всех ротаций.
    void generateKeys(bool needRelin = true, bool needGalois = true);
    /// Генерирует только Galois keys, необходимые для указанных шагов ротации.
    void generateKeys(const std::vector<int>& rotationSteps, bool needRelin = true);
    /// Возвращает физическую ёмкость слотов CKKS-контекста.
    std::size_t slotCount() const;

    /// Кодирует вещественные или комплексные значения на исходном уровне и scale профиля.
    Plain encode(const std::vector<double>&);
    Plain encodeComplex(const std::vector<std::complex<double>>&);
    /// Кодирует данные на уровне и со scale целевого ciphertext для последующей plaintext-операции.
    Plain encodeFor(const std::vector<double>&, const Cipher&);
    Plain encodeComplexFor(const std::vector<std::complex<double>>&, const Cipher&);
    /// Кодирует данные на уровне целевого ciphertext с явно заданным и проверенным scale.
    Plain encodeComplexAtScaleFor(const std::vector<std::complex<double>>&, double scale, const Cipher&);
    Plain encodeScalarFor(double, const Cipher&);
    Plain encodeScalarAtScaleFor(double, double scale, const Cipher&);
    /// Шифрует, расшифровывает и декодирует данные с загруженными ключами.
    Cipher encrypt(const Plain&);
    Plain decrypt(const Cipher&);
    std::vector<double> decode(const Plain&);
    std::vector<std::complex<double>> decodeComplex(const Plain&);

    /// Выполняет базовую CKKS-операцию без скрытой policy.
    Cipher add(const Cipher&, const Cipher&);
    Cipher sub(const Cipher&, const Cipher&);
    Cipher addPlain(const Cipher&, const Plain&);
    Cipher subPlain(const Cipher&, const Plain&);
    Cipher multiplyPlain(const Cipher&, const Plain&);
    Cipher multiply(const Cipher&, const Cipher&);
    Cipher relinearize(const Cipher&);
    Cipher rescaleToNext(const Cipher&);
    /// Снижает ciphertext до уровня target; повышение уровня невозможно.
    Cipher modSwitchTo(const Cipher&, const Cipher&);
    /// Согласует ciphertext с target для сложения, если их CKKS scale отличаются не более чем на один процент.
    Cipher alignForAddition(const Cipher&, const Cipher&);
    /// Выполняет ротацию CKKS-слотов с ранее сгенерированным Galois key.
    Cipher rotate(const Cipher&, int steps);

    /// Возвращает метаданные ciphertext и размеры сериализованных ключей/данных.
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

} // namespace m2424
