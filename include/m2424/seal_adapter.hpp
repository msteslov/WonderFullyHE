#pragma once

#include <cstddef>
#include <complex>
#include <cstdint>
#include <memory>
#include <vector>

namespace m2424 {

class CoeffToSlot;

/// Параметры CKKS-контекста, используемые при создании SealAdapter.
struct CkksProfile {
    /// Степень полиномиального модуля; должна быть степенью двойки.
    std::size_t polyModulusDegree{};
    /// Битовые размеры простых модулей в modulus chain.
    std::vector<int> coeffModulusBits{};
    /// Исходный scale для кодирования данных.
    double scale{};
    /// Логическое число используемых слотов; 0 означает всю доступную ёмкость.
    std::size_t slots{};
};

/// Снимок параметров ciphertext на текущем уровне modulus chain.
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

/// Ciphertext, успешно прошедший centered RNS ModRaise до верхнего уровня chain.
class RaisedCipher {
public:
    RaisedCipher(RaisedCipher&&) noexcept;
    RaisedCipher& operator=(RaisedCipher&&) noexcept;
    RaisedCipher(const RaisedCipher&) = delete;
    RaisedCipher& operator=(const RaisedCipher&) = delete;

private:
    explicit RaisedCipher(Cipher&& cipher, std::size_t sourceCoeffModulusSize);
    Cipher cipher_;
    std::size_t sourceCoeffModulusSize_{};

    friend class SealAdapter;
    friend class ::m2424::CoeffToSlot;
};

class SealAdapter {
public:
    /**
     * Создаёт CKKS-контекст и проверяет корректность профиля параметров.
     * @return Инициализированный адаптер без ключей.
     * @throws std::invalid_argument Если профиль не поддерживается Microsoft SEAL.
     */
    static SealAdapter create(const CkksProfile&);

    /// Генерирует public/secret и, при необходимости, evaluation keys для всех ротаций.
    /// @param needRelin Создавать ли relinearization keys.
    /// @param needGalois Создавать ли Galois keys для всех поддерживаемых ротаций.
    void generateKeys(bool needRelin = true, bool needGalois = true);
    /// Генерирует только Galois keys, необходимые для указанных шагов ротации.
    /// @param rotationSteps Разрешённые шаги ротации слотов.
    /// @param needRelin Создавать ли relinearization keys.
    void generateKeys(const std::vector<int>& rotationSteps, bool needRelin = true);
    /// Возвращает физическую ёмкость слотов CKKS-контекста.
    std::size_t slotCount() const;

    /// Кодирует вещественные значения на исходном уровне и scale профиля.
    /// @throws std::invalid_argument Если значения не конечны или превышают число слотов.
    Plain encode(const std::vector<double>&);
    /// Кодирует комплексные значения на исходном уровне и scale профиля.
    Plain encodeComplex(const std::vector<std::complex<double>>&);
    /// Кодирует данные на уровне и со scale целевого ciphertext для последующей plaintext-операции.
    Plain encodeFor(const std::vector<double>&, const Cipher&);
    /// Кодирует комплексные данные на уровне и со scale целевого ciphertext.
    Plain encodeComplexFor(const std::vector<std::complex<double>>&, const Cipher&);
    /// Кодирует данные на уровне целевого ciphertext с явно заданным и проверенным scale.
    Plain encodeComplexAtScaleFor(const std::vector<std::complex<double>>&, double scale, const Cipher&);
    /// Кодирует скаляр на уровне и со scale целевого ciphertext.
    Plain encodeScalarFor(double, const Cipher&);
    /// Кодирует скаляр на уровне target с явно заданным scale.
    Plain encodeScalarAtScaleFor(double, double scale, const Cipher&);
    /// Шифрует plaintext public key; ключ должен быть сгенерирован или загружен заранее.
    Cipher encrypt(const Plain&);
    /// Расшифровывает ciphertext secret key; ключ должен быть сгенерирован или загружен заранее.
    Plain decrypt(const Cipher&);
    /// Raw coefficients m + qI: inverse NTT и centered CRT по поднятой базе Q.
    std::vector<double> decryptRaisedCoefficientsAtRaisedModulus(const RaisedCipher&);
    /// Raw coefficients m mod q: inverse NTT и centered CRT по исходной базе q.
    std::vector<double> decryptRaisedCoefficientsAtSourceModulus(const RaisedCipher&);
    /// Декодирует вещественные CKKS-слоты plaintext.
    std::vector<double> decode(const Plain&);
    /// Декодирует комплексные CKKS-слоты plaintext.
    std::vector<std::complex<double>> decodeComplex(const Plain&);

    /// Складывает два ciphertext одного уровня и совместимого scale.
    Cipher add(const Cipher&, const Cipher&);
    /// Вычитает второй ciphertext из первого.
    Cipher sub(const Cipher&, const Cipher&);
    /// Складывает ciphertext и plaintext на совместимых уровне и scale.
    Cipher addPlain(const Cipher&, const Plain&);
    /// Вычитает plaintext из ciphertext на совместимых уровне и scale.
    Cipher subPlain(const Cipher&, const Plain&);
    /// Умножает ciphertext на plaintext без rescale.
    Cipher multiplyPlain(const Cipher&, const Plain&);
    /// Умножает два ciphertext без relinearization и rescale.
    Cipher multiply(const Cipher&, const Cipher&);
    /// Сводит ciphertext после умножения к двум компонентам.
    Cipher relinearize(const Cipher&);
    /// Снижает уровень modulus chain и нормализует scale после умножения.
    Cipher rescaleToNext(const Cipher&);
    /// Снижает ciphertext до уровня target; повышение уровня невозможно.
    Cipher modSwitchTo(const Cipher&, const Cipher&);
    /// Поднимает ciphertext с нижнего уровня текущей modulus chain на верхний уровень через centered RNS lift.
    /// Операция сохраняет residues исходной базы и не использует secret key.
    RaisedCipher modRaiseToTop(const Cipher&);
    /// Выполняет ротацию CKKS-слотов с ранее сгенерированным Galois key.
    /// @param steps Число позиций ротации; для него должен существовать Galois key.
    Cipher rotate(const Cipher&, int steps);
    /// Применяет CKKS automorphism комплексного сопряжения; требуется Galois key шага 0.
    Cipher conjugate(const Cipher&);

    /// Возвращает метаданные ciphertext и размеры сериализованных ключей/данных.
    std::size_t serializedSize(const Cipher&) const;
    /// Возвращает размер сериализованного plaintext.
    std::size_t serializedSize(const Plain&) const;
    /// Возвращает уровень, scale, размер и доступный modulus chain ciphertext.
    CipherInfo info(const Cipher&) const;
    /// Возвращает метаданные ciphertext после ModRaise.
    CipherInfo info(const RaisedCipher&) const;
    /// Возвращает scale ciphertext.
    double scale(const Cipher&) const;
    /// Возвращает суммарный логарифм modulus chain ciphertext по основанию два.
    double coeffModulusLog2(const Cipher&) const;
    /// Возвращает битовые размеры modulus chain из исходного профиля.
    std::vector<int> coeffModulusBits() const;
    /// Возвращает индекс ciphertext в modulus chain.
    std::size_t chainIndex(const Cipher&) const;
    /// Возвращает число активных модулей ciphertext.
    std::size_t coeffModulusSize(const Cipher&) const;
    /// Возвращает размер сериализованного public key.
    std::size_t publicKeySize() const;
    /// Возвращает размер сериализованных relinearization keys.
    std::size_t relinKeysSize() const;
    /// Возвращает размер сериализованных Galois keys.
    std::size_t galoisKeysSize() const;
    /// Проверяет наличие загруженных relinearization keys.
    bool hasRelinKeys() const noexcept;
    /// Проверяет наличие Galois keys для всех указанных шагов ротации.
    bool hasRotationKeys(const std::vector<int>& rotationSteps) const;
    /// Проверяет наличие Galois key для automorphism комплексного сопряжения.
    bool hasConjugationKey() const;

    /// Сериализует public key для передачи вычислительному контексту.
    SerializedBuffer savePublicKey() const;
    /// Сериализует secret key; вызывающий код отвечает за безопасное хранение результата.
    SerializedBuffer saveSecretKey() const;
    /// Сериализует relinearization keys.
    SerializedBuffer saveRelinKeys() const;
    /// Сериализует Galois keys.
    SerializedBuffer saveGaloisKeys() const;
    /// Сериализует ciphertext.
    SerializedBuffer saveCipher(const Cipher&) const;

    /// Загружает public key, сериализованный в контексте с теми же CKKS-параметрами.
    void loadPublicKey(const SerializedBuffer&);
    /// Загружает secret key, сериализованный в контексте с теми же CKKS-параметрами.
    void loadSecretKey(const SerializedBuffer&);
    /// Загружает relinearization keys.
    void loadRelinKeys(const SerializedBuffer&);
    /// Загружает Galois keys.
    void loadGaloisKeys(const SerializedBuffer&);
    /// Загружает ciphertext и проверяет его совместимость с контекстом.
    Cipher loadCipher(const SerializedBuffer&) const;

    SealAdapter();
    ~SealAdapter();
    SealAdapter(const SealAdapter&) = delete;
    SealAdapter& operator=(const SealAdapter&) = delete;
    SealAdapter(SealAdapter&&) noexcept;
    SealAdapter& operator=(SealAdapter&&) noexcept;

private:
    std::vector<double> decryptRaisedCoefficients(const RaisedCipher&, std::size_t modulusCount);

    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // namespace m2424
