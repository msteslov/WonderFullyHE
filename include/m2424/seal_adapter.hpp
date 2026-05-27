#pragma once

#include <cstddef>
#include <complex>
#include <cstdint>
#include <memory>
#include <vector>

namespace m2424 {

struct CkksProfile {
    std::size_t poly_modulus_degree{};
    std::vector<int> coeff_modulus_bits{};
    double scale{};
    std::size_t slots{};
};

struct CipherInfo {
    double scale{};
    std::size_t chain_index{};
    std::size_t coeff_modulus_size{};
    std::size_t ciphertext_size{};
    double coeff_modulus_log2{};
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
    static SealAdapter create(const CkksProfile&);

    void keygen(bool need_relin = true, bool need_galois = true);
    void keygen(const std::vector<int>& rotation_steps, bool need_relin = true);
    std::size_t slot_count() const;

    Plain encode(const std::vector<double>&);
    Plain encode_complex(const std::vector<std::complex<double>>&);
    Plain encode_like(const std::vector<double>&, const Cipher&);
    Plain encode_complex_like(const std::vector<std::complex<double>>&, const Cipher&);
    Plain encode_scalar_like(double, const Cipher&);
    Plain encode_scalar_at_scale_like(double, double scale, const Cipher&);
    Cipher encrypt(const Plain&);
    Plain decrypt(const Cipher&);
    std::vector<double> decode(const Plain&);
    std::vector<std::complex<double>> decode_complex(const Plain&);

    Cipher add(const Cipher&, const Cipher&);
    Cipher sub(const Cipher&, const Cipher&);
    Cipher add_plain(const Cipher&, const Plain&);
    Cipher sub_plain(const Cipher&, const Plain&);
    Cipher mul_plain(const Cipher&, const Plain&);
    Cipher mul_plain_rescale(const Cipher&, const Plain&);
    Cipher mul_relin_rescale(const Cipher&, const Cipher&);
    Cipher mod_raise_to_first(const Cipher&);
    Cipher multiply_decoded_value(const Cipher&, double multiplier);
    Cipher mod_switch_to(const Cipher&, const Cipher&);
    Cipher match_level_and_scale(const Cipher&, const Cipher&);
    Cipher rotate(const Cipher&, int steps);

    std::size_t serialized_size(const Cipher&) const;
    CipherInfo info(const Cipher&) const;
    double scale(const Cipher&) const;
    double coeff_modulus_log2(const Cipher&) const;
    double bootstrap_period_log2(const Cipher&) const;
    double bootstrap_period(const Cipher&) const;
    std::vector<int> coeff_modulus_bits() const;
    std::size_t chain_index(const Cipher&) const;
    std::size_t coeff_modulus_size(const Cipher&) const;
    std::size_t public_key_size() const;
    std::size_t relin_keys_size() const;
    std::size_t galois_keys_size() const;

    SerializedBuffer save_public_key() const;
    SerializedBuffer save_secret_key() const;
    SerializedBuffer save_relin_keys() const;
    SerializedBuffer save_galois_keys() const;
    SerializedBuffer save_cipher(const Cipher&) const;

    void load_public_key(const SerializedBuffer&);
    void load_secret_key(const SerializedBuffer&);
    void load_relin_keys(const SerializedBuffer&);
    void load_galois_keys(const SerializedBuffer&);
    Cipher load_cipher(const SerializedBuffer&) const;

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
