#pragma once

#include <cstddef>
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
};

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
    std::size_t slot_count() const;

    Plain encode(const std::vector<double>&);
    Cipher encrypt(const Plain&);
    Plain decrypt(const Cipher&);
    std::vector<double> decode(const Plain&);

    Cipher add(const Cipher&, const Cipher&);
    Cipher sub(const Cipher&, const Cipher&);
    Cipher mul_relin_rescale(const Cipher&, const Cipher&);
    Cipher rotate(const Cipher&, int steps);

    std::size_t serialized_size(const Cipher&) const;
    CipherInfo info(const Cipher&) const;
    double scale(const Cipher&) const;
    std::size_t chain_index(const Cipher&) const;
    std::size_t coeff_modulus_size(const Cipher&) const;
    std::size_t public_key_size() const;
    std::size_t relin_keys_size() const;
    std::size_t galois_keys_size() const;

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
