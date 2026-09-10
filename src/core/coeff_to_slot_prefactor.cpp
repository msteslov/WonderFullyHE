#include "m2424/coeff_to_slot_prefactor.hpp"

#include <seal/memorymanager.h>
#include <seal/util/uintarith.h>
#include <seal/util/uintcore.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace m2424 {
namespace {
using Words = std::vector<std::uint64_t>;
struct Dyadic { std::uint64_t mantissa; int exponent; };
Dyadic unpack(std::uint64_t bits) {
    const int exponent = static_cast<int>((bits >> 52) & 2047);
    return {(bits & 0xfffffffffffffULL) | (exponent ? 1ULL << 52 : 0),
            exponent ? exponent - 1023 - 52 : -1074};
}
Words shifted(Words value, int shift) {
    if (shift) {
        value.resize(value.size() + static_cast<std::size_t>(shift / 64) + 1);
        seal::util::left_shift_uint(value.data(), shift, value.size(), value.data());
    }
    return value;
}
int compare(const Words& a, const Words& b) {
    return seal::util::compare_uint(a.data(), a.size(), b.data(), b.size());
}
}
CoeffToSlotPrefactor CoeffToSlotPrefactor::sourceNormalization(const BootstrapInputContext& input) {
    double scale;
    std::memcpy(&scale, &input.scaleBinary64Bits, sizeof(scale));
    if (!std::isfinite(scale) || scale <= 0 || input.sourcePrimes.empty())
        throw std::invalid_argument("normalization requires exact positive source scale and modulus");
    CoeffToSlotPrefactor result;
    result.scaleBits_ = input.scaleBinary64Bits;
    result.primes_ = input.sourcePrimes;
    for (const auto prime : result.primes_) {
        if (prime < 2) throw std::invalid_argument("invalid source modulus factor");
        Words product(result.denominator_.size() + 1);
        seal::util::multiply_uint(result.denominator_.data(), result.denominator_.size(),
                                 prime, product.size(), product.data());
        result.denominator_ = std::move(product);
    }
    return result;
}
bool CoeffToSlotPrefactor::isIdentity() const noexcept {
    return primes_.empty() && scaleBits_ == 0x3ff0000000000000ULL;
}
bool CoeffToSlotPrefactor::operator==(const CoeffToSlotPrefactor& other) const noexcept {
    return scaleBits_ == other.scaleBits_ && primes_ == other.primes_;
}
double CoeffToSlotPrefactor::multiplyRounded(double operand) const {
    static_assert(std::numeric_limits<double>::is_iec559 && sizeof(double) == sizeof(std::uint64_t));
    if (!std::isfinite(operand)) throw std::invalid_argument("nonfinite prefactor operand");
    if (isIdentity() || operand == 0) return operand;
    std::uint64_t bits;
    std::memcpy(&bits, &operand, sizeof(bits));
    const auto a = unpack(scaleBits_), b = unpack(bits);
    Words numerator(2);
    seal::util::multiply_uint(&a.mantissa, 1, b.mantissa, numerator.size(), numerator.data());
    const int binaryExponent = a.exponent + b.exponent;
    int ratioExponent = seal::util::get_significant_bit_count_uint(numerator.data(), numerator.size())
        - seal::util::get_significant_bit_count_uint(denominator_.data(), denominator_.size());
    if (compare(shifted(numerator, std::max(0, -ratioExponent)),
                shifted(denominator_, std::max(0, ratioExponent))) < 0) --ratioExponent;
    const int exponent = binaryExponent + ratioExponent;
    if (exponent > 1023) throw std::overflow_error("prefactored diagonal overflows binary64");
    if (exponent < -1075) return std::copysign(0.0, operand);
    const int unitExponent = std::max(-1074, exponent - 52);
    Words n = shifted(numerator, std::max(0, binaryExponent - unitExponent));
    Words d = shifted(denominator_, std::max(0, unitExponent - binaryExponent));
    const auto count = std::max(n.size(), d.size()) + 1;
    n.resize(count); d.resize(count);
    Words quotient(count);
    auto pool = seal::MemoryManager::GetPool();
    seal::util::divide_uint_inplace(n.data(), d.data(), count, quotient.data(), pool);
    const int rounding = compare(shifted(n, 1), d);
    std::uint64_t significand = quotient[0];
    if (rounding > 0 || (rounding == 0 && (significand & 1))) ++significand;
    const double result = std::ldexp(static_cast<double>(significand), unitExponent);
    if (!std::isfinite(result)) throw std::overflow_error("rounded prefactored diagonal overflows binary64");
    return std::copysign(result, operand);
}
} // namespace m2424
