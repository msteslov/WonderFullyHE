#pragma once
#include <gmpxx.h>
#include <string>
#include <stdexcept>
namespace m2424::experimental {
inline mpq_class parseExactDecimal(const std::string& text) {
    if (text.empty()) throw std::invalid_argument("empty exact decimal");
    std::size_t position = 0;
    bool negative = false;
    if (text[position] == '+' || text[position] == '-') {
        negative = text[position] == '-';
        if (++position == text.size()) throw std::invalid_argument("invalid exact decimal");
    }
    const auto exponentPosition = text.find_first_of("eE", position);
    const auto mantissaEnd = exponentPosition == std::string::npos ? text.size() : exponentPosition;
    long exponent = 0;
    if (exponentPosition != std::string::npos) {
        std::size_t consumed = 0;
        try {
            exponent = std::stol(text.substr(exponentPosition + 1), &consumed, 10);
        } catch (const std::exception&) {
            throw std::invalid_argument("invalid exact decimal exponent");
        }
        if (consumed != text.size() - exponentPosition - 1
            || exponent < -100000 || exponent > 100000)
            throw std::invalid_argument("invalid exact decimal exponent");
    }
    std::string digits;
    std::size_t fractionalDigits = 0;
    bool decimalPoint = false;
    for (; position < mantissaEnd; ++position) {
        const char character = text[position];
        if (character == '.') {
            if (decimalPoint) throw std::invalid_argument("invalid exact decimal point");
            decimalPoint = true;
        } else if (character >= '0' && character <= '9') {
            digits.push_back(character);
            if (decimalPoint) ++fractionalDigits;
        } else {
            throw std::invalid_argument("invalid exact decimal digit");
        }
    }
    if (digits.empty()) throw std::invalid_argument("exact decimal has no digits");
    mpz_class numerator(digits, 10);
    if (negative) numerator = -numerator;
    mpz_class denominator = 1;
    if (fractionalDigits > 0)
        mpz_ui_pow_ui(denominator.get_mpz_t(), 10, static_cast<unsigned long>(fractionalDigits));
    if (exponent > 0) {
        mpz_class power;
        mpz_ui_pow_ui(power.get_mpz_t(), 10, static_cast<unsigned long>(exponent));
        numerator *= power;
    } else if (exponent < 0) {
        mpz_class power;
        mpz_ui_pow_ui(power.get_mpz_t(), 10, static_cast<unsigned long>(-exponent));
        denominator *= power;
    }
    mpq_class result(numerator, denominator);
    result.canonicalize();
    return result;
}
}
