#include "m2424/mod1_circuit.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool close(double lhs, double rhs, double tolerance) {
    return std::abs(lhs - rhs) <= tolerance;
}

bool close(m2424::Complex lhs, m2424::Complex rhs, double tolerance) {
    return std::abs(lhs - rhs) <= tolerance;
}

void require_cos_encrypted_matches_plain(const m2424::BootstrapMod1Model& model,
                                         double tolerance) {
    m2424::Mod1Circuit circuit(model);
    require(circuit.encrypted_evaluation_available(), "CosDiscrete encrypted evaluation should be available");

    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_deep_ckks());
    adapter.keygen(true, false);

    const m2424::ComplexVector input{{1.0e-5, 0.0}, {-1.5e-5, 0.5e-5}, {0.75e-5, -0.25e-5}};
    const auto expected = circuit.evaluate_plain(input);
    const auto encrypted = circuit.evaluate(adapter, adapter.encrypt(adapter.encode_complex(input)));
    const auto actual = adapter.decode_complex(adapter.decrypt(encrypted));

    for (std::size_t i = 0; i < input.size(); ++i) {
        if (!close(actual[i], expected[i], tolerance)) {
            throw std::runtime_error("encrypted CosDiscrete output does not match plaintext reference");
        }
    }
}

} // namespace

int main() {
    bool ok = true;
    try {
        m2424::Mod1Circuit legacy({m2424::BootstrapMod1Type::LegacySineP3, 3, 0, 8, 40.0});
        m2424::Mod1Circuit legacy_da({m2424::BootstrapMod1Type::LegacySineP3, 3, 1, 8, 40.0});
        m2424::Mod1Circuit cos3({m2424::BootstrapMod1Type::CosDiscrete, 3, 1, 8, 40.0});
        m2424::Mod1Circuit cos5({m2424::BootstrapMod1Type::CosDiscrete, 5, 0, 8, 40.0});
        m2424::Mod1Circuit cos7({m2424::BootstrapMod1Type::CosDiscrete, 7, 0, 8, 40.0});
        m2424::Mod1Circuit cos30({m2424::BootstrapMod1Type::CosDiscrete, 31, 3, 8, 60.0});

        const double x = 1.0e-4;
        const auto legacy_value = legacy.evaluate_plain(x);
        const auto legacy_da_value = legacy_da.evaluate_plain(x);
        const auto cos_value = cos30.evaluate_plain(x);
        const auto sine_reference = std::sin(2.0 * kPi * x) / (2.0 * kPi);

        const m2424::Complex z{x, -0.5 * x};
        const auto cos_complex = cos30.evaluate_plain(z);
        const auto complex_reference = std::sin(2.0 * kPi * z) / (2.0 * kPi);

        const std::vector<double> plain_input{x, -x, 0.5 * x};
        const auto plain_output = cos30.evaluate_plain(plain_input);
        const m2424::ComplexVector complex_input{{x, 0.0}, {-x, 0.25 * x}};
        const auto complex_output = cos30.evaluate_plain(complex_input);

        bool invalid_model_threw = false;
        try {
            m2424::Mod1Circuit invalid({m2424::BootstrapMod1Type::CosDiscrete, 30, 3, 8, 60.0});
            (void)invalid;
        } catch (const std::invalid_argument&) {
            invalid_model_threw = true;
        }

        const bool conditions =
            legacy.model().type == m2424::BootstrapMod1Type::LegacySineP3
            && legacy.estimated_levels() == 3
            && legacy_da.estimated_levels() == 6
            && legacy.encrypted_evaluation_available()
            && legacy_da.encrypted_evaluation_available()
            && cos3.encrypted_evaluation_available()
            && cos5.encrypted_evaluation_available()
            && cos7.encrypted_evaluation_available()
            && !cos30.encrypted_evaluation_available()
            && cos30.estimated_levels() == 8
            && close(legacy_value, x, 1e-8)
            && close(legacy_da_value, x, 1e-8)
            && close(cos_value, sine_reference, 1e-15)
            && close(cos_complex, complex_reference, 1e-15)
            && plain_output.size() == plain_input.size()
            && complex_output.size() == complex_input.size()
            && invalid_model_threw;

        ok = conditions;
        if (ok) {
            require_cos_encrypted_matches_plain(
                {m2424::BootstrapMod1Type::CosDiscrete, 3, 0, 8, 40.0},
                1e-6);
            require_cos_encrypted_matches_plain(
                {m2424::BootstrapMod1Type::CosDiscrete, 5, 0, 8, 40.0},
                1e-6);
            require_cos_encrypted_matches_plain(
                {m2424::BootstrapMod1Type::CosDiscrete, 7, 0, 8, 40.0},
                1e-6);
        }
    } catch (const std::exception& error) {
        ok = false;
        std::printf("[test_mod1_circuit] FAIL: %s\n", error.what());
    }

    if (ok) {
        std::printf("[test_mod1_circuit] PASS\n");
    }
    return ok ? 0 : 1;
}
