#include "m2424/bootstrap.hpp"
#include "m2424/bootstrap_prototype.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

namespace {

template <typename Fn>
double elapsed_ms(Fn&& fn) {
    const auto started = std::chrono::steady_clock::now();
    fn();
    const auto finished = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(finished - started).count();
}

m2424::ComplexVector make_expected(std::size_t slots, double amplitude) {
    m2424::ComplexVector expected;
    expected.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        expected.push_back({amplitude * std::sin(static_cast<double>(i) / 4.0), 0.0});
    }
    return expected;
}

std::vector<double> real_part(const m2424::ComplexVector& values) {
    std::vector<double> result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(value.real());
    }
    return result;
}

double max_error(m2424::SealAdapter& adapter,
                 const m2424::Cipher& cipher,
                 const m2424::ComplexVector& expected) {
    auto actual = adapter.decode_complex(adapter.decrypt(cipher));
    actual.resize(expected.size());
    double result = 0.0;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        result = std::max(result, std::abs(actual[i] - expected[i]));
    }
    return result;
}

double least_squares_real_gain(m2424::SealAdapter& adapter,
                               const m2424::Cipher& cipher,
                               const m2424::ComplexVector& expected) {
    auto actual = adapter.decode_complex(adapter.decrypt(cipher));
    actual.resize(expected.size());
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        numerator += actual[i].real() * expected[i].real();
        denominator += expected[i].real() * expected[i].real();
    }
    return denominator == 0.0 ? 0.0 : numerator / denominator;
}

std::size_t arg_size_or(int argc, char** argv, int index, std::size_t fallback) {
    if (argc <= index) {
        return fallback;
    }
    const auto parsed = std::strtoull(argv[index], nullptr, 10);
    return parsed == 0 ? fallback : static_cast<std::size_t>(parsed);
}

int arg_int_or(int argc, char** argv, int index, int fallback) {
    if (argc <= index) {
        return fallback;
    }
    const auto parsed = std::strtol(argv[index], nullptr, 10);
    return parsed <= 0 ? fallback : static_cast<int>(parsed);
}

double arg_double_or(int argc, char** argv, int index, double fallback) {
    if (argc <= index) {
        return fallback;
    }
    char* end = nullptr;
    const double parsed = std::strtod(argv[index], &end);
    return end == argv[index] ? fallback : parsed;
}

m2424::EvalModDegree arg_evalmod_or(int argc, char** argv, int index, m2424::EvalModDegree fallback) {
    if (argc <= index) {
        return fallback;
    }
    const std::string value = argv[index];
    if (value == "p3" || value == "P3") {
        return m2424::EvalModDegree::P3;
    }
    if (value == "da" || value == "double_angle" || value == "P3DoubleAngle") {
        return m2424::EvalModDegree::P3DoubleAngle;
    }
    return fallback;
}

std::vector<int> repeated_refresh_moduli(int work_bits, std::size_t work_levels) {
    std::vector<int> bits;
    bits.reserve(work_levels + 2);
    bits.push_back(60);
    for (std::size_t i = 0; i < work_levels; ++i) {
        bits.push_back(work_bits);
    }
    bits.push_back(60);
    return bits;
}

void sanitize(std::string& text) {
    std::replace(text.begin(), text.end(), ',', ';');
}

} // namespace

int main(int argc, char** argv) {
    const std::size_t slots = arg_size_or(argc, argv, 1, 4);
    const std::size_t cycles = arg_size_or(argc, argv, 2, 3);
    const int work_bits = arg_int_or(argc, argv, 3, 40);
    const double correction_factor = arg_double_or(argc, argv, 4, 1.0);
    const auto evalmod_degree = arg_evalmod_or(argc, argv, 5, m2424::EvalModDegree::P3);
    constexpr double tolerance = 1e-3;
    constexpr double amplitude = 1e-5;

    std::vector<int> rotations;
    const double rotation_plan_ms = elapsed_ms([&] {
        rotations = m2424::Bootstrapper::scalable_refresh_rotation_steps(slots);
    });

    const auto profile = work_bits == 40
        ? m2424::profiles::boot_deep_ckks()
        : m2424::CkksProfile{
            32768,
            repeated_refresh_moduli(work_bits, 13),
            std::exp2(static_cast<double>(work_bits)),
            16384
        };
    auto adapter = m2424::SealAdapter::create(profile);
    double keygen_ms = 0.0;
    keygen_ms = elapsed_ms([&] {
        adapter.keygen(rotations, true);
    });

    const auto expected = make_expected(slots, amplitude);
    auto current = adapter.encrypt(adapter.encode(real_part(expected)));
    const auto initial_info = adapter.info(current);

    std::printf("cycle,slots,evalmod,stage,chain_before,chain_after,scale_before_log2,scale_after_log2,continuation_levels,restore_level,preserve_value,max_abs_error,real_gain,cycle_ms,exception,status\n");
    std::printf("0,%zu,%s,setup,%zu,%zu,%.6f,%.6f,%zu,true,true,%.6e,%.6f,%.6f,,PASS\n",
                slots,
                m2424::to_string(evalmod_degree),
                initial_info.chain_index,
                initial_info.chain_index,
                std::log2(initial_info.scale),
                std::log2(initial_info.scale),
                initial_info.chain_index,
                max_error(adapter, current, expected),
                least_squares_real_gain(adapter, current, expected),
                rotation_plan_ms + keygen_ms);

    bool all_pass = true;
    double max_seen_error = 0.0;
    std::size_t min_continuation_levels = initial_info.chain_index;
    for (std::size_t cycle = 1; cycle <= cycles; ++cycle) {
        const auto before = adapter.info(current);
        m2424::BootstrapPrototypeReport report;
        std::string exception;
        double cycle_ms = 0.0;
        bool pass = false;
        try {
            cycle_ms = elapsed_ms([&] {
                m2424::BootstrapPrototype prototype(adapter, slots, tolerance);
                prototype.set_transform_backend(m2424::BootstrapTransformBackend::FftLike);
                prototype.set_circuit_order(m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst);
                prototype.set_evalmod_degree(evalmod_degree);
                prototype.set_plain_scale_log2(std::log2(adapter.info(current).scale));
                prototype.set_output_correction_factor(correction_factor);
                report = prototype.refresh_cipher_checked(current, expected);
            });
            current = report.result;
            const auto after = adapter.info(current);
            const double error = max_error(adapter, current, expected);
            max_seen_error = std::max(max_seen_error, error);
            min_continuation_levels = std::min(min_continuation_levels, report.continuation_levels);
            pass = report.preserve_value_criterion && error <= tolerance;
            all_pass = all_pass && pass;
            std::printf("%zu,%zu,%s,refresh,%zu,%zu,%.6f,%.6f,%zu,%s,%s,%.6e,%.6f,%.6f,,%s\n",
                        cycle,
                        slots,
                        m2424::to_string(evalmod_degree),
                        before.chain_index,
                        after.chain_index,
                        std::log2(before.scale),
                        std::log2(after.scale),
                        report.continuation_levels,
                        report.restore_level_criterion ? "true" : "false",
                        report.preserve_value_criterion ? "true" : "false",
                        error,
                        least_squares_real_gain(adapter, current, expected),
                        cycle_ms,
                        pass ? "PASS" : "FAIL");
        } catch (const std::exception& error) {
            exception = error.what();
            sanitize(exception);
            all_pass = false;
            std::printf("%zu,%zu,%s,refresh,%zu,0,%.6f,0,0,false,false,0,0,%.6f,%s,FAIL\n",
                        cycle,
                        slots,
                        m2424::to_string(evalmod_degree),
                        before.chain_index,
                        std::log2(before.scale),
                        cycle_ms,
                        exception.c_str());
            break;
        }
    }

    std::printf("summary,%zu,%s,multi_cycle,0,0,0,0,%zu,false,%s,%.6e,0,0,,%s\n",
                slots,
                m2424::to_string(evalmod_degree),
                min_continuation_levels,
                all_pass ? "true" : "false",
                max_seen_error,
                all_pass ? "PASS" : "FAIL");
    return all_pass ? 0 : 1;
}
