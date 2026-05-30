#include "m2424/accuracy.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <numeric>
#include <string>
#include <vector>

namespace {

std::vector<double> make_input(std::size_t n) {
    std::vector<double> values;
    values.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = static_cast<double>(i);
        values.push_back(0.22 + 0.17 * std::sin(x / 5.0) + 0.06 * std::cos(x / 11.0));
    }
    return values;
}

std::vector<double> head(const std::vector<double>& values, std::size_t n) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), n));
}

std::vector<int> chain_bits(int work_bits, std::size_t work_levels) {
    std::vector<int> bits;
    bits.reserve(work_levels + 2);
    bits.push_back(60);
    for (std::size_t i = 0; i < work_levels; ++i) {
        bits.push_back(work_bits);
    }
    bits.push_back(60);
    return bits;
}

double mean(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
}

double max_value(const std::vector<double>& values) {
    return values.empty() ? 0.0 : *std::max_element(values.begin(), values.end());
}

void run_case(const char* experiment,
              const char* variant,
              std::size_t poly_degree,
              int work_bits,
              int scale_log2,
              std::size_t work_levels,
              std::size_t depth,
              std::size_t trials) {
    const std::size_t payload_size = 64;
    std::vector<double> max_errors;
    std::vector<double> mean_errors;

    for (std::size_t trial = 0; trial < trials; ++trial) {
        try {
            const m2424::CkksProfile profile{
                poly_degree,
                chain_bits(work_bits, work_levels),
                std::pow(2.0, static_cast<double>(scale_log2)),
                poly_degree / 2
            };
            auto adapter = m2424::SealAdapter::create(profile);
            adapter.keygen(true, false);

            auto reference = make_input(payload_size);
            auto current = adapter.encrypt(adapter.encode(reference));

            for (std::size_t step = 0; step < depth; ++step) {
                current = adapter.mul_relin_rescale(current, current);
                for (double& value : reference) {
                    value *= value;
                }
            }

            const auto decoded = head(adapter.decode(adapter.decrypt(current)), payload_size);
            const auto accuracy = m2424::compare(reference, decoded, 0.0);
            max_errors.push_back(accuracy.max_abs_error);
            mean_errors.push_back(accuracy.mean_abs_error);
        } catch (const std::exception& error) {
            std::printf("%s,%s,%zu,%d,%d,%zu,%zu,%zu,failed,%s,0,0,0,0\n",
                        experiment,
                        variant,
                        poly_degree,
                        work_bits,
                        scale_log2,
                        work_levels + 2,
                        work_levels,
                        depth,
                        error.what());
            return;
        }
    }

    std::printf("%s,%s,%zu,%d,%d,%zu,%zu,%zu,ok,,%.6e,%.6e,%.6e,%.6e\n",
                experiment,
                variant,
                poly_degree,
                work_bits,
                scale_log2,
                work_levels + 2,
                work_levels,
                depth,
                mean(max_errors),
                max_value(max_errors),
                mean(mean_errors),
                max_value(mean_errors));
}

} // namespace

int main() {
    const std::size_t trials = 3;

    std::printf("experiment,variant,poly_modulus_degree,work_bits,scale_log2,chain_length,work_levels,depth,status,error,mean_max_abs_error,worst_max_abs_error,mean_mean_abs_error,worst_mean_abs_error\n");

    for (std::size_t work_levels = 1; work_levels <= 7; ++work_levels) {
        run_case("chain_length_fixed_depth", "depth1", 16384, 40, 40, work_levels, 1, trials);
        run_case("chain_length_fixed_depth", "depth2", 16384, 40, 40, work_levels, 2, trials);
        run_case("chain_length_fixed_depth", "depth3", 16384, 40, 40, work_levels, 3, trials);
    }

    for (int scale_log2 : {30, 35, 40, 45, 50}) {
        run_case("scale_fixed_chain", "depth2", 16384, 50, scale_log2, 4, 2, trials);
    }

    for (int work_bits : {30, 35, 40, 45, 50}) {
        run_case("work_bits_matched_scale", "depth2", 16384, work_bits, work_bits, 4, 2, trials);
    }

    return 0;
}
