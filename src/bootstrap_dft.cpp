#include "m2424/bootstrap_dft.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace m2424 {
namespace {

bool is_supported_bootstrap_slots(std::size_t slots) {
    return slots == 4 || slots == 8 || slots == 16;
}

ComplexMatrix scaled_matrix(ComplexMatrix matrix, double scale) {
    for (auto& row : matrix) {
        for (auto& value : row) {
            value *= scale;
        }
    }
    return matrix;
}

ComplexMatrix identity_matrix(std::size_t n) {
    ComplexMatrix matrix(n, ComplexVector(n, Complex{0.0, 0.0}));
    for (std::size_t i = 0; i < n; ++i) {
        matrix[i][i] = Complex{1.0, 0.0};
    }
    return matrix;
}

bool close_complex(Complex lhs, Complex rhs) {
    return std::abs(lhs - rhs) <= 1e-9;
}

std::vector<ComplexMatrix> factor_eval_layers(const ComplexVector& roots) {
    const std::size_t n = roots.size();
    if (n == 0 || (n & (n - 1)) != 0) {
        throw std::invalid_argument("bootstrap DFT roots must be a non-zero power of two");
    }
    if (n == 1) {
        return {};
    }

    const std::size_t half = n / 2;
    ComplexVector squared_roots;
    squared_roots.reserve(half);
    std::vector<std::size_t> group_for_row(n);
    for (std::size_t row = 0; row < n; ++row) {
        const Complex squared = roots[row] * roots[row];
        auto it = std::find_if(squared_roots.begin(), squared_roots.end(), [&](Complex value) {
            return close_complex(value, squared);
        });
        if (it == squared_roots.end()) {
            if (squared_roots.size() >= half) {
                throw std::invalid_argument("bootstrap DFT roots do not form square-root pairs");
            }
            group_for_row[row] = squared_roots.size();
            squared_roots.push_back(squared);
        } else {
            group_for_row[row] = static_cast<std::size_t>(std::distance(squared_roots.begin(), it));
        }
    }
    if (squared_roots.size() != half) {
        throw std::invalid_argument("bootstrap DFT roots must collapse into half as many squared roots");
    }

    ComplexMatrix split(n, ComplexVector(n, Complex{0.0, 0.0}));
    for (std::size_t i = 0; i < half; ++i) {
        split[i][2 * i] = Complex{1.0, 0.0};
        split[half + i][2 * i + 1] = Complex{1.0, 0.0};
    }

    std::vector<ComplexMatrix> layers;
    layers.push_back(std::move(split));
    for (const auto& child : factor_eval_layers(squared_roots)) {
        ComplexMatrix block(n, ComplexVector(n, Complex{0.0, 0.0}));
        for (std::size_t row = 0; row < half; ++row) {
            for (std::size_t col = 0; col < half; ++col) {
                block[row][col] = child[row][col];
                block[half + row][half + col] = child[row][col];
            }
        }
        layers.push_back(std::move(block));
    }

    ComplexMatrix combine(n, ComplexVector(n, Complex{0.0, 0.0}));
    for (std::size_t row = 0; row < n; ++row) {
        const std::size_t group = group_for_row[row];
        combine[row][group] = Complex{1.0, 0.0};
        combine[row][half + group] = roots[row];
    }
    layers.push_back(std::move(combine));
    return layers;
}

ComplexVector canonical_roots(std::size_t slots) {
    constexpr double kPi = 3.141592653589793238462643383279502884;
    const std::uint64_t polynomial_degree = static_cast<std::uint64_t>(2 * slots);
    const std::uint64_t cyclotomic_order = 2 * polynomial_degree;
    const Complex zeta = std::exp(Complex{0.0, 2.0 * kPi / static_cast<double>(cyclotomic_order)});
    ComplexVector roots;
    roots.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        std::uint64_t exponent = 1;
        for (std::size_t j = 0; j < i; ++j) {
            exponent = (exponent * 5U) % cyclotomic_order;
        }
        roots.push_back(std::pow(zeta, static_cast<double>(exponent)));
    }
    return roots;
}

void validate_plan_inputs(std::size_t slots,
                          double plain_scale_log2,
                          const std::vector<std::size_t>& levels,
                          double scaling_log2) {
    if (!is_supported_bootstrap_slots(slots)) {
        throw std::invalid_argument("bootstrap DFT v1 supports only slots=4,8,16");
    }
    if (!std::isfinite(plain_scale_log2) || plain_scale_log2 <= 0.0) {
        throw std::invalid_argument("bootstrap DFT plaintext scale log2 must be positive and finite");
    }
    if (!std::isfinite(scaling_log2)) {
        throw std::invalid_argument("bootstrap DFT scaling log2 must be finite");
    }
    if (levels.empty()) {
        throw std::invalid_argument("bootstrap DFT levels must not be empty");
    }
    for (std::size_t level : levels) {
        if (level == 0) {
            throw std::invalid_argument("bootstrap DFT levels must be positive");
        }
    }
}

std::vector<BootstrapDftLayer> make_layers(std::size_t slots,
                                           BootstrapDftType type,
                                           double plain_scale_log2,
                                           const std::vector<std::size_t>& levels,
                                           double scaling_log2) {
    std::vector<ComplexMatrix> matrices;
    if (type == BootstrapDftType::HomomorphicDecode) {
        matrices = factor_eval_layers(canonical_roots(slots));
    } else {
        matrices.push_back(invert_matrix(canonical_embedding_matrix(slots)));
    }
    if (matrices.empty()) {
        matrices.push_back(identity_matrix(slots));
    }
    const std::size_t layer_count = matrices.size();
    const double layer_scaling_log2 = scaling_log2 / static_cast<double>(layer_count);
    const double layer_scale = std::exp2(layer_scaling_log2);

    std::vector<BootstrapDftLayer> layers;
    layers.reserve(matrices.size());
    for (std::size_t i = 0; i < matrices.size(); ++i) {
        const auto level_budget_index = levels[std::min(i, levels.size() - 1)];
        layers.push_back(BootstrapDftLayer{
            type == BootstrapDftType::HomomorphicDecode ? "coeff_to_slots_fft_layer" : "slots_to_coeffs",
            DiagonalLinearTransform::from_matrix(scaled_matrix(std::move(matrices[i]), layer_scale)),
            plain_scale_log2,
            layer_scaling_log2,
            level_budget_index
        });
    }
    return layers;
}

} // namespace

const char* to_string(BootstrapCircuitOrder order) noexcept {
    switch (order) {
    case BootstrapCircuitOrder::ModRaiseFirst:
        return "ModRaiseFirst";
    case BootstrapCircuitOrder::SlotsToCoeffsFirst:
        return "SlotsToCoeffsFirst";
    }
    return "unknown";
}

const char* to_string(BootstrapTransformBackend backend) noexcept {
    switch (backend) {
    case BootstrapTransformBackend::DenseDiagonal:
        return "DenseDiagonal";
    case BootstrapTransformBackend::FftLike:
        return "FftLike";
    }
    return "unknown";
}

const char* to_string(BootstrapDftType type) noexcept {
    switch (type) {
    case BootstrapDftType::HomomorphicEncode:
        return "HomomorphicEncode";
    case BootstrapDftType::HomomorphicDecode:
        return "HomomorphicDecode";
    }
    return "unknown";
}

ComplexVector BootstrapDftPlan::apply_plain(const ComplexVector& input) const {
    if (input.size() != slots) {
        throw std::invalid_argument("bootstrap DFT input size must match plan slots");
    }
    auto current = input;
    for (const auto& layer : layers) {
        current = layer.transform.apply_plain(current);
    }
    return current;
}

std::vector<int> BootstrapDftPlan::rotation_steps() const {
    std::vector<int> steps;
    for (const auto& layer : layers) {
        auto layer_steps = layer.transform.rotation_steps();
        steps.insert(steps.end(), layer_steps.begin(), layer_steps.end());
    }
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

FactorizedLinearTransform::FactorizedLinearTransform(BootstrapDftPlan plan)
    : plan_(std::move(plan)) {
    if (plan_.layers.empty()) {
        throw std::invalid_argument("factorized linear transform must have at least one layer");
    }
}

const BootstrapDftPlan& FactorizedLinearTransform::plan() const noexcept {
    return plan_;
}

ComplexVector FactorizedLinearTransform::apply_plain(const ComplexVector& input) const {
    return plan_.apply_plain(input);
}

std::vector<int> FactorizedLinearTransform::rotation_steps() const {
    return plan_.rotation_steps();
}

Cipher FactorizedLinearTransform::apply(SealAdapter& adapter,
                                        const Cipher& input,
                                        std::vector<BootstrapDftLayerTrace>* trace,
                                        const std::string& stage_prefix) const {
    auto current = input;
    for (std::size_t i = 0; i < plan_.layers.size(); ++i) {
        const auto& layer = plan_.layers[i];
        current = layer.transform.apply_at_plain_scale(adapter, current, std::exp2(layer.plain_scale_log2));
        if (trace != nullptr) {
            const auto info = adapter.info(current);
            trace->push_back(BootstrapDftLayerTrace{
                stage_prefix.empty() ? layer.name : stage_prefix,
                i,
                info.chain_index,
                std::log2(info.scale),
                info.coeff_modulus_log2
            });
        }
    }
    return current;
}

BootstrapDftPlan make_bootstrap_dft_plan(std::size_t slots,
                                         BootstrapDftType type,
                                         double plain_scale_log2,
                                         std::vector<std::size_t> levels,
                                         double scaling_log2) {
    validate_plan_inputs(slots, plain_scale_log2, levels, scaling_log2);
    BootstrapDftPlan plan;
    plan.slots = slots;
    plan.type = type;
    plan.levels = std::move(levels);
    plan.scaling_log2 = scaling_log2;
    plan.plain_scale_log2 = plain_scale_log2;
    plan.layers = make_layers(slots, type, plain_scale_log2, plan.levels, scaling_log2);
    return plan;
}

} // namespace m2424
