#include "m2424/bootstrap_dft.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

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
    const std::size_t layer_count = 1;
    const double layer_scaling_log2 = scaling_log2 / static_cast<double>(layer_count);
    const double layer_scale = std::exp2(layer_scaling_log2);
    ComplexMatrix matrix = type == BootstrapDftType::HomomorphicDecode
        ? canonical_embedding_matrix(slots)
        : invert_matrix(canonical_embedding_matrix(slots));
    matrix = scaled_matrix(std::move(matrix), layer_scale);

    std::vector<BootstrapDftLayer> layers;
    layers.push_back(BootstrapDftLayer{
        type == BootstrapDftType::HomomorphicDecode ? "coeff_to_slots" : "slots_to_coeffs",
        DiagonalLinearTransform::from_matrix(matrix),
        plain_scale_log2,
        layer_scaling_log2,
        levels.front()
    });
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
