#pragma once

#include "m2424/diagonal_transform.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

enum class BootstrapCircuitOrder {
    ModRaiseFirst,
    SlotsToCoeffsFirst
};

enum class BootstrapTransformBackend {
    DenseDiagonal,
    FftLike
};

enum class BootstrapDftType {
    HomomorphicEncode,
    HomomorphicDecode
};

struct BootstrapDftLayer {
    std::string name;
    DiagonalLinearTransform transform;
    double plain_scale_log2{40.0};
    double scaling_log2{0.0};
    std::size_t level_budget_index{};
};

struct BootstrapDftLayerTrace {
    std::string stage;
    std::size_t layer{};
    std::size_t chain_index{};
    double scale_log2{};
    double coeff_modulus_log2{};
};

struct BootstrapDftPlan {
    std::size_t slots{};
    BootstrapDftType type{BootstrapDftType::HomomorphicDecode};
    std::vector<std::size_t> levels;
    double scaling_log2{};
    double plain_scale_log2{40.0};
    std::vector<BootstrapDftLayer> layers;

    ComplexVector apply_plain(const ComplexVector& input) const;
    std::vector<int> rotation_steps() const;
};

class FactorizedLinearTransform {
public:
    explicit FactorizedLinearTransform(BootstrapDftPlan plan);

    const BootstrapDftPlan& plan() const noexcept;
    ComplexVector apply_plain(const ComplexVector& input) const;
    std::vector<int> rotation_steps() const;
    Cipher apply(SealAdapter& adapter,
                 const Cipher& input,
                 std::vector<BootstrapDftLayerTrace>* trace = nullptr,
                 const std::string& stage_prefix = {}) const;

private:
    BootstrapDftPlan plan_;
};

const char* to_string(BootstrapCircuitOrder order) noexcept;
const char* to_string(BootstrapTransformBackend backend) noexcept;
const char* to_string(BootstrapDftType type) noexcept;

BootstrapDftPlan make_bootstrap_dft_plan(std::size_t slots,
                                         BootstrapDftType type,
                                         double plain_scale_log2,
                                         std::vector<std::size_t> levels = {1},
                                         double scaling_log2 = 0.0);

BootstrapDftPlan make_small_slots4_stc_plan(double plain_scale_log2,
                                            std::vector<std::size_t> levels = {1},
                                            double scaling_log2 = 0.0);

BootstrapDftPlan make_small_slots4_cts_plan(double plain_scale_log2,
                                            std::vector<std::size_t> levels = {1},
                                            double scaling_log2 = 0.0);

} // namespace m2424
