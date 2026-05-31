#pragma once

#include "m2424/bootstrap_plan.hpp"
#include "m2424/eval_mod.hpp"

#include <cstddef>
#include <vector>

namespace m2424 {

class Mod1Circuit {
public:
    explicit Mod1Circuit(BootstrapMod1Model model);

    const BootstrapMod1Model& model() const noexcept;
    std::size_t estimated_levels() const;
    bool encrypted_evaluation_available() const noexcept;

    double evaluate_plain(double input) const;
    Complex evaluate_plain(Complex input) const;
    std::vector<double> evaluate_plain(const std::vector<double>& input) const;
    ComplexVector evaluate_plain(const ComplexVector& input) const;

    Cipher evaluate(SealAdapter& adapter, const Cipher& input) const;

private:
    BootstrapMod1Model model_;
};

} // namespace m2424
