#pragma once

#include "m2424/bootstrap_plan.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/mod1_approximation.hpp"

#include <cstddef>
#include <vector>

namespace m2424 {

struct Mod1EncryptedEvaluation {
    Cipher result;
    PolynomialEvaluationStrategy strategy{PolynomialEvaluationStrategy::DirectOddPowers};
    std::size_t input_chain_index{};
    std::size_t output_chain_index{};
    std::size_t consumed_levels{};
    double output_scale_log2{};
};

class Mod1Circuit {
public:
    explicit Mod1Circuit(BootstrapMod1Model model);

    const BootstrapMod1Model& model() const noexcept;
    const Mod1Approximation& approximation() const noexcept;
    std::size_t estimated_levels() const;
    bool encrypted_evaluation_available() const noexcept;

    double evaluate_plain(double input) const;
    Complex evaluate_plain(Complex input) const;
    std::vector<double> evaluate_plain(const std::vector<double>& input) const;
    ComplexVector evaluate_plain(const ComplexVector& input) const;

    Cipher evaluate(SealAdapter& adapter, const Cipher& input) const;
    Mod1EncryptedEvaluation evaluate_with_report(SealAdapter& adapter, const Cipher& input) const;

private:
    BootstrapMod1Model model_;
    Mod1Approximation approximation_;
};

} // namespace m2424
