#include "m2424/approximation_lab.hpp"
#include "m2424/evalmod_cost_model.hpp"
#include "m2424/exact_modular_oracle.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>

int main() {
    const auto domain = m2424::analyzeEvalModDomain({4096, 1.0, 0.25, 0.1, 0.0005, -128.0});
    const auto oracle = m2424::exactCoefficientOracle({34, 23}, {101, 103}, 97,
                                                       m2424::OracleFloat(16, 256));
    // x itself is exact only on the central interval; this intentionally exercises K=0 certificate semantics.
    const auto central = m2424::analyzeEvalModDomain({1, 0.0, 0.0, 0.2, 0.0, -128.0});
    const auto certificate = m2424::certifyEvalModPolynomial({0.0, 1.0}, central, 33, 1e-3, 1e-12);
    const auto cost = m2424::estimateEvalModCost({31, 5, 8, 8, 5}, {2.0, 0.5, 0.25, 1024},
                                                 50, 4);

    bool discontinuityRejected = false;
    try {
        (void)m2424::analyzeEvalModDomain({1, 0.0, 0.0, 0.5, 0.0, -128.0});
    } catch (const std::domain_error&) {
        discontinuityRejected = true;
    }

    const bool ok = domain.integerBound >= 1 && std::abs(domain.normalizedResidualBound - 0.1005) < 1e-12
        && domain.discontinuityMargin > 0.0 && oracle.reconstructed == 5791
        && oracle.centeredSourceCoefficient == -29
        && abs(oracle.expectedValue - m2424::OracleFloat("-1.8125", 256)) < m2424::OracleFloat("1e-70", 256)
        && certificate.approximationMaxError < 1e-15 && certificate.derivativeMax == 1.0
        && certificate.complexNeighborhoodError < 1e-15
        && cost.latencyMs == 21.25 && cost.peakWorkingSetBytes == 4096
        && cost.requiredModulusBits == 300 && discontinuityRejected;
    std::printf("[test_evalmod_contract] %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
