#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/m2424.hpp"
#include "bootstrap_fixture.hpp"

#include <gmpxx.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

template <class Function>
double elapsedMs(Function&& function) {
    const auto start = Clock::now();
    function();
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

m2424::BootstrapBound bound(double value, const char* provenance) {
    return {value, m2424::BootstrapBoundKind::Deterministic, provenance, {}};
}

const char* statusName(m2424::BootstrapCertificationStatus status) {
    using Status = m2424::BootstrapCertificationStatus;
    switch (status) {
    case Status::Certified: return "Certified";
    case Status::InvalidInput: return "InvalidInput";
    case Status::MissingExactModulusContext: return "MissingExactModulusContext";
    case Status::InputScaleMismatch: return "InputScaleMismatch";
    case Status::DomainViolation: return "DomainViolation";
    case Status::SparseSecretCertificateUnavailable: return "SparseSecretCertificateUnavailable";
    case Status::KeySwitchBoundUnavailable: return "KeySwitchBoundUnavailable";
    case Status::ExtractionNotCertified: return "ExtractionNotCertified";
    case Status::CleaningBoundUnavailable: return "CleaningBoundUnavailable";
    case Status::SlotToCoeffGainUnavailable: return "SlotToCoeffGainUnavailable";
    case Status::ScaleScheduleInfeasible: return "ScaleScheduleInfeasible";
    case Status::HeadroomViolation: return "HeadroomViolation";
    case Status::InsufficientLevels: return "InsufficientLevels";
    case Status::MissingEvaluationKeys: return "MissingEvaluationKeys";
    case Status::ErrorBudgetExceeded: return "ErrorBudgetExceeded";
    case Status::FailureProbabilityExceeded: return "FailureProbabilityExceeded";
    case Status::SecurityBudgetExceeded: return "SecurityBudgetExceeded";
    case Status::RequiredBoundUnavailable: return "RequiredBoundUnavailable";
    case Status::ResourceBudgetExceeded: return "ResourceBudgetExceeded";
    case Status::LiftBoundUnavailable: return "LiftBoundUnavailable";
    case Status::UnsupportedEvalRoundDomain: return "UnsupportedEvalRoundDomain";
    case Status::TestOnlyAssumption: return "TestOnlyAssumption";
    }
    return "Unknown";
}

double maximumError(const m2424::ComplexVector& actual,
                    const m2424::ComplexVector& expected) {
    if (actual.size() != expected.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double maximum = 0.0;
    for (std::size_t index = 0; index < actual.size(); ++index) {
        maximum = std::max(maximum, std::abs(actual[index] - expected[index]));
    }
    return maximum;
}

} // namespace

int main() {
    try {
        constexpr std::size_t degree = 16;
        constexpr std::size_t slots = degree / 2;
        const m2424::CkksProfile profile{
            degree, std::vector<int>(15, 50), std::ldexp(1.0, 49), slots
        };
        m2424::BootstrapRequest request;
        request.target.targetAbsoluteError = 1e-6;
        request.upstream.raisedMagnitude = bound(
            degree * (std::ldexp(1.0, 31) + 1.0) / profile.scale,
            "Synthetic fixture raised magnitude bound");
        request.upstream.messageMagnitude = bound(
            std::ldexp(1.0, 31), "Exact test message coefficient bound");
        request.upstream.sourceNoiseMagnitude = bound(
            1.0, "Synthetic c=(m,1) source noise bound");
        request.lift = {
            1,
            m2424::BootstrapLiftEvidence::TestFixtureAssumption,
            "Explicit test-only K=1 fixture assumption"
        };

        m2424::SealAdapter adapter;
        double context_ms = elapsedMs([&] {
            adapter = m2424::test::BootstrapFixture::create(profile);
        });
        m2424::Bootstrapper bootstrap(degree, 2, 2);
        auto keys = bootstrap.rotationKeys();
        keys.push_back(0);
        double keys_ms = elapsedMs([&] { adapter.generateKeys(keys, true); });

        m2424::Cipher input;
        std::vector<double> coefficients(degree);
        double input_ms = elapsedMs([&] {
            const auto top = adapter.encrypt(adapter.encode({0.0}));
            const auto bottom = adapter.modSwitchToChainIndex(top, 0);
            for (std::size_t index = 0; index < degree; ++index) {
                coefficients[index] = static_cast<double>(static_cast<int>(index % 5) - 2)
                    * std::ldexp(1.0, 30);
            }
            std::vector<std::uint64_t> primes = adapter.dataModulusValues();
            primes.push_back(adapter.specialKeyModulusValue());
            std::vector<std::uint64_t> residues;
            residues.reserve(primes.size() * coefficients.size());
            for (auto prime : primes) {
                for (double coefficient : coefficients) {
                    mpz_class value(static_cast<long>(coefficient));
                    mpz_class modulus(std::to_string(prime));
                    mpz_class residue;
                    mpz_mod(residue.get_mpz_t(), value.get_mpz_t(), modulus.get_mpz_t());
                    residues.push_back(std::stoull(residue.get_str()));
                }
            }
            auto encoded = adapter.encodePolynomialRnsAtKeyScale(residues, profile.scale);
            auto source_plain = adapter.modSwitchPlainTo(encoded, bottom);
            input = m2424::test::BootstrapFixture::input(adapter, source_plain);
        });

        m2424::BootstrapPlan plan;
        double plan_ms = elapsedMs([&] { plan = bootstrap.prepare(adapter, input, request); });
        m2424::BootstrapResult result;
        double apply_ms = elapsedMs([&] { result = bootstrap.apply(adapter, input, plan); });

        double observed_error = std::numeric_limits<double>::infinity();
        double decode_ms = 0.0;
        if (result.output) {
            decode_ms = elapsedMs([&] {
                auto message = coefficients;
                for (double& coefficient : message) {
                    coefficient /= profile.scale;
                }
                const auto expected = m2424::coeffToSlotReference(message);
                const auto actual = adapter.decodeComplex(adapter.decrypt(*result.output));
                observed_error = maximumError(actual, expected);
            });
        }

        const bool ready = plan.executionReadiness().status
            == m2424::BootstrapCertificationStatus::Certified;
        const bool test_only = result.trace.result.status
            == m2424::BootstrapCertificationStatus::TestOnlyAssumption;
        const bool passed = ready && test_only && result.output
            && observed_error <= request.target.targetAbsoluteError;

        std::printf("scenario,poly_modulus_degree,slots,context_ms,keys_ms,input_ms,"
                    "plan_ms,apply_ms,decode_ms,total_ms,nodes,levels_used,"
                    "output_chain_index,observed_error,target_error,execution_status,"
                    "global_status,test_only,status\n");
        std::printf("test_fixture,%zu,%zu,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
                    "%zu,%zu,%zu,%.6e,%.6e,%s,%s,1,%s\n",
                    degree,
                    slots,
                    context_ms,
                    keys_ms,
                    input_ms,
                    plan_ms,
                    apply_ms,
                    decode_ms,
                    context_ms + keys_ms + input_ms + plan_ms + apply_ms + decode_ms,
                    result.trace.nodes.size(),
                    result.trace.plan.levelsUsed,
                    result.output ? adapter.chainIndex(*result.output) : 0,
                    observed_error,
                    request.target.targetAbsoluteError,
                    statusName(plan.executionReadiness().status),
                    statusName(result.trace.result.status),
                    passed ? "PASS" : "FAIL");
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 2;
    }
}
