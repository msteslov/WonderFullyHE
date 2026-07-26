#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include "m2424/abft.hpp"
#include "m2424/accuracy.hpp"
#include "m2424/bootstrap.hpp"
#include "m2424/bootstrap_scaling.hpp"
#include "m2424/checked_evaluator.hpp"
#include "m2424/diagonal_transform.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/linear_transform.hpp"
#include "m2424/polynomial.hpp"
#include "m2424/profiles.hpp"
#include "m2424/security_report.hpp"
#include "m2424/seal_adapter.hpp"

static std::vector<double> head(const std::vector<double>& values, std::size_t n) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), n));
}

static bool close_enough(const std::vector<double>& expected, const std::vector<double>& actual, double threshold) {
    return m2424::compare(expected, actual, threshold).ok;
}

static bool expect_runtime_error(void (*fn)()) {
    try {
        fn();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

static bool expect_invalid_argument(void (*fn)()) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

template <class Fn>
static bool expect_runtime_error_from(Fn&& fn) {
    try {
        fn();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

static void invalid_profile_case() {
    (void)m2424::SealAdapter::create(m2424::CkksProfile{0, {60, 40, 60}, std::pow(2.0, 40), 0});
}

static void empty_encode_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    (void)adapter.encode({});
}

static void encrypt_without_keys_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    auto p = adapter.encode({1.0, 2.0});
    (void)adapter.encrypt(p);
}

static void mul_without_relin_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    adapter.keygen(false, false);
    auto ct = adapter.encrypt(adapter.encode({1.0, 2.0}));
    (void)adapter.mul_relin_rescale(ct, ct);
}

static void rotate_without_galois_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    adapter.keygen(false, false);
    auto ct = adapter.encrypt(adapter.encode({1.0, 2.0}));
    (void)adapter.rotate(ct, 1);
}

static void accuracy_size_mismatch_case() {
    (void)m2424::compare({1.0, 2.0}, {1.0}, 1e-5);
}

static void unknown_profile_name_case() {
    (void)m2424::profiles::by_name("unknown_ckks_profile");
}

static bool has_stage(const m2424::BootstrapReport& report, const char* name) {
    return std::any_of(report.stages.begin(), report.stages.end(), [name](const m2424::BootstrapStage& stage) {
        return stage.name == name;
    });
}

static bool all_named_profiles_create_contexts() {
    for (const auto& entry : m2424::profiles::all()) {
        auto adapter = m2424::SealAdapter::create(entry.second);
        if (adapter.slot_count() != entry.second.poly_modulus_degree / 2) {
            return false;
        }
        const auto resolved = m2424::profiles::by_name(entry.first);
        if (resolved.poly_modulus_degree != entry.second.poly_modulus_degree
            || resolved.coeff_modulus_bits != entry.second.coeff_modulus_bits
            || resolved.scale != entry.second.scale
            || resolved.slots != entry.second.slots) {
            return false;
        }
    }
    return true;
}

static bool diagonal_transform_plan_ok() {
    const auto matrix = m2424::canonical_embedding_matrix(4);
    const auto inverse = m2424::invert_matrix(matrix);
    const auto transform = m2424::DiagonalLinearTransform::from_matrix(matrix);
    const auto inverse_transform = m2424::DiagonalLinearTransform::from_matrix(inverse);

    const m2424::ComplexVector input{
        {0.1, 0.01},
        {0.2, -0.02},
        {-0.05, 0.03},
        {0.07, -0.04}
    };
    const auto slots = transform.apply_plain(input);
    const auto roundtrip = inverse_transform.apply_plain(slots);
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (std::abs(input[i] - roundtrip[i]) > 1e-10) {
            return false;
        }
    }
    return transform.terms().size() == 4
        && inverse_transform.terms().size() == 4
        && transform.rotation_steps().size() == 6;
}

static bool eval_mod_polynomial_ok() {
    m2424::EvalModPolynomial eval_mod;
    const std::vector<double> input{
        -m2424::EvalModPolynomial::approximation_bound,
        -0.5 * m2424::EvalModPolynomial::approximation_bound,
        0.0,
        0.5 * m2424::EvalModPolynomial::approximation_bound,
        m2424::EvalModPolynomial::approximation_bound
    };
    const auto polynomial = eval_mod.evaluate_plain(input);
    const auto sine = eval_mod.sine_reference(input);
    return m2424::compare(sine, polynomial, 1e-12).ok;
}

int main() {
    const auto prof = m2424::profiles::basic_ckks();
    const auto slots = prof.slots;
    auto adapter = m2424::SealAdapter::create(prof);
    adapter.keygen(true, true);

    const std::size_t N = 64;
    std::vector<double> input; input.reserve(N);
    for (std::size_t i = 0; i < N; ++i) input.push_back(std::sin(static_cast<double>(i) / 10.0));

    auto p = adapter.encode(input);
    auto ct = adapter.encrypt(p);
    const auto public_key_bytes = adapter.save_public_key();
    const auto secret_key_bytes = adapter.save_secret_key();
    const auto relin_key_bytes = adapter.save_relin_keys();
    const auto galois_key_bytes = adapter.save_galois_keys();
    const auto cipher_bytes = adapter.save_cipher(ct);

    auto loaded_adapter = m2424::SealAdapter::create(prof);
    loaded_adapter.load_public_key(public_key_bytes);
    loaded_adapter.load_secret_key(secret_key_bytes);
    loaded_adapter.load_relin_keys(relin_key_bytes);
    loaded_adapter.load_galois_keys(galois_key_bytes);
    auto loaded_ct = loaded_adapter.load_cipher(cipher_bytes);
    auto loaded_out = head(loaded_adapter.decode(loaded_adapter.decrypt(loaded_ct)), N);

    auto public_only_adapter = m2424::SealAdapter::create(prof);
    public_only_adapter.load_public_key(public_key_bytes);
    public_only_adapter.load_galois_keys(galois_key_bytes);
    auto public_only_ct = public_only_adapter.load_cipher(cipher_bytes);
    const bool public_only_cannot_decrypt = expect_runtime_error_from([&] {
        (void)public_only_adapter.decrypt(public_only_ct);
    });
    auto ct2 = adapter.mul_relin_rescale(ct, ct);
    auto p2 = adapter.decrypt(ct2);
    auto out = adapter.decode(p2);

    std::vector<double> ref; ref.reserve(N);
    for (double x : input) ref.push_back(x * x);
    std::vector<double> out_head(out.begin(), out.begin() + std::min(out.size(), N));

    auto mul_accuracy = m2424::compare(ref, out_head, 1e-3);

    auto ct_add = adapter.add(ct, ct);
    std::vector<double> add_ref; add_ref.reserve(N);
    for (double x : input) add_ref.push_back(x + x);
    auto add_out = head(adapter.decode(adapter.decrypt(ct_add)), N);
    m2424::CheckedEvaluator checked(adapter, N, 1e-5);
    auto checked_add = checked.add(ct, ct, add_ref);

    auto ct_sub = adapter.sub(ct_add, ct);
    auto sub_out = head(adapter.decode(adapter.decrypt(ct_sub)), N);

    auto plain_shift = adapter.encode_like(std::vector<double>(N, 0.25), ct);
    auto ct_add_plain = adapter.add_plain(ct, plain_shift);
    auto add_plain_out = head(adapter.decode(adapter.decrypt(ct_add_plain)), N);
    std::vector<double> add_plain_ref; add_plain_ref.reserve(N);
    for (double x : input) add_plain_ref.push_back(x + 0.25);

    auto ct_sub_plain = adapter.sub_plain(ct_add_plain, plain_shift);
    auto sub_plain_out = head(adapter.decode(adapter.decrypt(ct_sub_plain)), N);

    auto scalar = adapter.encode_scalar_like(1.5, ct);
    auto ct_mul_plain = adapter.mul_plain_rescale(ct, scalar);
    auto mul_plain_out = head(adapter.decode(adapter.decrypt(ct_mul_plain)), N);
    std::vector<double> mul_plain_ref; mul_plain_ref.reserve(N);
    for (double x : input) mul_plain_ref.push_back(1.5 * x);

    auto ct_rot = adapter.rotate(ct, 1);
    auto rot_out = head(adapter.decode(adapter.decrypt(ct_rot)), N);
    std::vector<double> rot_ref; rot_ref.reserve(N);
    for (std::size_t i = 0; i < N; ++i) rot_ref.push_back(input[(i + 1) % N]);

    auto abft_payload = m2424::abft::append_checksum(input);
    auto abft_ct = adapter.add(adapter.encrypt(adapter.encode(abft_payload)), adapter.encrypt(adapter.encode(abft_payload)));
    auto abft_out = adapter.decode(adapter.decrypt(abft_ct));
    auto abft_check = m2424::abft::verify_appended_checksum(abft_out, N, 1e-5);

    const auto depth_profile = m2424::profiles::depth_ckks();
    auto depth_adapter = m2424::SealAdapter::create(depth_profile);
    depth_adapter.keygen(true, true);
    std::vector<double> depth_input;
    depth_input.reserve(16);
    for (std::size_t i = 0; i < 16; ++i) depth_input.push_back(0.25 + 0.25 * std::sin(static_cast<double>(i) / 7.0));
    m2424::Bootstrapper bootstrapper(depth_adapter);
    auto bootstrap_report = bootstrapper.analyze_depth(depth_input, 8);

    auto depth_ct = depth_adapter.encrypt(depth_adapter.encode(depth_input));
    m2424::LinearTransform transform({
        {0, {0.5}},
        {1, {0.25}}
    });
    m2424::LinearTransform zero_transform({
        {0, {0.0}},
        {1, {0.0, 0.0}}
    });
    const bool zero_transform_ok = zero_transform.rotation_steps().empty();
    auto linear_ct = transform.apply(depth_adapter, depth_ct);
    auto linear_out = head(depth_adapter.decode(depth_adapter.decrypt(linear_ct)), depth_input.size());
    std::vector<double> linear_ref; linear_ref.reserve(depth_input.size());
    for (std::size_t i = 0; i < depth_input.size(); ++i) {
        const double rotated = i + 1 < depth_input.size() ? depth_input[i + 1] : 0.0;
        linear_ref.push_back(0.5 * depth_input[i] + 0.25 * rotated);
    }

    auto sum_ct = m2424::sum_slots(depth_adapter, depth_ct, depth_input.size());
    auto sum_out = head(depth_adapter.decode(depth_adapter.decrypt(sum_ct)), depth_input.size());
    const double sum_ref_value = std::accumulate(depth_input.begin(), depth_input.end(), 0.0);

    m2424::PolynomialEvaluator polynomial({
        {1, 0.75},
        {2, 0.0},
        {3, -0.125}
    });
    auto polynomial_ct = polynomial.evaluate(depth_adapter, depth_ct);
    auto polynomial_out = head(depth_adapter.decode(depth_adapter.decrypt(polynomial_ct)), depth_input.size());
    std::vector<double> polynomial_ref; polynomial_ref.reserve(depth_input.size());
    for (double x : depth_input) polynomial_ref.push_back(0.75 * x - 0.125 * x * x * x);

    const auto basic_security = m2424::analyze_security("basic_ckks", prof);
    const auto balanced_security = m2424::analyze_security("balanced_ckks", m2424::profiles::balanced_ckks());
    const auto depth_security = m2424::analyze_security("depth_ckks", depth_profile);
    const auto high_precision_security = m2424::analyze_security("high_precision_ckks", m2424::profiles::high_precision_ckks());
    const auto boot_security = m2424::analyze_security("boot_ckks", m2424::profiles::boot_ckks());
    const auto boot_deep_security = m2424::analyze_security("boot_deep_ckks", m2424::profiles::boot_deep_ckks());
    const auto precision_boot_deep_security = m2424::analyze_security(
        "precision_boot_deep_ckks", m2424::profiles::precision_boot_deep_ckks());
    const auto precision_boot_ultra_security = m2424::analyze_security(
        "precision_boot_ultra_ckks_59", m2424::profiles::precision_boot_ultra_ckks_59());
    const std::vector<m2424::SecurityReport> security_reports{
        basic_security,
        balanced_security,
        depth_security,
        high_precision_security,
        boot_security,
        boot_deep_security,
        precision_boot_deep_security,
        precision_boot_ultra_security
    };
    const bool security_report_ok = basic_security.total_coeff_modulus_bits == 200
        && basic_security.tc128_limit == 218
        && basic_security.passes_tc128
        && !basic_security.passes_tc192
        && basic_security.effective_level == m2424::SecurityLevel::TC128
        && balanced_security.total_coeff_modulus_bits == 218
        && balanced_security.passes_tc128
        && !balanced_security.passes_tc192
        && balanced_security.effective_level == m2424::SecurityLevel::TC128
        && depth_security.total_coeff_modulus_bits == 280
        && depth_security.tc128_limit == 438
        && depth_security.tc192_limit == 305
        && depth_security.passes_tc128
        && depth_security.passes_tc192
        && !depth_security.passes_tc256
        && depth_security.effective_level == m2424::SecurityLevel::TC192
        && high_precision_security.total_coeff_modulus_bits == 270
        && high_precision_security.passes_tc128
        && high_precision_security.passes_tc192
        && !high_precision_security.passes_tc256
        && high_precision_security.effective_level == m2424::SecurityLevel::TC192
        && boot_security.total_coeff_modulus_bits == 400
        && boot_security.passes_tc128
        && !boot_security.passes_tc192
        && boot_security.effective_level == m2424::SecurityLevel::TC128
        && boot_deep_security.total_coeff_modulus_bits == 880
        && boot_deep_security.tc128_limit == 881
        && boot_deep_security.passes_tc128
        && !boot_deep_security.passes_tc192
        && boot_deep_security.effective_level == m2424::SecurityLevel::TC128
        && precision_boot_deep_security.total_coeff_modulus_bits == 870
        && precision_boot_deep_security.tc128_limit == 881
        && precision_boot_deep_security.passes_tc128
        && !precision_boot_deep_security.passes_tc192
        && precision_boot_deep_security.effective_level == m2424::SecurityLevel::TC128
        && precision_boot_ultra_security.total_coeff_modulus_bits == 840
        && precision_boot_ultra_security.tc128_limit == 881
        && precision_boot_ultra_security.passes_tc128
        && !precision_boot_ultra_security.passes_tc192
        && precision_boot_ultra_security.effective_level == m2424::SecurityLevel::TC128
        && m2424::project_minimum_security(security_reports) == m2424::SecurityLevel::TC128;
    const bool bootstrap_report_ok = bootstrap_report.input.available
        && bootstrap_report.depth_boundary.available
        && bootstrap_report.successful_multiplications == 4
        && bootstrap_report.next_exponent == 16
        && !bootstrap_report.stop_reason.empty()
        && !bootstrap_report.preserve_value_criterion
        && !bootstrap_report.restore_level_criterion
        && has_stage(bootstrap_report, "ModRaise")
        && has_stage(bootstrap_report, "CoeffToSlot")
        && has_stage(bootstrap_report, "EvalMod")
        && has_stage(bootstrap_report, "SlotToCoeff");
    const auto impossible_scale_plan = m2424::plan_bootstrap_scale_strategy(
        {60, 40, 40, 40, 40, 40, 40},
        m2424::CipherInfo{std::pow(2.0, 40), 6, 7, 2, 300.0},
        -256.0,
        100.0,
        60.0,
        3);
    const auto feasible_scale_plan = m2424::plan_bootstrap_scale_strategy(
        {60, 40, 40, 40, 40, 40, 40, 40, 40, 40},
        m2424::CipherInfo{std::pow(2.0, 40), 9, 10, 2, 420.0},
        -256.0,
        50.0,
        60.0,
        3);
    const bool scale_strategy_plan_ok = !impossible_scale_plan.feasible
        && impossible_scale_plan.blocker == "not_enough_levels_for_scale"
        && impossible_scale_plan.required_drop_log2 == 236.0
        && impossible_scale_plan.available_drop_log2 == 120.0
        && feasible_scale_plan.feasible
        && feasible_scale_plan.total_levels_needed == 6
        && feasible_scale_plan.scale_after_squash_log2 <= 60.0;
    const auto evalmod_capacity_plan = m2424::plan_evalmod_first_multiply_capacity(
        {60, 40, 40, 40, 40, 40, 40, 40, 40, 40},
        feasible_scale_plan,
        2.0);
    const bool evalmod_capacity_plan_ok = evalmod_capacity_plan.first_multiply_ready
        && evalmod_capacity_plan.blocker == "none";

    const bool arithmetic_ok = std::isfinite(mul_accuracy.max_abs_error) && std::isfinite(mul_accuracy.mean_abs_error) && mul_accuracy.ok
        && close_enough(add_ref, add_out, 1e-5)
        && checked_add.ok
        && checked_add.operation == "add"
        && checked_add.info.ciphertext_size > 0
        && close_enough(input, sub_out, 1e-5)
        && close_enough(add_plain_ref, add_plain_out, 1e-5)
        && close_enough(input, sub_plain_out, 1e-5)
        && close_enough(mul_plain_ref, mul_plain_out, 1e-5)
        && close_enough(rot_ref, rot_out, 1e-5)
        && close_enough(input, loaded_out, 1e-5)
        && public_only_cannot_decrypt
        && !public_key_bytes.empty()
        && !secret_key_bytes.empty()
        && !relin_key_bytes.empty()
        && !galois_key_bytes.empty()
        && !cipher_bytes.empty();
    const bool bootstrap_parts_ok = close_enough(linear_ref, linear_out, 1e-4)
        && std::abs(sum_out.front() - sum_ref_value) < 1e-4
        && close_enough(polynomial_ref, polynomial_out, 1e-3);

    bool ok = arithmetic_ok
        && bootstrap_parts_ok
        && abft_check.ok
        && adapter.slot_count() == slots
        && adapter.serialized_size(ct) > 0
        && adapter.scale(ct) > 0.0
        && adapter.coeff_modulus_size(ct) > 0
        && adapter.public_key_size() > 0
        && adapter.relin_keys_size() > 0
        && adapter.galois_keys_size() > 0
        && expect_invalid_argument(invalid_profile_case)
        && expect_invalid_argument(empty_encode_case)
        && expect_runtime_error(encrypt_without_keys_case)
        && expect_runtime_error(mul_without_relin_case)
        && expect_runtime_error(rotate_without_galois_case)
        && expect_invalid_argument(accuracy_size_mismatch_case)
        && expect_invalid_argument(unknown_profile_name_case)
        && zero_transform_ok
        && all_named_profiles_create_contexts()
        && diagonal_transform_plan_ok()
        && eval_mod_polynomial_ok()
        && security_report_ok
        && bootstrap_report_ok
        && scale_strategy_plan_ok
        && evalmod_capacity_plan_ok;
    std::printf("[test_smoke] max=%.6e mean=%.6e arithmetic=%s bootstrap_parts=%s security=%s bootstrap_report=%s => %s\n",
               mul_accuracy.max_abs_error,
               mul_accuracy.mean_abs_error,
               arithmetic_ok ? "PASS" : "FAIL",
               bootstrap_parts_ok ? "PASS" : "FAIL",
               security_report_ok ? "PASS" : "FAIL",
               bootstrap_report_ok ? "PASS" : "FAIL",
               ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
