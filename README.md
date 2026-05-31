# WonderFullyHE

WonderFullyHE is a C++17 CKKS research library built on top of Microsoft SEAL.
The stable part of the project is the `m2424` library layer: SEAL adapter,
accuracy metrics, ABFT checks, operation-budget tracking, parameter planning,
linear transforms, polynomial evaluation, and profile/security reports.

Bootstrapping is still experimental. The repository contains implementation
pieces and diagnostic harnesses for `ModRaise -> CoeffToSlot -> EvalMod ->
SlotToCoeff`, but full production bootstrapping is not yet a stable API.

## Build

Clone with the SEAL submodule:

```bash
git clone --recurse-submodules <repo-url>
```

Configure and build the library, tests, and stable apps:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Run the smallest demo:

```bash
./build/demo_basic
```

## Build Options

`M2424_BUILD_APPS=ON` builds stable command-line demos and measurement tools.
It is enabled by default.

`M2424_BUILD_RESEARCH_APPS=OFF` keeps experimental bootstrap harnesses out of
the default build. Enable it only when working on bootstrap internals:

```bash
cmake -S . -B build -DM2424_BUILD_RESEARCH_APPS=ON
cmake --build build --target bench_bootstrap_multi_cycle -j
```

## Stable Apps

- `demo_basic` - minimal encode/encrypt/evaluate/decrypt demo.
- `demo_abft` - checksum-style ABFT checks for encrypted operations.
- `demo_checked_pipeline` - checked `add -> mul -> rotate -> sum_slots` flow.
- `demo_client_compute_roundtrip` - client/server-style key and ciphertext flow.
- `demo_eval_mod_polynomial` - standalone EvalMod polynomial check.
- `demo_galois_key_optimization` - rotation-key size/runtime comparison.
- `demo_mod_raise` - low-level ModRaise diagnostic.
- `demo_noise_growth` - depth and error growth without bootstrapping.
- `demo_precision_profiles` - profile comparison.
- `demo_profile_report` - profile table.
- `demo_secure_stats` - encrypted sum/mean example.
- `demo_security_report` - SEAL security-limit report.
- `bench_ckks` - core CKKS operation timings and error.
- `bench_chain_accuracy` - chain length, scale, work-bit accuracy sweep.
- `bench_parameter_planner` - planner sanity check on SEAL runs.
- `bench_parallel_throughput` - independent-ciphertext throughput.

## Tests

Default CTest runs invariant tests only:

- `test_smoke`
- `test_accuracy`
- `test_adapter_failures`
- `test_bootstrap_dft`
- `test_bootstrap_scalable_refresh`
- `test_bootstrap_scale_design`
- `test_operation_budget`
- `test_parameter_planner`

## Current Bootstrap Status

The codebase has reusable bootstrap building blocks:

- FFT-like `CoeffToSlot` / `SlotToCoeff` plans.
- `ModRaise` wrapper over SEAL-level primitives.
- EvalMod P3/P5/P7 plus a diagnostic `P3DoubleAngle` path.
- Scale and period planning objects.
- A public experimental `SlotsToCoeffsFirst/FftLike/P3` refresh path.

The current blocker is the period/output-scale model around EvalMod. Existing
research harnesses are kept behind `M2424_BUILD_RESEARCH_APPS` until that model
is correct enough for a stable full bootstrap API.

## API Sketch

```cpp
#include "m2424/m2424.hpp"

auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
adapter.keygen(true, true);

auto encrypted = adapter.encrypt(adapter.encode({1.0, 2.0, 3.0}));
auto squared = adapter.mul_relin_rescale(encrypted, encrypted);
auto decoded = adapter.decode(adapter.decrypt(squared));
```

Useful documentation:

- `docs/architecture.md`
- `docs/ckks_parameters.md`
- `docs/project_status.md`
