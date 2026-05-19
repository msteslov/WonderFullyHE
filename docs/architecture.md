# Architecture

The project implements a small protected-computation library on top of Microsoft SEAL.

```text
Application demos and experiments
        |
        v
m2424 library API
        |
        +-- SealAdapter: CKKS context, keys, encode/encrypt/evaluate/decrypt
        +-- accuracy: max/mean error reports
        +-- abft: checksum-based correctness checks
        +-- Bootstrapper: depth diagnostics and bootstrapping pipeline status
        +-- profile_report: reproducible CKKS parameter tables
        |
        v
Microsoft SEAL CKKS backend
```

## Public API layers

### SealAdapter

`SealAdapter` is the backend boundary. It hides direct SEAL types behind `Plain` and `Cipher` and exposes the CKKS operations used by the library:

- CKKS context creation from `CkksProfile`;
- key generation;
- encode/decode;
- encrypt/decrypt;
- add/subtract;
- multiply + relinearize + rescale;
- vector rotation;
- serialized object sizes;
- ciphertext diagnostics: `scale`, `chain_index`, `coeff_modulus_size`.

This keeps Microsoft SEAL as the backend dependency and leaves the public API under `m2424`.

### accuracy

The `accuracy` module defines the common numerical correctness criteria:

```text
max_abs_error
mean_abs_error
compare(expected, actual, tolerance)
```

The demos and tests use this module as the common accuracy criterion.

### abft

The `abft` module implements checksum-based correctness checks. The current checks cover:

- appended checksum for addition;
- appended checksum for subtraction;
- reference checksum for elementwise multiplication;
- sum preservation under rotation.

The ABFT layer is used to validate numerical consistency of protected computations.

### Bootstrapper

`Bootstrapper` is the bootstrapping entry point. The current implementation provides:

- multiplication-depth diagnostics;
- status of the CKKS bootstrapping pipeline;
- explicit stages: `ModRaise`, `CoeffToSlot`, `EvalMod`, `SlotToCoeff`.

The module separates depth analysis and bootstrapping pipeline state from the lower-level SEAL adapter.

### profile_report

`profile_report` makes CKKS parameter reporting reproducible:

```bash
./build/demo_profile_report
```

The output table contains `N`, slots, modulus chain, total modulus bits, scale, and estimated multiplication depth.

## Current demos

- `demo_secure_stats` shows protected cloud-style aggregation: sum and mean over encrypted data.
- `demo_abft` shows correctness checks for homomorphic operations.
- `demo_noise_growth` shows depth consumption and failure without bootstrapping.
- `demo_bootstrap_pipeline` shows current bootstrapping-module status.
- `bench_ckks` provides operation timing, numerical error, and serialized object sizes.
- `demo_profile_report` prints the CKKS parameter table.

## Next implementation steps

The next implementation step is to add reusable bootstrapping building blocks:

1. polynomial evaluation helpers for `EvalMod`;
2. rotation-based linear transform helpers for `CoeffToSlot` and `SlotToCoeff`;
3. richer ciphertext diagnostics for parameter and depth analysis.
