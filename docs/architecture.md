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

`SealAdapter` is the backend boundary. It hides direct SEAL types behind `Plain` and `Cipher` and exposes only the operations currently needed by the project:

- CKKS context creation from `CkksProfile`;
- key generation;
- encode/decode;
- encrypt/decrypt;
- add/subtract;
- multiply + relinearize + rescale;
- vector rotation;
- serialized object sizes;
- ciphertext diagnostics: `scale`, `chain_index`, `coeff_modulus_size`.

This keeps Microsoft SEAL as an implementation detail instead of making the whole SEAL API part of the project API.

### accuracy

The `accuracy` module defines the common numerical correctness criteria:

```text
max_abs_error
mean_abs_error
compare(expected, actual, tolerance)
```

All demos and tests should use this module instead of duplicating their own error logic.

### abft

The `abft` module implements checksum-based correctness checks. The current checks cover:

- appended checksum for addition;
- appended checksum for subtraction;
- reference checksum for elementwise multiplication;
- sum preservation under rotation.

ABFT does not improve cryptographic security. Its role is to detect incorrect computation results or unacceptable numerical drift.

### Bootstrapper

`Bootstrapper` is the current bootstrapping entry point. It does not yet return a refreshed ciphertext. It currently provides:

- multiplication-depth diagnostics;
- status of the planned CKKS bootstrapping pipeline;
- explicit stages: `ModRaise`, `CoeffToSlot`, `EvalMod`, `SlotToCoeff`.

This lets the project present bootstrapping as an active module while the mathematically complete implementation is developed stage by stage.

### profile_report

`profile_report` makes CKKS parameter reporting reproducible. Instead of copying parameters manually into slides, run:

```bash
./build/demo_profile_report
```

The output table is the source for `N`, slots, modulus chain, total modulus bits, scale, and estimated multiplication depth.

## Current demos

- `demo_secure_stats` shows protected cloud-style aggregation: sum and mean over encrypted data.
- `demo_abft` shows correctness checks for homomorphic operations.
- `demo_noise_growth` shows depth consumption and failure without bootstrapping.
- `demo_bootstrap_pipeline` shows current bootstrapping-module status.
- `bench_ckks` provides operation timing, numerical error, and serialized object sizes.
- `demo_profile_report` prints the CKKS parameter table used in documentation and slides.

## Development direction

The next implementation step is to add bootstrapping building blocks that do not depend on the final full model:

1. polynomial evaluation helpers for future `EvalMod`;
2. rotation-based linear transform helpers for future `CoeffToSlot` and `SlotToCoeff`;
3. richer ciphertext diagnostics for parameter and depth analysis.
