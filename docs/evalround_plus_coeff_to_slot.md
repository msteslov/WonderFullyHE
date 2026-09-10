# EvalRound+ CoeffToSlot (PR-1)

This change implements only the two linear branches from model v9 sections 3.2
and 3.5. It does not implement EvalRound, a nonlinear cleaner or a full bootstrap.

Given `y = sigma(u / Delta0)` after ModRaise, HP computes the original `F` and
`F D_phi` operators followed by conjugate/add. LP folds `alpha = Delta0/qSource`
into both `factors.front()` and `secondFirstFactor` during constant preparation.
The remaining factors, BSGS decomposition and rescale schedule are shared as
plan definitions; execution of HP and LP is independent from the same raised
ciphertext. There is no standalone alpha ciphertext multiplication or rescale.

## API and exact normalization

Include `m2424/evalround_plus_coeff_to_slot.hpp`. Resolve the source context before
ModRaise with `resolveBootstrapInput(adapter, sourceCipher)`. Pass that context,
the raised ciphertext and HP/LP contracts to `EvalRoundPlusCoeffToSlot::prepare`.
Then call `apply` with the same raised input and the prepared plan. Input magnitude
may be supplied as a deterministic bound; an omitted bound remains Unknown.

`CoeffToSlotPrefactor` retains the exact input scale bits and source prime factors.
Its default value is exactly one. The denominator product and rational rounding
use arbitrary-length integer limbs through the already vendored SEAL utilities;
there is no new GMP, MPFR, Boost or other runtime dependency. Each binary64 real
or imaginary diagonal component is multiplied by the exact rational and rounded
once to nearest, ties to even, at preparation. This does not claim that the
existing floating-point twiddle generation or SEAL encoding is exact.

A prefactor is part of prepared-plan identity. Old prepare/apply overloads still
mean prefactor one and cannot accidentally consume a prepared LP plan. The new
executor checks the actual ModRaise source base, context fingerprint and exact
binary64 scale before any homomorphic operation, along with existing degree,
levels and rotation/conjugation-key requirements.

## Trace and certificate limitations

Only `CoeffToSlotHP` and `CoeffToSlotLP` evidence is represented. No complete
`BootstrapCertificate` is constructed. Input compatibility success is distinct
from branch certification.

For each half, every factor and the final conjugate/add projection records the
actual scale, level and active primes, and bounds for kappa, diagonal error,
ideal magnitude, propagated semantic error and local added error. The recurrence
is v9 equation (6):

```
M_next <= kappa * M
E_next <= kappa * E + delta * (M + E) + B_local
```

Ideal kappa bounds follow from radix-2 row sums, first-factor normalization and
unit phase. Input implementation error is exactly zero relative to the supplied
semantic input `u`: as stipulated by v9, `nu_b` is part of `u`, not an extra linear
transform error. This does not assert that fresh ciphertext error is zero.

Remaining rigorous-certificate blockers:

- ideal-twiddle/composition and encoded-diagonal perturbation bounds;
- complete double-hoisted BSGS bounds for decomposition/key noise, inner/final
  ModDown, rescale and scale representation;
- the final conjugation key-switch/addition bound.

Until these are available, propagated and final HP/LP errors are Unknown,
`evidence.verified` is false, and each branch returns `RequiredBoundUnavailable`
with provenance. No observed CKKS error is promoted to an analytical certificate.
In particular, observed-versus-certified total error and headroom acceptance
cannot yet be established for either branch.

## Validation

Starting revision: `18cd31d121cbe4b4273e1dc4c1fedfc7f88c3c7b`.

Commands:

```
cmake --build build-pr0 -j
build-pr0/bin/test_evalround_plus_coeff_to_slot
ctest --test-dir build-pr0 --output-on-failure
cmake -S . -B build-analysis -DBUILD_TESTING=ON -DM2424_ENABLE_EVALMOD_ANALYSIS=ON
cmake --build build-analysis -j
ctest --test-dir build-analysis --output-on-failure
```

Baseline: 11/11 PASS. Analysis: 16/16 PASS. The existing
`tests/test_coeff_to_slot.cpp` was not modified. The new test covers:

- independent coefficient-half references for N=8,16,32, every basis vector and
  mixed vectors, multiple factorizations and both phase-corrected halves;
- exact rational oracle vectors, multi-prime products beyond binary64 range,
  subnormal tie-to-even behavior and avoiding premature scalar underflow;
- debug identity `u/qSource = I + (m+nu_b)/qSource`, with nonzero integer I;
- bitwise equality of HP and the legacy CoeffToSlot execution;
- exact source/scale/context binding, prefactor/factorization binding and missing
  rotation/conjugation-key rejection before execution;
- factor-level modulus traces, Unknown-bound propagation and no extra rescale.

The baseline test run observed HP max error `9.323388e-11` and LP max error
`9.236878e-11` (random encryption may vary these values). These are diagnostic
errors against independent decryption oracles, both below the unchanged
`2e-10` regression threshold used by the existing CoeffToSlot backend test.
Certified total HP/LP bounds remain Unknown, not zero and not these measurements.

At N=16384 with the test's 60-bit primes, both branches consume four levels
(240 modulus bits), with eight rescale operations per branch across two halves.
No tolerance or security parameter was changed. The initial input lowering is a
test fixture, separate from this consumption. The stale default build cache is
avoided using the existing clean `build-pr0` directory. The linker still emits
the environment's nonfatal warning about a missing Qt library search directory.
