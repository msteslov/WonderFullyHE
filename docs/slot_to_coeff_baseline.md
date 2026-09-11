# Standalone SlotToCoeff baseline (PR-4)

The target is v9 §5.2: `H[j,k]=zeta^(r_j*k)` and
`T(a0,a1)=H*a0+D_plus*H*a1`, with `D_plus[j,j]=zeta^(r_j*S)`.
There is no normalization after the transform and no connection to EvalRound
outputs in this change. PR-1 code, numerical behavior and outstanding certificates
are unchanged.

`SlotToCoeffPlan` (also exposed as `SlotToCoeff`) provides mathematical planning,
`prepare`, immutable `PreparedSlotToCoeffPlan`, `preflight` and `apply`. Preparation
is implemented in the optional `m2424_evalmod_analysis` library. The core plan,
plaintext path, certificate recurrence, preflight and executor do not acquire a
GMP/MPFR runtime dependency. Link the optional library to prepare rigorous plans.

Prepared data binds degree, factorization, context fingerprint, actual active
primes/chain index, both exact binary64 input scales, input component counts,
prepared diagonal plaintexts and key requirements. The baseline requires equal
input scales and active primes, rejecting mismatches rather than rewriting scales.
All diagonal preparation occurs before execution. Apply neither encodes nor
generates keys; it verifies the actual scale bits, chain index and primes after
every double-hoisted kernel, rescale and final addition. Tests can observe these
ciphertexts through the callback; production does not decrypt.

Two-component inputs need rotation keys only. Three-component inputs additionally
require explicit input relinearization, performed and certified before the linear
factors. Missing relin keys are rejected before preparation/execution. Pure linear
factors themselves do not perform ciphertext multiplication or relinearization.
`requirements()` describes the mathematical linear plan; prepared
`certificate().keys` additionally reflects the actual input component counts.

## Mathematical factors and the existing oracle

The factorization is the reverse sequence of adjoints of the unnormalized
canonical inverse FFT: bit swaps followed by conjugate-transposed radix-2
butterflies, grouped by the requested positive radices. The phase `D_plus` is
folded into the output rows of the last factor of the second branch.

Factors are represented symbolically: every nonzero entry is an exact power of
`zeta`; a separate sentinel denotes zero. Composition checks that each resulting
entry has one FFT path. Thus `kappa_r` is exactly the maximum nonzero row count,
not a floating-point estimate of unit-modulus coefficients. Grouping is explicit;
no new factorization optimizer, AKS, SCORE, Tuple-S2C, Grafting or Thrifty path is
introduced. The BSGS split uses the existing power-of-two baby/giant strategy.

The repository's existing oracle names describe the opposite direction from this
bootstrap-stage name: `slotToCoeffReference(slots)` returns real coefficients,
while `coeffToSlotReference(coefficients)` returns their canonical embedding.
Neither oracle was renamed or duplicated. Tests use the former to recover the
input coefficient halves from StC output, and the latter to check forward output.
Complex halves use linearity and the two existing real-coefficient oracle calls.
Small tests cover N=8,16,32, every basis column of both halves, mixed complex
vectors, several groupings, and `StC(CtS(a))` independently of EvalRound.

## Certified diagonal encoding and operator perturbation

`CertifiedRootDiagonalEncoder` is reusable outside StC. It takes exact root
exponents, a rational prefactor and a binary64 scale, interpreted as an exact
dyadic rational. The rational prefactor can describe normalized/folded factors
when this helper is later applied to CtS. No such PR-1 integration occurs here.

The inverse embedding is computed with a radix-2 FFT at MPFR precision p=192.
Let `u=2^-p`. Correctly rounded pi, multiplication by an integer, exact power-of-two
angle division and correctly rounded sin/cos imply complex root error `w<=64u`.
In detail, the angle error is below 16u for angles in [0,2*pi]; each trigonometric
component adds at most u, and complex modulus is bounded by the sum of component
errors. This is comfortably below 64u.

For a butterfly with ideal input magnitude M and accumulated complex error E,
the point-arithmetic error is bounded conservatively by
`64u*(M+E)*(1+w)`: four rounded real products, two sums for complex multiplication
and two final butterfly additions fit this bound. Therefore preparation propagates

```
M_next = 2*M
E_next = 2*E + w*(M+E) + 64u*(M+E)*(1+w).
```

The starting root error is w. The final unit-root twist gets the corresponding
multiplication bound; division by N is exact in MPFR. A wide MPFR exponent range
is required, so underflow/overflow cannot invalidate these estimates.

The resulting coefficient midpoint is extracted as an exact rational. After
multiplying by the rational prefactor and exact dyadic encoding scale, integer
rounding and its distance from that midpoint are computed with GMP. For every
coefficient, its error is bounded by this exact rounding distance plus the proven
FFT error. The sum of coefficient errors bounds the error in every canonical
slot. This yields a certificate for the actual encoded diagonal, without relying
on SEAL's floating-point encoder or measured decoding error.

Integers are written directly to reduced RNS limbs and transformed by exact NTT
at key level. `modSwitchPlainTo` permits exact prefix restriction of such a
prepared plaintext when needed for test inputs; this avoids SEAL's public
plaintext mod-switch rejection of pure key-level plaintexts. No rescale or
re-encoding is involved.

For a factor, the sum of its diagonal sup-error bounds is `delta_r`, a rigorous
bound on `||Bhat_r-B_r||_inf`. It includes constant preparation/encoding only.
It is never merged into the backend-local bound.

## Double-hoisted arithmetic proof

PR-3 and PR-4 use the same extracted finite-support arithmetic helpers. With N
coefficients, ternary secret support 1, evaluation-key coefficient support B=21,
active primes q_i, special prime P and scale s, define

```
K(s) = N^2 * B * sum(q_i-1) / (P*s)
R(s) = N*(1+N) / (2*s).
```

K is key-switch noise before component rounding. R is two-component divide-and-
round noise; the generic helper also supports the other component counts. CBD
support follows from the difference of two 21-bit Hamming weights. Unknown,
probabilistic or insufficient support assumptions are rejected. Gaussian builds
remain unsupported by this deterministic certificate. Ciphertexts and evaluation
keys must share the expected ternary secret and finite-support key-generation
model; preflight cannot establish these distribution assumptions by inspection.

The native double-hoisted kernel extends baby rotations to QP, multiplies and
accumulates prepared diagonals there, performs an inner ModDown per giant group,
extends/rotates each group again, sums in QP and performs one final ModDown.

Let G be the number of nonempty groups and g the number of nonzero giant steps.
Let beta be the maximum row count of ideal diagonals with nonzero baby step plus
the sum of those diagonals' perturbation bounds. All baby noises have bound K(s),
and giant permutations restore the original rows, so their weighted contribution
is bounded by `beta*K(s)`, without assuming independence or counting a zero-step
rotation as a key switch.

For plaintext scale t and runtime product scale v=fl(s*t), the kernel-local bound is

```
B_double = beta*K(s) + (G+1)*R(v) + g*K(v) + B_product_scale,
B_product_scale = abs(s*t/v-1)*((kappa+delta)*(M+E)+beta*K(s)).
```

Every ratio is evaluated as an exact rational. After dropping the actual integer
prime q and obtaining runtime scale v_next, rescale adds

```
R(v_next) + abs(v/(q*v_next)-1)*(kappa*M + E_double).
```

This includes the binary64 conversion of q and division rounding. The sum of the
kernel and rescale additions is `B_loc`. Identity RNS operations can have a proven
zero local bound; key switches, ModDown and rescale are never replaced by synthetic
zero bounds.

Preparation checks positive exact centered headroom for input, baby extensions,
every group accumulator/inner ModDown/giant extension, the outer QP accumulator,
final ModDown, rescale and pair addition. QP checks use the equivalent
P-normalized inequality; `Q-2*ceil(scale*norm)>0` is sufficient for the extended
QP check as well. The public trace records the internal proofs in addition to
observable kernel/rescale stages. Group envelopes conservatively dominate every
group rather than assuming cancellation between groups.

## Reusable recurrence and gain

`propagateLinearTransform` retains separate contributions:

```
M_next <= kappa*M
E_next <= kappa*E + delta*(M+E) + B_loc.
```

All nonnegative floating-point bound operations round outward. Unknown required
bounds prevent certification. At each factor, the direct-sum two-branch operator
uses the maximum branch delta. The final pair addition has kappa=2, delta=0 and
zero arithmetic error at equal scales. `linearTransformGain` computes the product
of these actual `(kappa+delta)` bounds. Thus the factor of two is accounted for
explicitly, and N is only the ideal reference gain, not the implementation
certificate. Additive runtime errors stay outside Gamma and inside the recurrence.

The result exposes `certificate.factors` with per-factor and runtime traces,
`outputError`, `gamma`, `metrics.depth` (levels consumed per branch), total rotation
and rescale counts, exact input/output scales and required keys. There is no
bootstrap-wide Certified result and no EvalRound-to-StC connection.

## Remaining scope limits

Rigorous preparation requires the optional analysis library. The selected CBD
baseline has no Unknown required bounds, but finite-support key-switch estimates
are conservative and do not claim that a later full-bootstrap accuracy budget is
met. Unknown upstream errors, unsupported noise models and insufficient headroom
or levels are explicit rejection reasons. Outstanding PR-1 bounds are untouched.
PR-5 and subsequent stages are not implemented.

## Recorded validation

Build and test commands:

```
cmake --build build-pr0 -j 4
ctest --test-dir build-pr0 --output-on-failure
cmake --build build-analysis -j 4
ctest --test-dir build-analysis --output-on-failure
```

Baseline: 14/14 passed. Analysis: 22/22 passed, including the complete old PR-0
through PR-3 suites and both new StC tests. The new tests also passed independently.
No existing test source/tolerance was changed. Negative tests cover unknown
required bounds, insufficient levels, missing rotations and input relin keys,
changed factorization/context/input scales/active primes, invalid rounding mode,
and headroom violations. The size-three input path is executed as well as checked
for missing-key rejection.

Representative size-two profile: N=16384, six 60-bit data primes and a 60-bit
special prime, input scale 2^55. Grouping: `[5,5,5,4]` over six bit-swap and thirteen
butterfly stages. There are 314 rotation operations, 8 rescale operations across
the two branches, 4 levels consumed per branch, and no relinearization for these
size-two inputs. The additional size-three test performs two input
relinearizations without consuming extra levels.

| Factor | kappa (each branch) | delta (each branch) | B_local (each branch) |
| --- | ---: | ---: | ---: |
| 0 | 1 | 4.2419277866743739e-13 | 9.42498672883e-7 |
| 1 | 16 | 1.1047031066973692e-13 | 1.25207009205e-5 |
| 2 | 32 | 1.3898735138957541e-14 | 1.75274910816e-5 |
| 3 | 16 | 5.6932836635238736e-14 | 5.63636444892e-6 |
| Pair addition | 2 | 0 (no encoding) | 0 (exact RNS addition) |

Gamma: `16384.00000000716`. Certified output error: `0.0340233987753`.
Observed maximum output error in the full-suite run: `7.98328499403e-10`.
Observed errors vary with freshly generated keys and are never used to obtain
these bounds. The gap is due to deterministic worst-case key-switch supports,
not an empirical tolerance. This certificate alone does not establish a future
full-bootstrap accuracy target.

The actual data chain index is `5 -> 4 -> 3 -> 2 -> 1`. Input scale is
`36028797018963968`; successive post-rescale scales are
`36028797018969088`, `36028797018977280`, `36028797019004928`, and
`36028797019063296`. Both branches and their final addition have identical scale
bits. The prepared trace records the actual primes and exact rational arithmetic
scale at every stage, including the intermediate plaintext product scales.

Changed files:

* `include/m2424/slot_to_coeff.hpp`
* `include/m2424/linear_transform_contract.hpp`
* `include/m2424/experimental/evalmod_analysis/certified_diagonal.hpp`
* `include/m2424/experimental/evalmod_analysis/finite_support_arithmetic.hpp`
* `include/m2424/seal_adapter.hpp`
* `include/m2424/m2424.hpp`
* `src/math/slot_to_coeff.cpp`
* `src/math/certified_diagonal.cpp`
* `src/core/slot_to_coeff_internal.hpp`
* `src/core/linear_transform_contract.cpp`
* `src/planning/slot_to_coeff_preparation.cpp`
* `src/planning/evalround_execution.cpp` (shared helpers only)
* `src/ckks/seal_adapter.cpp`
* `src/CMakeLists.txt`
* `tests/test_slot_to_coeff.cpp`
* `tests/test_slot_to_coeff_contract.cpp`
* `tests/CMakeLists.txt`
* `docs/slot_to_coeff_baseline.md`
