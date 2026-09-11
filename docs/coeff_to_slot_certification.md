# Deterministic EvalRound+ CoeffToSlot certificate

Scope: follow-up to PR-1, using the PR-4 linear-transform preparation and
finite-support arithmetic. Mathematical source: v9 sections 3.2–3.3, especially
(4) and the factor error recurrence. No PR-5, probabilistic tail assumption,
new factorization optimizer, or full bootstrap certificate is introduced.

## API and compatibility

`EvalRoundPlusCoeffToSlot::prepareCertified` prepares a
`CertifiedEvalRoundPlusCoeffToSlot` consumed by the certified `apply` overload.
Preparation lives in the optional `m2424_evalmod_analysis` library. The core
executor, symbolic factors, contracts and traces do not depend on GMP/MPFR.

The original `prepare`/`apply` overloads remain diagnostic and retain their
original floating-point diagonal preparation and Unknown bounds. Existing
PR-0–PR-4 test sources and tolerances remain unchanged. Ordinary CoeffToSlot
continues to have prefactor one. The new preparation uses the same CtS
factorization, baby/giant rotations, level schedule and ideal operator.

The certified branch traces contain:

- `certifiedFactors`: the existing PR-4 factor/runtime trace types, with
  `LinearTransformFactorBound`, `LinearTransformPropagation`, exact scales,
  active primes, centered headroom and an arithmetic-term breakdown;
- `evidence`: verification of the available HP or LP error proof;
- `errorBudget`: a separate comparison with the requested branch accuracy;
- `outputError`: the maximum of the two half errors, including projection.

A valid error proof does not imply the requested error budget is feasible.
`ErrorBudgetExceeded` is retained even when the branch proof is verified.

The preparation binds the actual source primes, raised primes, special prime,
context fingerprint, exact binary64 scale and factorization. Missing rotation
or conjugation keys reject before arithmetic. Execution checks actual scale
bits, active primes and chain index after every kernel, rescale, conjugation
and addition. No metadata scale alignment is performed.

## Exact operators and shared arithmetic

The symbolic CtS factors are reverse adjoints of the existing PR-4 symbolic
forward factors, with reversed grouping. No root exponent is inferred from a
floating-point value. The second half folds the input phase
`zeta^(-r_j S)` into the first factor.

Both first HP factors have gain `1/N`. Both first LP factors have gain
`Delta0 / (N * product(sourcePrimes))`. `Delta0` is reconstructed as an exact
dyadic rational from its binary64 bits and the denominator is multiplied with
arbitrary-precision integers. The exact rational is passed directly to
`CertifiedRootDiagonalEncoder`; it is never first rounded to a scalar double.
Subsequent factors have gain one. There is no normalization ciphertext
multiply or extra rescale.

`prepareRootLinearTransform` is shared by CtS and StC. It computes:

- kappa = maximum row support times exact rational gain;
- delta = sum of rigorous encoded-diagonal perturbation bounds;
- B_local = weighted baby key noise + group inner ModDown + giant key noise
  + final ModDown + rescale rounding + exact scale representation errors.

The baby-noise multiplier is `kappa_baby + delta_baby`, including the exact
first-factor gain. The finite-support bounds and the internal QP/headroom
inequalities are exactly those described in
[the PR-4 baseline](slot_to_coeff_baseline.md). All factor propagation uses
`propagateLinearTransform`:

```
M_next <= kappa * M
E_next <= kappa * E + delta * (M + E) + B_local
```

The initial semantic input is exactly `y = sigma(u/Delta0)`; its source noise
is already in `u`. Hence the transform's initial semantic error is exactly
zero by definition, while its raised magnitude must be supplied with a known
analytical bound. This is not a replacement for an unknown backend bound.

The final real projection has kappa 2, delta 0 and nonzero local error
`finiteSupportKeyNoise + finiteSupportDivideRound` for conjugation. The
same-scale RNS addition is exact. Both conjugation and addition have their
own runtime state/error/headroom entries.

## LP domain

`CoeffToSlotCertificationInput` accepts bounds on the unscaled coefficient
magnitudes of `m` and `nu_b`. These must be deterministic, finite,
nonnegative, and have provenance. `certifyCoeffToSlotDomain` computes an
outward bound on

```
rho_cert = (messageMagnitude + sourceNoiseMagnitude) / q_src + E_LP
```

The denominator is the exact product of source primes. Missing upstream
bounds or missing E_LP yield `RequiredBoundUnavailable` and an Unknown rho;
`rho_cert >= 1/2` yields `DomainViolation`. This domain result is separate
from the HP/LP arithmetic proofs. It does not establish a lift bound K or a
complete EvalRound execution certificate.

## Measured baseline

Profile: N=16384, S=8192, seven actual 60-bit SEAL primes including the
special prime, input scale `2^59.5`. Both branches retain factorization
`[6,5,5,3]` of 19 raw butterfly/permutation stages. Each branch consumes four
levels, eight rescales over two halves, 134 ordinary rotations and two
conjugations. HP and LP execute independently.

The test encrypts the exact zero polynomial at the one-prime source level.
This still exercises nonzero ModRaise integer lifts. The upstream certificate
uses coefficient noise support `21(2N+1)+(N+1)/2`, obtained from asymmetric
SEAL encryption with ternary secret/ephemeral support and CBD support 21,
conservatively ignoring reduction by the encryption prime. The raised
embedding magnitude is bounded by `N(N+1)q_src/(2 Delta0)` from centered
ciphertext components. No decryption enters these bounds.

The following factor numbers are identical for both halves of each branch;
the second half's phase is separately encoded and certified. Row 4 is the
final conjugation/add projection.

| Branch | Factor | kappa | delta | B_local |
|---|---:|---:|---:|---:|
| HP | 0 | 0.00390625 | 2.27307162948250e-13 | 3.17064463871e-10 |
| HP | 1 | 32 | 3.48863026950771e-15 | 9.68780784531e-7 |
| HP | 2 | 4 | 2.61889922246048e-13 | 1.18004569897e-7 |
| HP | 3 | 1 | 1.07490855238859e-14 | 2.23487380403e-8 |
| HP | 4 | 2 | 0 | 1.39940808431e-8 |
| LP | 0 | 0.0027621358640285572 | 2.27484992937481e-13 | 2.72422180813e-10 |
| LP | 1 | 32 | 3.48863026950771e-15 | 9.68618020365e-7 |
| LP | 2 | 4 | 2.61889922246048e-13 | 1.15894436407e-7 |
| LP | 3 | 1 | 1.07490855238859e-14 | 2.19269894188e-8 |
| LP | 4 | 2 | 0 | 1.39940808431e-8 |

After factor rescales the chain indices are 4, 3, 2, 1, and the scales are
`8.1523861408341478e17`, `8.1523861408360013e17`,
`8.1523861408422579e17`, `8.152386140855465e17` respectively.
The stored schedule also contains the exact prime vectors and rational
arithmetic scales for the intermediate product/rescale states.

The final full-suite run produced:

| Branch | Certified error | Observed final coefficient-oracle error | Certified / observed |
|---|---:|---:|---:|
| HP | 0.0110686011341 | 1.01358921257e-10 | 1.09202041585e8 |
| LP | 0.0110729811405 | 8.42490522918e-11 | 1.31431521653e8 |

Observed values vary with encryption/key randomness. Every tested runtime
stage satisfies observed <= certified. Plaintext factor traces also match
independent canonical embedding basis references. The final backend error is
checked directly against raised/source coefficient debug oracles, with the
existing 2e-10 CtS diagnostic threshold. Stage reference computations are
floating-point diagnostics and are never used to produce a certificate.

`rho_cert = 0.0110729811411`, derived from the upstream coefficient noise
bound and E_LP. It passes rho < 1/2. Neither branch meets an error target of
1e-10; the deterministic baseline explicitly reports the accuracy failure.

## Comparison with PR-4

The PR-4 StC bound is `0.0340233987753`. The full CtS bounds are about
0.32532 and 0.32545 of that value, but their dominant error sources differ.
This is a comparison of the recorded profiles, not a claim that the operators
have the same input scales or gains: StC used scale 2^55.

Within CtS B_local, the largest term is factor 1 weighted baby key noise,
`9.68060426732e-7` for either branch. Factor 2 contributes
`1.10635477341e-7` from the same source. Final conjugation key noise is
`1.38294346675e-8`, plus `1.64646175587e-10` of ModDown rounding. Other
rescale rounding terms are `1.64646175587e-10` each. Inner/final ModDown at
the product scale and giant key noise are negligible here (about 1e-27 and
1e-25). Product scale representation is exactly zero for this dyadic schedule;
rescale representation is nonzero and explicitly included.

Weighting only the recorded B_local values by subsequent ideal kappas gives
about `8.12612e-6` (HP) and `8.10832e-6` (LP), approximately 1/4187 and
1/4196 of the recorded StC output bound. These are explanatory component
sums, not replacements for the certified recurrence, which also propagates
delta. The full CtS error is dominated by the first-factor encoded operator
perturbation applied to the conservative raised magnitude, then amplified by
subsequent factors. In contrast, the recorded PR-4 StC baseline was dominated
by propagated baby key noise. No term was reduced to make a target pass.

## Validation and remaining limitations

The new `test_coeff_to_slot_cert` covers exact-root plaintext basis cases,
every backend kernel/rescale/conjugation/add stage, exact scale/prime
schedules, folded normalization, provenance, known final bounds, altered
source/scale/context/factorization, missing keys, upstream Unknown states,
and exact domain division with a source product beyond binary64 range.

Validation commands and results:

- `cmake --build build-pr0 -j 4`: passed.
- `ctest --test-dir build-pr0 --output-on-failure`: 14/14, 15.04 s.
- `cmake --build build-analysis -j 4`: passed.
- Standalone new test via CTest: 1/1, 49.86 s.
- `ctest --test-dir build-analysis --output-on-failure`: 23/23, 148.31 s;
  includes the final new exact-domain and upstream-noise assertions.
- `git diff --check`: passed.

The full suite reproduced the PR-4 StC certificate `0.0340233987753` exactly.

Remaining limitations are the explicit deterministic accuracy failure,
required upstream certificates for arbitrary callers, and the optional
analysis dependency for rigorous preparation. All required CtS local and
operator bounds are known on the tested finite-support backend. No new
probabilistic tails or later optimization stages are part of this change.
