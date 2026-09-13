# K=64 direct multi-interval Remez follow-up

Base revision: `2873ddff786ef53300afc3726bbb98360725fd04`.
This is the narrow PR-2/PR-3 follow-up requested after the bounded Chebyshev
diagnostic. It does not change `K=h_b=64`, security parameters, tolerances, the
production modulus chain, or the ciphertext backend.

## Fixed target and implementation audit

The domain uses the exact binary64 radius emitted by the unchanged selected
N=16384 sparse CtS certificate:

    rho binary64 bits = 0x3f0b460edc2fc0f3
    rho interpreted as binary64 = 5.2020386213659767e-5

For every `I=-64,...,64`, the closed cell is `[I-rho,I+rho]`; digit `j` has
the constant target `bit_j(I+64)`. Thus the direct target consists of 129
disjoint closed intervals and eight independently synthesized binary digits.

The pre-change `remezOdd` implementation used 384-bit MPFR, a 2048-point grid
per nonnegative interval, odd powers only, Gaussian elimination with partial
pivoting, and at most 24 exchange iterations. Its convergence diagnostic
required alternating signs, amplitude ratio at most 1.25, stable exchange
points, and stabilization of the dense-grid maximum. The refactor extracts
those mechanics into a target-agnostic exchange core. The legacy wrapper still
constructs the same nearest-integer residual target, odd basis, grid, scale,
and iteration limit; a regression test locks its exact degree-7 decimal
coefficient vector.

The new public experimental request supplies ordered disjoint intervals with
an exact-decimal constant target for each interval, a full degree basis, a
variable scale, grid density, and iteration bound. No grid point or exchange
point is inserted in an excluded gap. Candidate generation and its convergence
flag remain diagnostics: only the existing outward whole-cell verifier creates
evidence.

## Bounded direct search

Every one of the 32 direct Remez attempts uses degrees
`{64,128,192,256}`, the full Chebyshev basis in `u=x/64`, eight samples per
interval including both endpoints, at most four exchange iterations, and
384-bit MPFR. Exact Chebyshev-to-monomial conversion is followed by the existing
384-bit outward Horner verifier on 16 complete closed subdivisions of every
cell. `grid` below is the maximum on the generator grid only; `rigorous` is the
authoritative outward whole-domain bound. `conv` records the exchange
heuristic, and every final exchange point was checked to be inside the union.

| digit | degree | conv | iterations | grid | rigorous | `|c|max` | `a<=1` |
|---:|---:|:---:|---:|---:|---:|---:|:---:|
| 0 | 64 | no | 4 | 3.2533603438026103e2 | 6.2912160215034224e7 | 4.9999999999999989e-1 | no |
| 0 | 128 | yes | 3 | 5.0000000000000011e-1 | 3.8026446126520688e21 | 4.9999999999999994e-1 | no |
| 0 | 192 | no | 4 | 3.2963055693767237e39 | 3.8808997656195933e73 | 2.4674011130470905 | no |
| 0 | 256 | no | 4 | 4.5825178144617754e44 | 1.8066658370477226e95 | 2.4674011170836114 | no |
| 1 | 64 | no | 4 | 5.0e-1 | 4.3434017547545503e4 | 4.9999999999999994e-1 | no |
| 1 | 128 | no | 4 | 8.3984213964768277e25 | 2.1536438776193502e47 | 7.8535402906658780e-1 | no |
| 1 | 192 | no | 4 | 1.2224609483815205e43 | 1.1678660748397883e77 | 7.8536146240858906e-1 | no |
| 1 | 256 | no | 4 | 2.4299366217084288e56 | 1.6151787810516606e107 | 2.6635298600619994 | no |
| 2 | 64 | no | 4 | 3.7371596675247156e14 | 3.2272585697530552e19 | 5.5552025405354355e-1 | no |
| 2 | 128 | no | 4 | 3.3870603593307725e26 | 8.0292532389339120e48 | 7.7096007924611476e-1 | no |
| 2 | 192 | no | 4 | 9.7324618230550034e47 | 4.3950827555499572e81 | 1.3172464216152486 | no |
| 2 | 256 | no | 4 | 1.6171830305838514e62 | 7.8966113690148319e112 | 3.3835907597467787 | no |
| 3 | 64 | no | 4 | 2.0256620854255381e14 | 6.6557956058117104e19 | 6.1247139623522351e-1 | no |
| 3 | 128 | no | 4 | 3.2513787938994960e28 | 1.2250518191254292e51 | 7.5894845318492410e-1 | no |
| 3 | 192 | no | 4 | 7.6058559851664490e49 | 4.1402296949544553e83 | 1.1659091676667146 | no |
| 3 | 256 | no | 4 | 8.8249711228095742e62 | 4.0198011316886945e113 | 2.8648953754427149 | no |
| 4 | 64 | no | 4 | 1.1213413187094935e15 | 7.0102255109930353e20 | 6.3271826446626067e-1 | no |
| 4 | 128 | no | 4 | 2.9889356934351709e28 | 1.0060579563037141e51 | 6.8915001174595214e-1 | no |
| 4 | 192 | no | 4 | 1.0288541584821227e49 | 8.1565404593686265e82 | 1.1879913927167900 | no |
| 4 | 256 | no | 4 | 4.9215325695960019e62 | 2.2495567295914771e113 | 2.5787165614763876 | no |
| 5 | 64 | no | 4 | 1.0567621634115525e14 | 4.5172305064149071e20 | 6.3292265913838319e-1 | no |
| 5 | 128 | no | 4 | 2.0968821778324997e28 | 8.0345063354862726e50 | 6.8893446044632345e-1 | no |
| 5 | 192 | no | 4 | 3.1817316747952843e48 | 4.2210197955120540e82 | 1.6853513495088119 | no |
| 5 | 256 | no | 4 | 4.8691616331980730e62 | 2.2284374730004153e113 | 2.5779546844259338 | no |
| 6 | 64 | no | 4 | 1.0721604207601384e14 | 4.5115245586003827e20 | 8.2859829974668564e-1 | no |
| 6 | 128 | no | 4 | 2.0968178369111288e28 | 8.0344733912617469e50 | 9.6434318059292667e-1 | no |
| 6 | 192 | no | 4 | 3.1817316708390605e48 | 4.2210197958171173e82 | 1.6853513495538084 | no |
| 6 | 256 | no | 4 | 4.8691616331944403e62 | 2.2284374729994779e113 | 2.5779546844259222 | no |
| 7 | 64 | yes | 3 | 6.2967015571384372e-2 | 1.2260097184914792e4 | 3.9381532810206875e-19 | no |
| 7 | 128 | no | 4 | 5.6873231029276063e-4 | 9.5773958568301420e15 | 4.0106277130404384e-38 | no |
| 7 | 192 | yes | 4 | 2.9242507118296329e-7 | 2.0394312770179371e29 | 4.0933209965213033e-55 | no |
| 7 | 256 | yes | 4 | 2.8084032784527216e-7 | 2.9417790095985345e41 | 2.0603661661945365e-74 | no |

The earlier cosine-bridge family remains only as a separately certified
diagnostic comparator; it is never used as the Remez target. Selection uses the
smallest rigorous bound, not convergence or grid error. Direct degree-64 Remez
wins digits 0, 1, and 7; degree-64 cosine diagnostics remain smaller for digits
2 through 6. The selected rigorous bounds are respectively
`6.2912160215034224e7`, `4.3434017547545503e4`,
`2.2420002301215485e15`, `5.023928320119301e15`,
`1.9117534474022958e15`, `2.7048893808792710e15`,
`1.4533897409317772e15`, and `1.2260097184914792e4`.

## Feasibility result and strict gate

All 32 direct candidates have finite outward certificates, and several Remez
runs report convergence, but every direct bound exceeds `a_j<=1`. Every selected
cross-family bound also exceeds that cleaner precondition. Convergence and small
sampled error therefore do not certify the extractor.

The model's cleaner recurrence `a_next <= 5 a^2 + B` is applicable only after
the `a<=1` input-domain gate (and its zero-local-error contraction threshold is
stricter, `a<1/5`). Because all eight selected digits do not pass `a_j<=1`, the
implementation does not construct an `EvalRoundCandidate`, does not run the
cleaning planner for any local budget, and does not compile a ciphertext DAG.
Consequently there is no valid cleaning schedule, final `E_I`, depth, rescale
count, reconstruction level count, key list, or scale range to report.

The exact first blocker is `CleaningDomainViolation`. This is a negative result
only for the explicitly bounded search, not a global impossibility theorem.

**No Certified direct multi-interval Remez extractor was found in this bounded
degree search.**

## Regression coverage

`test_evalmod_synthesis` locks the legacy degree-7 `remezOdd` coefficients.
`test_evalround_k64_synthesis` checks all 129 target rows, all 32 direct
digit/degree attempts, union-only exchange points, convergence metadata,
outward proof provenance and exact K/rho/digit binding, rigorous cross-family
selection, and the strict no-candidate/no-planner gate. The existing generic
polynomial tests retain missing-digit, wrong digit/radix/target/K/rho, grid-only,
unknown-bound, and coefficient-mutation rejection coverage.
