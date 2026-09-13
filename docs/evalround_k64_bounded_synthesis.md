# K=64 binary-digit bounded synthesis follow-up

Base revision: `fc6576edc6c3659a6b23e910ea8be939a8aab60f`.
This is the narrow PR-2/PR-3 feasibility follow-up. It does not change
`K=h_b=64`, the full `Bootstrapper` path, security parameters, tolerances, the
production modulus chain, or any probabilistic model.

## Fixed domain and search

The input radius is the exact binary64 value emitted by the unchanged selected
N=16384 sparse CtS certificate:

    rho binary64 bits = 0x3f0b460edc2fc0f3
    rho interpreted as binary64 = 5.2020386213659767e-5

The target table has 129 rows, one for every `I=-64,...,64`, and eight columns
equal to `bit_j(I+64)`. The deterministic degree sequence is exactly
`{64,128,192,256}`. The generic compiler limit remains degree 256.

The existing `remezOdd`/`MultiIntervalMinimax` generator cannot be applied to
this target: it constructs an odd polynomial for `x-round(x)`. Every one of its
32 digit/degree records is therefore `NotApplicableToTarget`, with no
coefficients and no certificate. Remez convergence is not reused as evidence.

The applicable `MultiIntervalChebyshev` path uses the existing Chebyshev and
exact-decimal infrastructure. Its explicitly documented candidate-generation
surrogate is the periodic C1 cosine bridge between consecutive exact binary
digit values. It equals the required target at every integer center and has
zero center derivative. This surrogate is only a generator: diagnostics are
evaluated against the constant digit target on the real cells. Chebyshev
coefficients are retained as decimal strings, converted exactly under the
substitution `u=x/64`, and reparsed as exact rationals. Only
`certifyEvalRoundDigitPolynomial` supplies evidence.

Each certificate uses 384-bit outward MPFR Horner arithmetic on 16 complete
closed subdivisions of each of all 129 intervals, including both endpoints.
The following is the complete applicable-family record. `|c|max` is the largest
absolute exact-decimal monomial coefficient; grid error is diagnostic only.

| digit | degree | generator | grid error | rigorous whole-domain error | `|c|max` | `a<=1` |
|---:|---:|---|---:|---:|---:|:---:|
| 0 | 64 | completed | 1.27619 | 1.8844606471484064e16 | 0.228389 | no |
| 0 | 128 | completed | 1.10059 | 1.5464366934178115e40 | 0.195721 | no |
| 0 | 192 | completed | 0.721592 | 1.3777042387361996e65 | 2.13078 | no |
| 0 | 256 | completed | 6.67715e-9 | 6.7951851457391055e77 | 2.46740 | no |
| 1 | 64 | completed | 1.34310 | 4.232628634043383e15 | 0.290463 | no |
| 1 | 128 | completed | 0.193873 | 9.5739985816692474e38 | 0.622421 | no |
| 1 | 192 | completed | 0.226696 | 5.0978976387014202e63 | 0.734392 | no |
| 1 | 256 | completed | 0.221052 | 4.0601689047061750e88 | 1.00025 | no |
| 2 | 64 | completed | 0.533444 | 2.2420002301215485e15 | 0.499098 | no |
| 2 | 128 | completed | 0.527137 | 1.0858889407516437e40 | 0.333837 | no |
| 2 | 192 | completed | 0.207767 | 6.1847956920313482e63 | 0.715018 | no |
| 2 | 256 | completed | 0.0774937 | 2.1071195281809190e88 | 1.21871 | no |
| 3 | 64 | completed | 0.842614 | 5.023928320119301e15 | 0.0823752 | no |
| 3 | 128 | completed | 0.377317 | 9.1769293489761857e38 | 0.480515 | no |
| 3 | 192 | completed | 0.219266 | 5.3569116916447001e63 | 0.750391 | no |
| 3 | 256 | completed | 0.123828 | 1.3331939453866094e88 | 1.14919 | no |
| 4 | 64 | completed | 0.644585 | 1.9117534474022958e15 | 0.292309 | no |
| 4 | 128 | completed | 0.375192 | 2.4161914126843178e39 | 0.395585 | no |
| 4 | 192 | completed | 0.191137 | 7.5417905308807853e63 | 0.758305 | no |
| 4 | 256 | completed | 0.111617 | 9.7259730140792264e87 | 1.16086 | no |
| 5 | 64 | completed | 0.712557 | 2.704889380879271e15 | 0.211220 | no |
| 5 | 128 | completed | 0.341207 | 2.3498949240093578e39 | 0.458832 | no |
| 5 | 192 | completed | 0.110145 | 6.2757116650570843e62 | 0.757122 | no |
| 5 | 256 | completed | 0.0727073 | 1.6474402244707921e88 | 1.16044 | no |
| 6 | 64 | completed | 0.701979 | 1.4533897409317772e15 | 1.00000 | no |
| 6 | 128 | completed | 0.336940 | 2.3561800699241836e39 | 1.00000 | no |
| 6 | 192 | completed | 0.0291100 | 4.9526575968336431e63 | 1.00000 | no |
| 6 | 256 | completed | 0.0524358 | 9.0519801489744206e87 | 1.16007 | no |
| 7 | 64 | completed | 0.0156735 | 3.5157608294825902e13 | 3.18821e-4 | no |
| 7 | 128 | completed | 0.00627599 | 5.3890723201993052e36 | 3.20691e-5 | no |
| 7 | 192 | completed | 0.00115074 | 1.0272788327292737e61 | 4.52037e-5 | no |
| 7 | 256 | completed | 0.000229096 | 8.6383513893953445e83 | 2.93271e-6 | no |

The best rigorous record for every digit is therefore degree 64 in the
MultiIntervalChebyshev family. Their maximum-coefficient magnitudes range from
`3.1882067136223852e-4` through `1.0000000000000029`.

The sharp gap between several small grid diagnostics and enormous rigorous
bounds is intentional evidence, not a tolerance failure. Exact conversion to
the compiler's monomial representation creates cancellation-heavy Horner
forms. The existing outward verifier correctly carries interval dependency
through those forms; a grid observation cannot replace or tighten that result.

## Cleaning and reconstruction diagnostics

All eight selected polynomials have finite deterministic outward bounds bound
to binary radix, the correct digit, K, exact rho bits, and exact coefficients.
Nevertheless, every selected bound exceeds the cleaner precondition `a<=1`.
Thus the eight *digit interval proofs* exist, but the K=64 extractor is **not
Certified**.

The unchanged PR-2 planner was run for local diagnostic budgets `1e-2`, `1e-4`,
and `1e-6`. All three runs return `CleaningDomainViolation`; there is no valid
cleaning schedule, no final finite `E_I`, and no certified post-cleaning weighted
contribution. Before cleaning, the eight weighted bounds `2^j a_j` are,
respectively, approximately:

    1.88446e16, 8.46526e15, 8.96800e15, 4.01914e16,
    3.05881e16, 8.65565e16, 9.30169e16, 4.50017e15.

These budgets are diagnostics only. None is the exhausted production backward
budget, and no production `requiredIntegerError` is fabricated.

## Ciphertext compilation feasibility

Because all eight interval certificates exist, the candidate was supplied to
the generic execution compiler in an analysis/test context. Its prerequisite
mathematical plan is rejected at `CleaningDomainViolation`, so the compiler
returns `ExtractionNotCertified` before constructing a ciphertext DAG. There
is consequently no honest extraction depth, cleaning depth, rescale count,
reconstruction level count, key list, or scale range to report for this
candidate.

The selected sparse N=16384 profile and its four-level CtS consumption are
unchanged. A level comparison cannot be made before a mathematical candidate
passes the cleaner-domain gate; the production chain is not enlarged. The
exact first remaining blocker is the cancellation-sensitive rigorous monomial
interval bound, which leaves every synthesized digit outside `a<=1`.

This is only a negative result inside the explicitly enumerated search. It is
not a global impossibility theorem and makes no production `1e-10` feasibility
claim.

## Regression coverage

`test_evalround_k64_synthesis` checks the eight-digit count, all 129 target
rows, all 64 digit/family/degree records, complete-cell outward proof metadata,
exact K/rho/digit binding, explicit bounded-search provenance, three local
planner budgets, and the analysis compiler's first blocker. The existing
generic polynomial test additionally checks grid-only rejection, missing digit,
wrong digit/radix/target/K/rho, Unknown bounds, and coefficient mutation with
certificate recomputation. K=1, sparse CtS, security inputs, and tolerances are
unchanged.

Final validation on this tree: baseline CTest `15/15` in 10.25 seconds;
analysis CTest `29/29` in 204.33 seconds. `git diff --check` is clean.
