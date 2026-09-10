# EvalRound reference and planner (PR-2)

Source: `ckks_bootstrapping_model_v9.pdf`, section 4. Starting revision:
`133fb5953c12567e5fdd775aa780fbdab28e1710`.

## Scope and API

`m2424/evalround.hpp` defines the public problem, radix, extraction description,
per-digit bounds, configured cost, plan/certificate and rejection reasons.
`EvalRoundProblem::requiredIntegerError` must be supplied by the caller. There is
no automatic choice of `1e-10`. The input contract assumes a certified integer
bound K and radius rho supplied by the upstream stage; PR-2 does not certify that
a ciphertext satisfies this domain.

`planEvalRound(problem, candidates)` selects the lowest configured cost among
certified candidates. Certification here is a conditional mathematical/reference
error certificate, not ciphertext executability or a complete bootstrap certificate.
No EvalRound ciphertext operations, key generation, scale/modulus scheduling,
Grafting, Thrifty evaluation or other PR-3 work is implemented.

The core planner depends only on existing core types. The optional
`m2424_evalmod_analysis` target supplies the GMP/MPFR reference through
`experimental/evalmod_analysis/evalround_reference.hpp`, reusing `ExactInteger`
and `ExactScale`. GMP/MPFR are not new mandatory core runtime dependencies.

## Representations and extraction

Binary uses J=I+K and the minimum digit count with 2^k >= 2K+1. Reconstruction is
sum(2^j b_j)-K. It never uses two's complement. Balanced ternary uses the minimum
count with 3^k >= 2K+1 and digits in {-1,0,1}. Integer count/reconstruction code
also supports the endpoints of the full uint32 K range without floating-point
logarithms or signed overflow. K=0 has zero digits.

Whole-domain extraction certificates bind both K and rho. Factory reference
candidates include:

- `PiecewiseReference`, for arbitrary K: the intervals are disjoint for rho<1/2,
  so the exact center uniquely determines each binary digit or root-of-unity trit.
  Exact mathematical target error is zero.
- `BinaryQuadraticK1`: b0(x)=1-x^2 and b1(x)=(x^2+x)/2. At I=-1,0,1 these encode
  I+1. Expansion of x=I+xi on all three intervals gives errors at most
  2*rho+rho^2 and (3*rho+rho^2)/2 respectively.
- `TernaryPhaseReferenceK1`: exp(2*pi*i*x/3), the one-trit phase route from v9 4.2.
  The whole-domain error is at most 2*pi*rho/3, using |exp(i*t)-1|<=|t|.
  This is an exact-function reference, not a polynomial approximation to exp.

Factory local cleaning/reconstruction errors are zero only for the explicitly
exact-arithmetic mathematical reference. They do not assert zero backend error
or zero finite-precision MPFR rounding. The reference evaluator records MPFR
values as exact dyadic rationals after rounding at each stage; these observations
are diagnostics and never become analytical bounds. Exact digit reconstruction
is tested separately from finite-precision complex-root reconstruction.

`ExternalPolynomial` and `DigitExtract` describe externally supplied extractors.
They require verified whole-domain evidence and known per-digit bounds. This PR
does not synthesize their approximation polynomials or implement the literature's
DigitExtract reference circuit. Without evidence they are diagnostic/rejected,
never Certified. DigitExtract additionally requires precisely
0<rho<=1/8 and 2*rho<epsilon_de<=1/4. These restrictions do not apply to other
methods. In particular, rho=1/8 cannot admit epsilon_de under the strict inequality.

## Cleaning, reconstruction and minimality

The exact GMP reference implements f2(x)=3*x^2-2*x^3 and
f3(x,y)=(y^2+4*x-2*x^2*y)/3 with y=conjugate(x).

The planner uses only the hard analytical recurrences from v9 4.3:

```
binary:  a_next <= 5*a^2 + B_local
ternary: a_next <= 3*a^2 + B_local
```

Each reachable cleaning step requires 0<=a<=1 and a known local bound. Unknown
bounds are never zero. Unused future cleaning rounds do not need a bound.
Reconstruction follows v9 4.5:

```
E_binary  <= sum_j 2^j * (a_j + B_rec_j)
E_ternary <= sum_j 3^j * ((2/sqrt(3))*a_j + B_rec_j)
```

Positive arithmetic is rounded outward, retaining exact zero and exact normal
power-of-two products. The ternary gain is a fixed upward binary64 enclosure of
2/sqrt(3). Unsupported rounding modes and fast-math reject certification.

The search enumerates per-digit error schedules up to the configured maximum
(default 16, maximum 64 rounds per digit). Dynamic programming retains the lowest
weighted error for each total iteration count. This is sufficient because future
contributions are independent and nonnegative. The first feasible total count is
minimal within this search; counts can differ across digits. Each candidate's
cost is extraction + reconstruction + cleaningPerDigitIteration*totalIterations;
configured stage costs can include all applicable stage work. There is no claim
of global impossibility when this bounded search fails.

Selected plans retain the input budget/domain, extraction evidence, reachable
local bounds with provenance, error after each round, reconstruction contributions
and fixed candidate failure assumptions. Primary failure events use the existing
ID-deduplicated union bound; derived errors are not counted as new events.

## Validation

Commands run:

```
cmake --build build-analysis -j
ctest --test-dir build-analysis -R 'test_evalround_(planner|reference)$' --output-on-failure
cmake -S . -B build-pr0 -DBUILD_TESTING=ON -DM2424_ENABLE_EVALMOD_ANALYSIS=OFF
cmake --build build-pr0 -j
ctest --test-dir build-pr0 --output-on-failure
ctest --test-dir build-analysis --output-on-failure
```

New tests separately: 2/2 PASS. Baseline: 12/12 PASS. Analysis: 18/18 PASS.
The reference test is intentionally optional/analysis-only; the planner test is
part of both suites. `otool -L build-pr0/bin/test_evalround_planner` lists only
system C++/system libraries, not GMP or MPFR.

Tests cover exhaustive integer centers and both interval endpoints for several
small K, exact zero-error reconstruction, interior and out-of-domain inputs,
exact rational stress tests for both cleaning inequalities (including local
error), MPFR phase/cleaning diagnostics against analytical bounds, all rejection
gates, reconstruction weights, exhaustive verification of minimum cleaning counts,
nonuniform per-digit counts, and cost-driven binary/ternary selection at K=1.

Representative cases (rho=0.01 and requiredIntegerError=1e-8 for K=1):

| Case | Selected radix | Total cleaning iterations | Certified reference E_I | Cost |
|---|---:|---:|---:|---:|
| Binary extraction cost 1, ternary cost 10 | 2 | 6 | 2.492669878e-9 | 1.8 |
| Binary extraction cost 10, ternary cost 1 | 3 | 3 | 9.349471465e-11 | 1.5 |
| K=8, rho=0.2, piecewise target, ternary lower cost | 3 | 0 | 0 (exact target) | 2 |

For requiredIntegerError=1e-10, the reference test selects eight total binary
cleanings with E_I=2.208425134e-17, or three ternary cleanings with
E_I=9.349471465e-11. Every tested MPFR stage stays below its analytical bound;
measurements are not used to derive that bound.

No existing tolerance, security parameter, old EvalMod implementation, PR-1
HP/LP computation or unresolved PR-1 backend bound was changed. No ciphertext
levels or modulus bits are consumed by PR-2.

## Remaining limitations

No blocker remains for this reference/planner scope. Polynomial approximation
certificates for general practical extraction/DigitExtract remain external;
unknown candidates are retained only as diagnostics/rejections. Ciphertext
execution and its nonzero local arithmetic bounds, headroom, exact scale/modulus
schedule and key requirements are deliberately not supplied in PR-2. Thus a
Certified reference plan is not a production-executable bootstrap plan.
