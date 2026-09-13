# PR-3 follow-up: generic certified polynomial extraction

Base revision: dd2ee4721045053560214e6789d7ce1d89d78277.
See [the pre-edit contract audit](evalround_polynomial_contract_audit.md).
No general-K extraction construction, full sparse Bootstrapper execution change,
security parameter change, tail bound, optimization or PR-7 work is included.

## Three separate contracts

A. **Mathematical/reference extraction.** The integer center I determines the
   offset digits of J=I+K on each interval [I-rho,I+rho]. PR-2's reference plan
   expresses mathematical error/cleaning assumptions. PiecewiseReference and
   exact ternary phase remain reference-only and cannot execute as ciphertexts.
B. **Polynomial approximation.** Each digit carries its actual polynomial and
   a proof for that polynomial minus bit_j(I+K) over the entire union D_K,rho.
   Unknown coefficients, Unknown error, wrong radix/index/target/K/rho, changed
   analytic coefficients, and grid-only evidence fail compilation. The compiler
   rechecks the analytic identity or reruns the outward interval verifier.
C. **Ciphertext arithmetic.** Only an immutable EvalRoundExecutionPlan with the
   actual DAG, keys, primes, binary64 scale schedule and local error/headroom
   bounds may execute. Reference local zeros never become backend local bounds.

A does not imply B, and B does not imply C. No observed error implies any of
these certificates. The generic reference evaluator can evaluate supplied
polynomials without asserting that they approximate digits.

## Shared exact representation

`EvalRoundDigitPolynomial` retains radix, digit index, target, certified K/rho,
whole-domain approximation bound, proof kind, provenance, primary failure-event
IDs and an existing `experimental::EvalModPolynomial`. The latter's storage
(basis and decimal coefficient strings in ascending order) moved into a small
shared dependency-free header; it is the same type, not a second polynomial
format. The old EvalMod exact-decimal parser was moved unchanged into an inline
optional-analysis utility and is reused by both compilers. Decimal strings are
exact rationals, never binary64 coefficient approximations.

The existing K=1 candidate now carries:

| Digit | Basis | Ascending exact decimal coefficients | Polynomial |
|---|---|---|---|
| 0 | Monomial | ["1", "0", "-1"] | 1-x^2 |
| 1 | Monomial | ["0", "0.5", "0.5"] | (x^2+x)/2 |

Its existing all-interval analytic bounds, 2*rho+rho^2 and
(3*rho+rho^2)/2, are rechecked against those coefficients. The identity verifier
is not an arithmetic implementation: all executable extraction goes through
`PolynomialCompiler`, including BinaryQuadraticK1 and ExternalPolynomial.
The independent K=1 reference formula remains available as an oracle.

## Generic DAG and propagation

The compiler constructs powers recursively, scales integer coefficient
numerators, sums monomial terms and applies their exact common denominator.
Constants use Builder's exact RNS encoding. Modulus alignment and scale
alignment use real ModSwitch/MultiplyPlain operations. The fast compatible
scale case retains the old K=1 runtime schedule; a binary64 scale is never
rewritten to fake alignment. Builder charges encoding and scale-representation
errors even for MultiplyPlain(1).

Each actual node receives active primes, exact dyadic scale bits, ideal
magnitude, propagated input error, local arithmetic error and centered headroom.
The generic power-basis baseline accepts nonconstant monomial polynomials of
bounded degree (<=256); unsupported basis, exhausted levels/scales/headroom and
constant-only transparent-ciphertext paths reject. Existing exact Chebyshev
conversion remains an explicit separate utility; no new factorization optimizer
is introduced.

For digit j, the extraction error supplied to the backend planner is

    a_j = verified whole-domain approximation error
        + actual polynomial DAG semantic error.

The latter includes upstream input error through every polynomial operation,
real projection, constant encoding, multiplication, relinearization and rescale.
The two errors are neither substituted for one another nor counted twice.

Cleaning remains the original Builder implementation of f2=x^2*(3-2x), with
hard recurrence a_next<=5*a^2+B_local. Available arithmetic rounds are compiled;
PR-2's dynamic program chooses the minimum total iterations meeting the caller's
requiredIntegerError. There is no fixed round count or automatic 1e-10 target.

Reconstruction is generated for the digit count as sum_j(2^j*b_j)-K. Every
weight, constant shift and scale/modulus alignment is actual certified
arithmetic. The two-digit uniform reservation preserves the old certificate;
the multi-digit baseline replays Builder reconstruction with zero incoming
arithmetic error and full digit envelopes over retained cleaning combinations,
reserving the maximum local error. Search is bounded and can reject excessive
work/scales; it never assumes a general-K polynomial exists. Unreachable search
nodes are removed before execution.

The execution baseline remains deterministic: primary event IDs are retained,
but a polynomial requiring probabilistic events is not promoted to this
finite-support execution certificate.

## Interval evidence reuse

`certifyEvalRoundDigitPolynomial` reuses approximation_lab's MPFR Real/Interval
storage and outward interval multiplication. For every integer I in [-K,K],
it covers the complete closed interval with complete closed cells. Exact
rational endpoints use the supplied binary64 rho without decimal loss.
Outward-rounded coefficient intervals and Horner arithmetic enclose p(x), then
subtract the **constant** offset digit bit_j(I+K). The maximum endpoint absolute
value encloses the error on every cell, including both domain endpoints.

This changes and proves the target explicitly. It does not reuse the old EvalMod
subtraction p(x)-(x-I), its target derivative correction, or its complex target
bound. The old EvalMod verifier's arithmetic is unchanged apart from factoring
out the shared interval primitive. No new mathematical assumption is needed:
v9 already defines the constant digit on each disjoint interval. The verifier
uses 384-bit outward arithmetic and explicit work limits; it does not synthesize
coefficients, infer a certificate from grid samples, or implement DigitExtract.

## K=1 regression and remaining K=64 blocker

The existing secure backend test profile remains N=32768, nine 60-bit data
primes plus one 60-bit special prime, input scale 2^59, rho=1/128 and required
E_I=1e-4. It still uses 54 reachable nodes and two cleanings per digit.
The certified E_I remains 1.23363394493e-5, matching the previous certificate.
Five levels are consumed on the critical path, from chain index 8 to 3.
There are **ten** rescale operations across both digit paths:
two extraction rescales plus two per cleaner iteration, 2+2*(2+2).

| Boundary | Chain index | b0 scale | b1 scale |
|---|---:|---:|---:|
| Input | 8 | 5.7646075230342349e17 | same |
| Extraction | 7 | 1.1529215046149734e18 | 2.3058430092299469e18 |
| First cleaner | 5 | 1.1529215046519357e18 | 9.223372037215486e18 |
| Second cleaner | 3 | 1.1529215047667548e18 | 5.9029581044057847e20 |
| Reconstruction | 3 | 6.8056473403066278e38 | combined |

Maximum observed final error is approximately 1.68105738e-6 over all old backend
centers/endpoints and slots; fresh keys change only observations. Every runtime
node is checked against its certified polynomial semantic error by the unchanged
backend test.

The new tiny-N test independently exercises the external-polynomial interval
route for the existing K=1 polynomials, checks exact schedules and every node,
and evaluates the existing cubic cleaner through generic coefficient machinery
as an arithmetic-only test. Its interval-route certificate is
1.2331869627679e-5, with observed error approximately 1.68105728e-6. This fixture
is not a security claim or a general-K extraction construction.

There is still **no concrete usable whole-domain-certified general-K digit
extractor in v9/repository**. The first blocker toward K=64 is the missing
polynomials and per-digit proof of sufficiently small errors on all 129
intervals for bit_j(I+64), within the cleaner and supplied final-error budgets.
If labelled DigitExtract, the candidate must additionally satisfy its existing
rho/epsilon gate. Merely supplying degrees, a reference target or an EvalMod
polynomial/certificate is insufficient. No fake production K=64 candidate is
introduced. Full sparse Bootstrapper execution and K=h_b are untouched.

## Final validation

Both builds and complete suites passed:

- `cmake --build build-pr0 -j4`; baseline CTest: **15/15**, 15.60 s.
- `cmake --build build-analysis -j4`; analysis CTest: **28/28**, 199.97 s.
- `test_evalround_polynomial` also passed independently.
- `git diff --check`: clean.

The full-suite secure K=1 observation was 1.68105736631e-6, below the
unchanged certificate 1.23363394493e-5 and required 1e-4. All existing test
sources, tolerances and security parameters are unchanged.

Changed files comprise the shared decimal polynomial and exact-parser headers;
EvalRound public/optional extraction and reference headers; polynomial compiler,
execution compiler, reference evaluator and PR-2 candidate factory; the reused
approximation interval implementation and EvalMod parser include; the new
polynomial test and its CMake registration; and this report, the contract audit
and the historical ciphertext-baseline documentation link.
