> The [generic polynomial follow-up](evalround_generic_polynomial.md) replaces
> the original hard-coded extraction compiler while preserving this K=1
> regression schedule and certificate. This document records the original PR-3.

# EvalRound ciphertext baseline (PR-3)

`EvalRoundPlanStatus::Certified` remains a mathematical/reference certificate.
It is not an execution capability. The core executor accepts only an immutable
`EvalRoundExecutionPlan`, produced by the optional analysis compiler. There is no
conversion from the mutable reference plan and no metadata scale rewrite opcode.
No full bootstrap certificate is produced.

The first executable candidate is `BinaryQuadraticK1/SEAL-baseline`. Its input
contract is a two-component ciphertext representing a real z in D_(1,rho), with
a separately supplied deterministic semantic error bound. Both coefficient halves
can use the same plan through `executeEvalRoundPair`. Both inputs are preflighted
before either half executes. Domain membership and the upstream error bound are
conditional assumptions, not inferred by decrypting the input.

## Polynomial and scale schedule

A real projection `(x+conjugate(x))/2` preserves the ideal real input and includes
conjugation/key-switching noise. Extraction implements exactly
`b0=1-x*x`, `b1=(x*x+x)/2`, with the shifted binary reconstruction `b0+2*b1-1`.
To add x to its square, a real plaintext multiplication by 1 raises its scale to
the square's scale. Half constants use plaintext scale 2 and exact integer 1.

Each cleaner evaluates `a*a`, relinearizes and rescales, then multiplies by
`3-2*a`, relinearizes and rescales again. Its linear operand is modulus-switched
to the square's level. This is the baseline two-rescale circuit, without deferred
tensor products or Thrifty evaluation.

Reconstruction aligns levels and multiplies the two operands by plaintext 1 and
2 encoded at the other operand's actual scale. Both scale products therefore have
identical binary64 bits. This is an arithmetic operation with an encoding and scale
representation bound; it is not a metadata correction. There is no extra rescale
for this alignment. The resulting high scale must pass centered headroom.

Every node records its operation, parents, key requirement, actual SEAL chain
index, exact active primes, input and output dyadic scales, exact arithmetic scale
before binary64 rounding, magnitude, propagated/local/total semantic errors,
constant encoding error, scale representation error, and exact positive centered
headroom with provenance. `binary64Bits` on an arithmetic scale identifies its
corresponding rounded runtime representation. Index zero is the bottom SEAL data
level, not the number of levels consumed.

## Arithmetic proof

All compiler bound calculations use exact GMP rationals, with outward conversion
to binary64 after each node. These dependencies remain confined to
`m2424_evalmod_analysis`; the core immutable plan and executor do not require GMP
or MPFR. Plaintext scalar integers are rounded from exact rational coefficients
and dyadic scales, then encoded directly as exact residues at the correct level.
The constant error is `abs(encodedInteger/encodingScale - idealConstant)`.

For addition, errors add. For multiplication, propagated error is
`M_a E_b + M_b E_a + E_a E_b`. RNS addition and ring multiplication themselves are
exact under the separately checked no-wrap gate. If an operation's exact scale
and its runtime dyadic scale differ by ratio r, the certificate adds
`abs(r-1)*(M+propagatedError+localError)`. This includes conversion of a rescale
prime to binary64, as well as multiplication/division rounding. No actual prime
product is approximated in double.

For N coefficients, ternary secret support 1 and t ciphertext components, a
componentwise divide-and-round contributes at most
`N*(1+N+...+N^(t-1))/(2*outputScale)` in slot norm. With CBD evaluation-key noise
support 21, key switching adds at most
`N^2*21*sum(q_i-1)/(P*scale)` plus its two-component ModDown rounding bound.
Relinearization multiplies this bound by the number of eliminated components.
Conjugation uses one switch. These are deterministic support bounds, not observed
noise estimates. Unknown or smaller-than-backend support claims are rejected.

At every node the compiler checks the exact positive margin
`Q/2-ceil(scale*(M+E))`, including the destination of every modulus switch.
The runtime verifies actual scales, primes and chain indices after every operation.
Preflight checks the initial context, exact scale, chain, component count, rounding
mode and required keys before arithmetic starts.

The per-node semantic error is relative to the exact polynomial DAG on the ideal
input. The final integer error is a separate certificate using the contracting
v9 recurrence. To compute a cleaner's local B, its block is re-evaluated with zero
incoming *arithmetic* error and the bound `|actual digit|<=1+a`. All internal
arithmetic remains included. PR-2 then receives extraction error including backend
arithmetic, these B values and conservative reconstruction errors. It selects the
minimum cleaning count for the caller's `requiredIntegerError`. The final bound is
exactly the PR-2 weighted reconstruction bound. A larger naive propagated DAG error
must not be mistaken for the contracting integer-error certificate, or vice versa.

## Scope and remaining blockers

* Executable extraction is limited to binary K=1. Piecewise reference and exact
  ternary phase targets are rejected. A ternary backend requires a separately
  certified whole-domain polynomial approximation and is not implemented here.
* The compiler searches the requested PR-2 cleaning range up to feasible baseline
  levels/scales. Reconstruction local bounds are uniform over retained counts;
  this conservative choice can reject budgets that another arithmetic schedule
  could satisfy. This is not a new schedule/factorization optimizer.
* A rigorous upstream input bound remains required. Unknown CoeffToSlot backend
  bounds from PR-1 are not solved or replaced. Probabilistic upstream certificates
  and unsupported Gaussian key-noise builds are rejected by this first path.
* Matching ciphertext/evaluation secret keys and the supplied domain/error
  contract are caller responsibilities. The plan checks key availability and
  context parameters; it does not decrypt to establish these assumptions.
* PR-4 and subsequent bootstrap stages are absent.

Tests use exact scalar RNS inputs at all centers and both interval endpoints,
an independent exact-rational polynomial oracle with 256-bit MPFR diagnostics,
and test-only decryption of every node and every CKKS slot. The upstream test
bound is derived from ternary/CBD encryption support and encryption ModDown
rounding. Reference observations never become production bounds.

## Recorded validation

The baseline suite passed 12/12 and the optional analysis suite passed 20/20,
including both new tests. Existing PR-0/PR-1/PR-2 test sources and tolerances were
unchanged. Commands: `cmake --build build-pr0 -j 4`,
`ctest --test-dir build-pr0 --output-on-failure`,
`cmake --build build-analysis -j 4`, and
`ctest --test-dir build-analysis --output-on-failure`.

For N=32768, nine 60-bit data primes plus one 60-bit special prime, input scale
2^59, K=1 and rho=1/128, the selected plan has 54 reachable nodes and two cleaning
iterations per digit. For the looser external budget 0.01 the planner instead
selects one iteration per digit.

| Boundary | Actual SEAL chain index | b0 scale | b1 scale |
| --- | ---: | ---: | ---: |
| Input | 8 | 5.7646075230342349e17 | same input |
| Extraction | 7 | 1.1529215046149734e18 | 2.3058430092299469e18 |
| First cleaner | 5 | 1.1529215046519357e18 | 9.223372037215486e18 |
| Second cleaner | 3 | 1.1529215047667548e18 | 5.9029581044057847e20 |
| Reconstruction | 3 | 6.8056473403066278e38 | combined |

Required integer error: `1e-4`. Certified integer error:
`1.23363394493e-5`. Maximum observed error over all slots at all nine input centers
and interval endpoints in the recorded run: `1.68105737919e-6`. Observations vary
slightly with newly generated keys; they are not used to derive the certificate.
No required arithmetic bound is Unknown in this selected execution plan.

Changed files:

* `include/m2424/evalround_execution.hpp`
* `include/m2424/experimental/evalmod_analysis/evalround_execution.hpp`
* `include/m2424/m2424.hpp`
* `src/core/evalround_execution_internal.hpp`
* `src/ckks/evalround_execution.cpp`
* `src/planning/evalround_execution.cpp`
* `src/CMakeLists.txt`
* `tests/evalround_test_support.hpp`
* `tests/test_evalround_backend.cpp`
* `tests/test_evalround_certificates.cpp`
* `tests/CMakeLists.txt`
* `docs/evalround_ciphertext_baseline.md`
