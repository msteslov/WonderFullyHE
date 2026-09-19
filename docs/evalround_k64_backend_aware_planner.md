# K=64 backend-aware bounded EvalRound planner

## Scope and fixed inputs

This planner closes only the current binary direct-polynomial candidate space.
It consumes the 64 records already produced by the bounded search: eight digits,
the `MultiIntervalChebyshev` and `MultiIntervalMinimax` families, and degrees
64, 128, 192, and 256. It does not synthesize another polynomial. The target
table, canonical decimal coefficient strings, direct-x and centered interval
certificates, `K=64`, rho bits `0x3f0b460edc2fc0f3`, and required integer error
`1e-2` are unchanged.

Every record is checked before compilation. Its scaled-Chebyshev polynomial
`q(t)` must satisfy

    convertScaledChebyshevToMonomial(q, "64") == p

coefficient by coefficient in exact rational arithmetic, where `p` is the
canonical certified polynomial. A sampled or binary64-tolerance comparison is
never used for this gate.

## Why approximation-only selection was wrong

The previous path selected one record per digit using `E_approx`, then compiled
only those eight records. Backend error can change the ordering. In the fixed
fixture, digit 7 demonstrates this directly:

| selection rule | family | degree | `E_approx` | `E_backend` | `a0` |
|---|---|---:|---:|---:|---:|
| approximation only | MultiIntervalMinimax | 256 | 2.8084032784527216e-7 | not used for selection | not used for selection |
| backend-aware | MultiIntervalMinimax | 128 | 5.6873231029276063e-4 | 9.4593436569510261e-2 | 9.5162168879803027e-2 |

The degree-128 record has a worse approximation certificate but a smaller
total backend-certified initial error. The planner therefore compares
`a0=E_approx+E_backend`, never `E_approx` alone.

## Actual compilation and trajectory rules

All 64 extraction DAGs use one real projection, the real
`MultiplyPlain(1/64)` normalization, and a shared memoized Chebyshev basis when
runtime state is compatible. Every coefficient encoding, level alignment,
scale alignment, multiplication, relinearization, rescale, modulus switch, and
plaintext operation is a real exact-arithmetic Builder transition. There is no
metadata-only scale normalization.

For each executable extraction, the planner records `E_approx`, the actual DAG
error `E_backend`, and `a0`. It retains round zero as an option when `a0<=1`.
Before every cleaner application it requires `a<=1`, compiles the actual
baseline `f2` block, computes its local error `B_cln`, and applies exactly

    a_next <= 5*a^2 + B_cln.

A result above one is a `CleaningDomainViolation`. A result with
`a_next>=a` is deterministically non-contractive and exhausts that trajectory;
there is no future contraction mechanism in this baseline. No later round is
compiled after either condition. Backend failures of optional later rounds are
stored on that candidate only.

The bounded state includes the candidate identity, cleaning count, certified
error, output level and scale, exact centered headroom, and deterministic
resource counts. If reconstruction can still be feasible, dominance is valid
only when a state is no worse in error, critical levels, scale feasibility,
headroom, ciphertext multiplications, rescales, and reachable nodes. Feasible
plans use lexicographic cost: critical levels, ciphertext multiplications,
rescales, reachable nodes, then candidate identity. Optimality is only within
the enumerated candidate space.

## Bounded negative result

In the fixed `N=16`, 48 by 50-bit-prime, input-scale `2^49` fixture, the best
executable candidate by total `a0` and its reachable minimum are:

| digit | family | degree | `E_approx` | `E_backend` | `a0` | reachable rounds | minimum digit error | first trajectory stop |
|---:|---|---:|---:|---:|---:|---:|---:|---|
| 0 | MultiIntervalChebyshev | 256 | 6.6771536303215967e-9 | 7.7465649189410929e-2 | 7.7465655866564564e-2 | 6 | 3.3832507614222893e-12 | optional round 6 headroom |
| 1 | MultiIntervalChebyshev | 128 | 1.9387252977226169e-1 | 1.8146416825708209e-1 | 3.7533669802934377e-1 | 0 | 3.7533669802934377e-1 | non-contractive round 0 |
| 2 | MultiIntervalChebyshev | 256 | 7.7493658741219940e-2 | 1.9863877514609588e-1 | 2.7613243388731584e-1 | 0 | 2.7613243388731584e-1 | non-contractive round 0 |
| 3 | MultiIntervalChebyshev | 256 | 1.2382838484582197e-1 | 3.4421609355339466e-1 | 4.6804447839921665e-1 | 0 | 4.6804447839921665e-1 | round 0 leaves `a<=1` |
| 4 | MultiIntervalChebyshev | 256 | 1.1161747729148501e-1 | 2.8296923804282936e-1 | 3.9458671533431439e-1 | 0 | 3.9458671533431439e-1 | non-contractive round 0 |
| 5 | MultiIntervalChebyshev | 256 | 7.2707326724283164e-2 | 1.9144396504228156e-1 | 2.6415129176656471e-1 | 0 | 2.6415129176656471e-1 | non-contractive round 0 |
| 6 | MultiIntervalChebyshev | 192 | 2.9110040933540242e-2 | 2.3219823641549006e-1 | 2.6130827734903028e-1 | 0 | 2.6130827734903028e-1 | non-contractive round 0 |
| 7 | MultiIntervalMinimax | 128 | 5.6873231029276063e-4 | 9.4593436569510261e-2 | 9.5162168879803027e-2 | 6 | 3.3832084563092870e-12 | optional round 6 headroom |

Even before reconstruction arithmetic, the exact weighted lower bound

    sum_j 2^j * min_error_j

is approximately `37.089517491455155`, far above `1e-2`. Reconstruction-local
error is nonnegative, so no Cartesian choice can satisfy the required integer
error. The global failure is therefore `ErrorBudgetExceeded` from cleaner
non-convergence/reconstruction budget, not the optional digit-0 or digit-7
headroom failures. No immutable K=64 DAG is published.

The diagnostic construction covers 32,626 nodes, 413 ciphertext
multiplications, 413 relinearizations, 680 rescales, 9,607 modulus switches,
and 10,936 plaintext multiplications. Maximum constructed multiplication depth
is 20, maximum level consumption is 32, and runtime scales range from
`140737488359043.38` to `2.5108406942623581e58`.

This is NOT a global impossibility theorem for v10 EvalRound. It only closes
the current binary direct-polynomial candidate space.

## K=1 certificate drift

K=1 remains the same 54-node common-denominator monomial DAG, with no
scaled-Chebyshev execution and no metadata scale normalization. The historical
`1.23363394493e-5` value belongs to the analytic binary-quadratic extraction
certificate. The current `1.2331869627679e-5` value belongs to the later
whole-domain outward-interval certificate for the same exact polynomial.
Replacing the conservative analytic extraction bound by the independently
recomputed interval bound tightens the composed result; it does not change the
polynomial, DAG, cleaner, runtime arithmetic, target, tolerance, or parameters.
The tighter deterministic value and the 54-node count are locked by regression
tests, and observed errors must remain below it with exact runtime metadata
agreement.
