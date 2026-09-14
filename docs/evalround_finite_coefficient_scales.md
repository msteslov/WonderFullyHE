# PR-6 follow-up: finite coefficient plaintext scales

Base revision: `6e3d18405f2b6abd056b30c9ab05d9e5cfc83b1d`.

This follow-up changes only the backend encoding/scheduling of an already
selected exact EvalRound polynomial. It does not change any coefficient string,
whole-domain interval certificate, candidate-selection rule, polynomial power
DAG, `K=64`, rho, target, security setting, tolerance, sparse CtS semantics, or
legacy EvalMod Remez vector.

## Why the common denominator stopped compilation

The original generic power-basis path forms the least common denominator `D`
of all nonconstant rational coefficients, evaluates the integer-numerator
polynomial, and finally performs a certified `MultiplyPlain(1/D)` at plaintext
scale `D`. This is exact and remains the preferred path. It is still used when
all required runtime scales are finite and satisfy the existing scale,
plaintext-in-modulus, and exact centered-headroom gates.

The selected K=64 coefficients are exact terminating decimal strings with a
common denominator larger than finite binary64 range. Although the mathematical
polynomial is finite and certified, `D` cannot be a real SEAL plaintext scale.
That representation restriction caused the old
`ScaleScheduleInfeasible: Exact polynomial common denominator has no finite
binary64 plaintext scale` result.

## Direct finite-scale fallback

If the common-denominator path cannot obtain a finite valid runtime plaintext
scale, each nonzero coefficient is compiled directly. For the unchanged exact
mathematical coefficient `c` and selected finite positive runtime scale `S`, the
backend encodes

    n = round(c*S)
    c_hat = n/S
    delta_c = |c_hat-c|.

All three values are retained exactly in the node trace: the ideal coefficient
`c`, integer `n`, represented rational `n/S`, and rational `delta_c`. A
coefficient that rounds to zero remains an explicit `MultiplyPlain` node and
charges its full `|c|` error. The constant term uses the same existing
`AddPlain` encoding semantics.

For an input state `(M,E)`, the coefficient-encoding contribution is the
existing certified Builder term

    B_enc = (M+E)*delta_c.

The ideal propagation remains `M' = M*|c|`, `E_prop = E*|c|`. Exact runtime
scale-ratio error is then added as before. Thus the arithmetic certificate is
against the original exact polynomial, never against a silently substituted
`c_hat` polynomial.

## Deterministic bounded scale search

The fallback considers only exactly representable powers of two. It starts at
the largest power of two strictly below the smallest actual active SEAL prime
and searches downward, with an explicit cap of 2048 candidates and a lower
bound at the minimum normal binary64 power of two. This is a baseline schedule,
not Grafting or another PR-7 optimizer.

Every candidate is checked using the actual node state and exact arithmetic:

- `S` and the binary64 output ciphertext scale are finite and positive;
- the output scale is below the exact product `Q` of active primes;
- `n=round(c*S)` fits the exact active modulus;
- coefficient error and scale-representation error are included in `(M',E')`;
- `Q - 2*ceil(outputScale*(M'+E'))` is positive.

Builder construction repeats the normal hard gates, so the scheduler cannot
bypass execution-DAG validation. Scale alignment remains real
`MultiplyPlain(1)` and `ModSwitch` arithmetic. No operation or executor path
rewrites ciphertext scale metadata.

## Regression and K=64 result

K=1 continues to use the exact common-denominator path. Its reachable DAG
remains 54 nodes with the same exact prime/scale schedule, required
relinearization and conjugation keys, cleaner structure, and certified integer
error. The ciphertext tests continue to require every observed node and final
integer error to be at or below its certificate.

For the unchanged K=64 candidate, the analysis fixture remains `N=16`, 48
actual 50-bit coefficient primes, input scale `2^49`, and both required keys.
The mathematical `1e-2` plan is unchanged:

    rounds = [0,8,3,4,4,3,2,0]
    mathematical E_I = 9.3076657406671488e-3

The common-denominator failure is no longer the first blocker. Direct finite
coefficient encoding reaches the next concrete certified gate:

    ErrorBudgetExceeded:
    Extraction digit 0 exact arithmetic error exceeds cleaner domain

The exact coefficient-encoding and propagated power-basis error is larger than
the mandatory `a<=1` cleaner-domain gate. It is not reported as unavailable,
clamped, or replaced by an empirical tolerance. The failure occurs before an
immutable reachable plan is published (`nodes=0`), so K=64 is not executed and
no reachable-DAG depth, level, rescale, scale-range, headroom minimum, backend
`E_I`, or observed error is claimed. Work on this newly exposed blocker is
outside this focused follow-up.
