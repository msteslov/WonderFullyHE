# K=64 exact centered-certificate follow-up

Base revision: `8c8033c77682cd81d05025408d98d959ab018621`.
The hard mathematical contract for this follow-up is the repository-designated
v9 EvalRound model plus `codex_bootstrap_implementation_spec_v1.md`. No v10
text changes that contract. Candidate generation, `K=h_b=64`, rho bits
`0x3f0b460edc2fc0f3`, cleaner, executor, security parameters, tolerances, and
production modulus chain are unchanged.

## Exact centered proof

For the existing exact-decimal monomial polynomial

    p(x) = sum_{n=0}^d c_n x^n

the verifier parses every `c_n` as an exact `mpq_class`. At each integer center
`I`, it computes the exact coefficients of the same polynomial under `x=I+y`:

    d[I,k] = sum_{n=k}^d c_n * binom(n,k) * I^(n-k)
    p(I+y) = sum_{k=0}^d d[I,k] y^k.

The implementation performs this Taylor shift by exact rational Horner
composition with `(I+y)`, which is algebraically identical to the displayed
formula and avoids any binary64 coefficient operation. It subtracts the exact
integer target `bit_j(I+K)` from `d[I,0]` and evaluates the result using
384-bit outward MPFR interval Horner arithmetic on a complete closed partition
of the exact dyadic interval `[-rho,+rho]`.

The original direct-x verifier remains present and independently evaluates the
same exact coefficients and target on every closed subcell of
`[I-rho,I+rho]`. A digit certificate stores `E_direct_x`, `E_centered`, the
subdivision count, precision, and selected proof method. Its production bound
is exactly `min(E_direct_x,E_centered)`. Grid error is never part of that
minimum. The execution compiler recomputes both bounds from the supplied exact
coefficients and rejects stale method, precision, K, rho, digit, target,
coverage, or bound metadata.

## Selected K=64 candidates

The Remez and cosine-diagnostic candidates were not regenerated or changed.
All 64 existing family/digit/degree records were recertified. The report
executable prints grid, direct-x, centered, selected, ratio, and proof method
for every record. Selection below uses only the rigorous minimum.

| digit | family | degree | grid (diagnostic) | old direct-x | centered | selected rigorous | direct/centered | `a_j<=1` |
|---:|---|---:|---:|---:|---:|---:|---:|:---:|
| 0 | MultiIntervalChebyshev | 256 | 6.6771538298482369e-9 | 6.7951851457391055e77 | 6.6771536303215967e-9 | 6.6771536303215967e-9 | 1.0176769207288442e86 | yes |
| 1 | MultiIntervalChebyshev | 128 | 1.9387252977225700e-1 | 9.5739985816692474e38 | 1.9387252977226169e-1 | 1.9387252977226169e-1 | 4.9382955867526146e39 | yes |
| 2 | MultiIntervalChebyshev | 256 | 7.7493658741220106e-2 | 2.1071195281809190e88 | 7.7493658741219940e-2 | 7.7493658741219940e-2 | 2.7190863903037181e89 | yes |
| 3 | MultiIntervalChebyshev | 256 | 1.2382838484581304e-1 | 1.3331939453866094e88 | 1.2382838484582197e-1 | 1.2382838484582197e-1 | 1.0766464789527553e89 | yes |
| 4 | MultiIntervalChebyshev | 256 | 1.1161747729148175e-1 | 9.7259730140792264e87 | 1.1161747729148501e-1 | 1.1161747729148501e-1 | 8.7136649654607398e88 | yes |
| 5 | MultiIntervalChebyshev | 256 | 7.2707326724283372e-2 | 1.6474402244707921e88 | 7.2707326724283164e-2 | 7.2707326724283164e-2 | 2.2658517355728492e89 | yes |
| 6 | MultiIntervalChebyshev | 192 | 2.9110040933535808e-2 | 4.9526575968336431e63 | 2.9110040933540242e-2 | 2.9110040933540242e-2 | 1.7013571393255068e65 | yes |
| 7 | MultiIntervalMinimax | 256 | 2.8084032784527216e-7 | 2.9417790095985345e41 | 2.8084032784527216e-7 | 2.8084032784527216e-7 | 1.0474916591107584e48 | yes |

All eight selected rigorous bounds satisfy the cleaner domain `a_j<=1`.
Except for digit 1, they are already below the stricter zero-local-error
contraction threshold `1/5`; digit 1 is also below `1/5` and contracts slowly.

## Unchanged cleaning and reconstruction planner

The initial vector for every budget is:

    [6.6771536303215967e-9, 1.9387252977226169e-1,
     7.7493658741219940e-2, 1.2382838484582197e-1,
     1.1161747729148501e-1, 7.2707326724283164e-2,
     2.9110040933540242e-2, 2.8084032784527216e-7]

The unchanged exact-reference recurrence is `a_next=5a^2` here because these
diagnostic mathematical plans have zero backend-local cleaner error. Each
sequence below includes the initial extraction error; `weighted` is
`2^j*a_final`.

### Budget `1e-2`

| j | rounds | certified sequence | weighted |
|---:|---:|---|---:|
| 0 | 0 | 6.6771536303e-9 | 6.6771536303e-9 |
| 1 | 8 | .19387253, .18793279, .17659367, .15592661, .12156555, .07389091, .02729933, .00372627, 6.9425352e-5 | 1.3885070464e-4 |
| 2 | 3 | .07749366, .03002634, .00450790, 1.0160600e-4 | 4.0642400291e-4 |
| 3 | 4 | .12382838, .07666734, .02938941, .00431869, 9.3255273e-5 | 7.4604218257e-4 |
| 4 | 4 | .11161748, .06229231, .01940166, .00188212, 1.7711906e-5 | 2.8339050161e-4 |
| 5 | 3 | .07270733, .02643178, .00349319, 6.1012026e-5 | 1.9523848292e-3 |
| 6 | 2 | .02911004, .00423697, 8.9759676e-5 | 5.7446192806e-3 |
| 7 | 0 | 2.8084032785e-7 | 3.5947561964e-5 |

Certified `E_I = 9.3076657406671488e-3`.

### Budget `1e-4`

Rounds are `[0,9,4,5,5,4,3,0]`. The complete per-digit sequences are:

    j0: 6.6771536303e-9
    j1: .19387253, .18793279, .17659367, .15592661, .12156555,
        .07389091, .02729933, .00372627, 6.9425352e-5, 2.4099398e-8
    j2: .07749366, .03002634, .00450790, 1.0160600e-4, 5.1618897e-8
    j3: .12382838, .07666734, .02938941, .00431869, 9.3255273e-5,
        4.3482730e-8
    j4: .11161748, .06229231, .01940166, .00188212, 1.7711906e-5,
        1.5685581e-9
    j5: .07270733, .02643178, .00349319, 6.1012026e-5, 1.8612337e-8
    j6: .02911004, .00423697, 8.9759676e-5, 4.0283997e-8
    j7: 2.8084032785e-7

Weighted contributions are
`[6.6771536303e-9, 4.8198795445e-8, 2.0647558768e-7,
3.4786183636e-7, 2.5096930126e-8, 5.9559476895e-7,
2.5781758343e-6, 3.5947561964e-5]`. Certified
`E_I = 3.9755642870705469e-5`.

### Budget `1e-6`

Rounds are `[0,9,4,5,5,5,4,1]`. Relative to the `1e-4` schedule, digit 5 adds
`1.7320953555e-15`, digit 6 adds `8.1140022371e-15`, and digit 7 adds
`3.9435644872e-13`; the complete sequences for digits 0--4 are exactly those
listed for `1e-4`, while the changed complete sequences are:

    j5: .07270733, .02643178, .00349319, 6.1012026e-5, 1.8612337e-8,
        1.7320954e-15
    j6: .02911004, .00423697, 8.9759676e-5, 4.0283997e-8, 8.1140022e-15
    j7: 2.8084032785e-7, 3.9435644872e-13

Weighted contributions are
`[6.6771536303e-9, 4.8198795445e-8, 2.0647558768e-7,
3.4786183636e-7, 2.5096930126e-8, 5.5427051376e-14,
5.1929614318e-13, 5.0477625436e-11]`. Certified
`E_I = 6.3436135559179875e-7`.

## Backend feasibility

Because the mathematical plan is Certified, the shallowest `1e-2` schedule
was supplied to the unchanged generic ciphertext compiler. The analysis-only
test context deliberately provides 48 50-bit coefficient primes at `N=16`,
scale `2^49`, and both relinearization and conjugation keys. Security is
disabled only for this test fixture; the production chain is untouched.

Compilation stops with `RequiredBoundUnavailable: Arithmetic bound overflow`.
This occurs before an immutable reachable DAG is produced (`nodes=0`), despite
the deliberately ample level count. Therefore extraction nodes/depth, cleaning
depth, critical-path levels, rescale count, reconstruction operations, scale
range, and headroom are not available and are not claimed. The exact first
remaining blocker is the existing compiler arithmetic-bound representation,
not the centered polynomial certificate, cleaner domain, planner budget, or
available analysis levels.

## Regression coverage

Exact-shift tests cover degrees 0, 1, 2, 3, a sparse degree-128 polynomial,
positive and negative centers, exact rational evaluation points, both complete
interval endpoints, and the existing K=1 b0/b1 polynomials. K=1 selected
certificates are constrained never to exceed their independent direct-x bound.
Negative coverage rejects grid-only evidence, wrong K/rho/digit/target,
modified and non-exact coefficients, stale centered metadata, incomplete
coverage, and a missing digit. The locked legacy EvalMod degree-7 Remez vector
and unchanged K=64 generator tests remain in place.
