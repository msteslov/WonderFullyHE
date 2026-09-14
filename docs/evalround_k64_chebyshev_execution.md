# K=64 scaled-Chebyshev execution follow-up

Base revision: `7d8304403c2555387bb0b426cb15fd4214fd8af0`.

This PR-6 follow-up changes only the ciphertext execution representation of the
already selected exact K=64 extraction polynomials. It does not change any
canonical coefficient string, target, degree, family selection, whole-domain
certificate, `K`, rho bits, cleaner polynomial, security input, tolerance, or
sparse `K=h_b=64` contract.

## Certified polynomial and execution representation

The whole-domain certificate remains attached to the canonical exact monomial
polynomial `p_j(x)`. Direct evaluation of that representation forms large
intermediate powers and coefficients whose cancellation is invisible to
triangle/product magnitude propagation. Its former `a_0>1` result was therefore
a failure of that execution DAG, not a failure of the certified polynomial.

Generation now also retains the natural exact execution representation

    q_j(t) = sum_k a[j,k] T_k(t),    t = x/64,

including the pre-conversion Chebyshev result from the existing Remez run for
digit 7. Before any ciphertext arithmetic, compilation parses the decimal
coefficients and variable scale as exact rationals, converts `q_j(x/64)` to
monomial form, trims only exact trailing zeros, and requires coefficient-wise
equality with canonical `p_j(x)`. A changed coefficient or variable scale
returns `ExtractionNotCertified`. All eight selected K=64 representations pass
this algebraic gate. The canonical centered certificates are independently
recomputed against `p_j` as before.

Normalization is a real `MultiplyPlain(1/64)` at plaintext scale 64, for which
the encoded integer is exactly one. There is no scale metadata rewrite and no
normalization-only rescale.

## Fast-doubling DAG and magnitude certificate

The compiler memoizes `T_k` and uses

    T_(2m)   = 2*T_m^2 - 1,
    T_(2m+1) = 2*T_m*T_(m+1) - T_1.

Every ciphertext product is followed by the baseline relinearize/rescale
sequence. Constants, coefficient products, level alignment, and scale alignment
are real Builder operations. Unit tests establish the exact `T_2`, `T_3`, and
`T_8` identities, compare `T_255` and `T_256` with an independent exact
recurrence, and bound their critical multiplication depth by 8 rather than a
256-step chain.

For the exact binary64 rho from the problem,

    R = (64 + rho)/64 >= 1.

The ideal magnitude of every completed `T_k` node is tightened, without
tightening its semantic error, to the exact rational recurrence

    B_0=1, B_1=R, B_(k+1)=2*R*B_k-B_(k-1).

Thus `|T_k(t)| <= T_k(R)` on `[-R,R]`. Tests cover both endpoints through
degree 256. Exact centered headroom is recomputed after every tightening.

## Backend error composition

Finite coefficient encoding retains exact `a_k`, runtime
`round(a_k*S)/S`, and the exact encoding delta. Zero-rounded nonzero
coefficients are not removed. The Builder charges input propagation,
multiplication interaction, key switching, divide-round, constant encoding,
binary64 scale representation, real level alignment, and real scale alignment
on the actual fast-doubling graph.

At extraction boundaries the implementation records separately
`E_approx`, the actual Chebyshev-DAG `E_backend`, and
`a_0=E_approx+E_backend`. It then runs the unchanged real binary cleaner and
uses `a_next <= 5*a^2+B_cln`. Exact errors are rounded only upward to a
384-bit dyadic certificate at stable stage boundaries, avoiding unbounded GMP
expression growth without weakening any gate. Reconstruction reservation
replays the same exact scalar/add transitions without copying the complete
published Builder trace for every Cartesian cleaning choice.

## K=64 result in the fixed analysis fixture

The fixture is unchanged: `N=16`, 48 50-bit coefficient primes, input scale
`2^49`, and present relinearization/conjugation keys. Every extraction passes
the exact cleaner-domain gate:

| digit | selected family | degree | `E_approx` | `E_backend` | `a_0` | real cleaner rounds reached |
|---:|---|---:|---:|---:|---:|---:|
| 0 | MultiIntervalChebyshev | 256 | 6.6771536303215967e-9 | 7.7465649189410929e-2 | 7.7465655866564564e-2 | 6 |
| 1 | MultiIntervalChebyshev | 128 | 1.9387252977226169e-1 | 1.8146416825708209e-1 | 3.7533669802934377e-1 | 2 |
| 2 | MultiIntervalChebyshev | 256 | 7.7493658741219940e-2 | 1.9863877514609588e-1 | 2.7613243388731584e-1 | 3 |
| 3 | MultiIntervalChebyshev | 256 | 1.2382838484582197e-1 | 3.4421609355339466e-1 | 4.6804447839921665e-1 | 1 |
| 4 | MultiIntervalChebyshev | 256 | 1.1161747729148501e-1 | 2.8296923804282936e-1 | 3.9458671533431439e-1 | 2 |
| 5 | MultiIntervalChebyshev | 256 | 7.2707326724283164e-2 | 1.9144396504228156e-1 | 2.6415129176656471e-1 | 3 |
| 6 | MultiIntervalChebyshev | 192 | 2.9110040933540242e-2 | 2.3219823641549006e-1 | 2.6130827734903028e-1 | 3 |
| 7 | MultiIntervalMinimax | 256 | 2.8084032784527216e-7 | 1.1313778939100920e-1 | 1.1313807023133705e-1 | 6 |

These values supersede the cancellation-heavy monomial backend errors but do
not replace the unchanged approximation bounds. A real cleaning schedule is
not selected: the first concrete later gate is

    HeadroomViolation:
    cleaning digit 0 round 6: centered no-wrap proof unavailable.

Consequently no complete immutable K=64 execution DAG is published and, per
the execution contract, K64 is not decrypted or empirically claimed. The
bounded diagnostic construction before fail-closed publication contains 7697
nodes, 307 ciphertext multiplications, 307 relinearizations, 548 rescales,
1949 ModSwitches, and 2446 plaintext multiplications. Its maximum constructed
multiplication depth is 20, maximum level consumption is 32, and runtime scales
range from `140737488359463.38` to `2.5108406942623581e58`. Every published
intermediate headroom before the failing operation is positive.

## K=1 regression

K=1 has no execution-representation metadata and remains on the exact
common-denominator monomial compiler. Its reachable immutable DAG is still 54
nodes. The certified integer error is `1.2331869627679e-5`; the regression
fixture observes at most `1.6810572769366e-6`, with every node checked against
its certified error, modulus, chain index, scale bits, and component count.

The next production-enablement step is therefore the concrete K64 cleaner
headroom/schedule gate above, not coefficient regeneration or a new
approximation method.
