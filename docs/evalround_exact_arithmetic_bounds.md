# PR-3 exact arithmetic-bound representation follow-up

Base revision: `20494f101e02f57404c57cb2271fc19f0b1cfff6`.
This change follows the v9 EvalRound contract and does not alter extraction
polynomials, candidate generation, the polynomial DAG/evaluation order,
cleaner mathematics, `K=h_b=64`, rho, security parameters, tolerances, or the
production modulus chain.

## Old representation loss

`arithmetic::State` and `Calculation` already computed their formulas in
`mpq_class`, including product propagation, constant encoding, finite-support
key switching/divide-round, and exact scale-ratio error. The old `calculate()`
then replaced the exact state with:

    state.M = q(up(M))
    state.E = q(up(E + local))

`publish()`, `localBlock()`, `tightenMagnitude()`, and stage-bound construction
also required `up(mpq_class)`. Thus a finite rational was converted to an
outward binary64 value and immediately reparsed as a dyadic rational. If it
exceeded `DBL_MAX`, the calculation stopped as
`RequiredBoundUnavailable: Arithmetic bound overflow` before the independent
exact headroom comparison could be authoritative.

For the selected K=64 plan, the first node value that cannot be projected to a
finite binary64 is exactly:

    EvalRound polynomial power node 39: ideal magnitude

The node is not removed or reordered by this change.

## Exact internal state and trace API

`State::M`, `State::E`, and all `Calculation` fields now remain exact
nonnegative `mpq_class` values through propagation. `localBlock()` returns an
exact rational, digit paths retain exact errors, and cleaner recurrence uses
the unchanged exact formula `5*a^2+B_cln` until a PR-2 planner boundary.
Runtime SEAL scales remain their actual binary64 values and enter proofs via
their exact dyadic `mpq_class(scale)` representation.

Each `EvalRoundExecutionNode` now has six `EvalRoundExactBound` fields for:

- ideal magnitude;
- propagated semantic error;
- local arithmetic error;
- total semantic error;
- constant encoding error;
- scale-representation error.

`EvalRoundExactBound` stores exact numerator and denominator, deterministic
kind, provenance, and an optional outward binary64 projection. A finite value
above `DBL_MAX` remains a known deterministic exact bound: its optional
projection is absent, while the compatibility `BootstrapBound` is deterministic
with infinity rather than zero, `DBL_MAX`, or `Unknown`. The exact field is the
authority.

The centered headroom gate remains

    Q/2 - ceil(scale * (M+E))

and is evaluated from exact active-prime product `Q`, exact dyadic runtime
scale, and exact `M,E`, independently of binary64 projection. Tests exercise a
finite `2^1100` magnitude that has no binary64 projection but passes exact
headroom with ample `Q`, and `2^2500` which fails through the concrete
`HeadroomViolation`. A hand-computed multiplication test retains exactly
`M=4` and `E=41/100` from inputs `M=2,E=1/10`.

## Finite planner boundaries

At the extraction and cleaning boundaries, exact errors are first compared
against the cleaner domain and error gates. Only a value that passes the
relevant exact gate is projected outward for the legacy PR-2 planner API.
A known oversized value therefore produces `ErrorBudgetExceeded` rather than
`RequiredBoundUnavailable`. Genuinely absent upstream or key-switch evidence
continues to produce `RequiredBoundUnavailable`.

## K=64 backend result

The unchanged centered certificates select the same polynomials and the
shallowest mathematical plan remains Certified:

    requiredIntegerError = 1e-2
    rounds = [0,8,3,4,4,3,2,0]
    mathematical E_I = 9.3076657406671488e-3

Compilation was attempted with the required analysis-only fixture:

    N = 16
    coefficient primes = 48 x 50 bits
    input scale = 2^49
    relinearization key = present
    conjugation key = present

The artificial arithmetic-bound overflow is gone. Exact propagation and
headroom pass the old node 39 point. The first new real backend gate is:

    ScaleScheduleInfeasible:
    Exact polynomial common denominator has no finite binary64 plaintext scale

This is the actual runtime metadata restriction of the unchanged generic
polynomial evaluation strategy: its exact common denominator cannot be used as
a finite SEAL plaintext scale. It is not an unavailable arithmetic proof.
The failure occurs before an immutable reachable execution DAG is published,
so the public result has `nodes=0`. Consequently reachable depth, consumed
levels, rescale count, minimum reachable-DAG headroom, output scale range, and
backend-certified `E_I` are unavailable and are not claimed. The 48-prime
fixture establishes that ordinary level shortage is not the preceding gate.
Execution/decryption is therefore not permitted for K=64.

## K=1 regression

K=1 retains the same 54-node reachable DAG, exact prime/scale schedule,
required relinearization and conjugation keys, two rescale nodes per cleaning
round, and final backend-certified bound `1.23363394493e-5` for the `1e-4`
budget. Every exact trace value has the same finite outward compatibility
projection as before, and all observed per-node/final errors remain below the
certificate.

Legacy EvalMod Remez coefficients, K=64 Remez/Chebyshev generation and
selection, centered approximation bounds, sparse CtS semantics, and all
security/tolerance inputs are unchanged.
