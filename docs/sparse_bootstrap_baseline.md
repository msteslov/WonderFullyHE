# PR-6: fixed-weight encapsulation and all-family security

Base: bb04c5c05120c164de85698418480b1536a07cf6.
Part A passed; see [the source audit](seal_backend_bound_audit.md). No PR-3/4/5
finite-support formula, tolerance, security parameter or accuracy target changed.

## Executable stage

`SealAdapter::generateSparseBootstrapKeys(h)` generates a uniform h-subset using
partial Fisher–Yates and SEAL's cryptographic PRNG with unbiased bounded draws.
Each selected coefficient has an independent uniform sign. The coefficient
alphabet is {-1,0,1}, with exactly h nonzeros. The API accepts 2<=h<=N; it does
not choose h from the K=1 executor or claim that any accepted h is secure.

Both switch directions use separate identity-automorphism key containers.
`identitySwitchKey` duplicates the **key construction**, not the switching
arithmetic, from `KeyGenerator::generate_one_kswitch_key`: encrypt zero under the
destination secret at Q_top*P and add P times the source secret in decomposition
limb j. `Evaluator::apply_galois_inplace(g=1)` performs the exact identity
permutation and calls the audited `switch_key_inplace`. There is no vendored
SEAL change, ciphertext scale rewrite, rescale, or decryption in this path.

`SparseCipher` and `SparseRaisedCipher` cannot be passed to ordinary arithmetic.
Their key-generation identity is private. `restoreSparse` returns a regular
`RaisedCipher` under the main secret. The prepared stage binds the generation,
weight, context, exact source/raised primes, binary64 scale, and public-material
inventory. Missing either switching key fails preflight. Imported ordinary
keys have unknown distribution and invalidate this generated-key-only path.

The optional analysis library prepares the certificate using the shared
`finiteSupportKeyNoise` and `finiteSupportDivideRound` helpers. It never imports
GMP/MPFR into the core runtime. Both directions use K+R at the unchanged scale;
R uses conservative ordinary coefficient support H=1 also for the sparse key.
Thus no new system of local arithmetic bounds is introduced.

Write E_enc and E_rest for those canonical semantic errors. In coefficient units:

    M_nu_b <= M_nu + Delta0 E_enc.

This conservatively uses the inverse-embedding coefficient norm <= canonical
sup norm. Preparation requires M_m+M_nu_b < q_src/2. Each centered ciphertext
component has coefficient magnitude at most (q_src-1)/2. Since ||s_b||_1=h,

    |u_j| <= (h+1)(q_src-1)/2,
    |q_src I_j| < (h+2)q_src/2.

For integer I and h>=2, K=h is a valid v9 section 3.1 baseline. K is independent
of decryptions and observed maxima; a tighter deterministic K is not needed.
The canonical raised magnitude is N(h+1)(q_src-1)/(2 Delta0). The restored
magnitude adds E_rest. A separate centered-headroom check uses the exact raised Q.
Restoration error is **not** silently folded into the pre-lift nu_b or I.

`prepareSparseBootstrap`/`executeSparseBootstrap` expose this local stage.
Arithmetic readiness and global security are separate, as for earlier stage
executors. The stage does not issue a full BootstrapCertificate.

With generated sparse keys, the PR-5 planner ignores caller `BootstrapLiftBound`
and derives K from this certificate. Every supported weight currently gives
K>1, so the complete planner returns `UnsupportedEvalRoundDomain` before any
EvalRound/CtS execution. It records the separate security failure as well.
Caller `Analytical` K without generated sparse keys is rejected. The legacy
`TestFixtureAssumption` execution is compiled only in testing builds.

## Security contract

`RlweSecurityEvidence` is a trusted external estimator-result interface, **not**
an estimator implementation. Evidence must identify the family, exact statement,
estimator/version, artifact, assumptions, and finite security lower bound.
Every statement binds the entire public-key graph and publication counts.
Missing, conflicting or stale evidence fails. All required family estimates are
combined by minimum; any Unknown keeps lambda_boot Unknown. A known minimum
below target also fails. The code cannot establish the truth of an external
artifact; accepting its estimator and assumptions is an explicit trust boundary.
No such concrete evidence is provided by this commit.

Ordinary public/evaluation keys and both switching systems are at the actual
key modulus Q_top*P. Encryption records its actual predecessor modulus (SEAL
samples there before divide-round), including multiple encryption levels.
The table also includes deterministic encapsulated/raised/restored images of
one invocation. These are not fresh independent samples; their error law is a
derived convolution, not CBD. Relations explicitly include shared ordinary s,
ephemeral encryption u, gadget messages, automorphisms, and the s <-> s_b cycle.
All generated material is conservatively treated as public. External imported
key material is outside this generated-key certificate and is rejected.

The exact signed fixed-weight search space is binomial(N,h)*2^h. GMP integer
arithmetic plus upward MPFR log2 gives a **search-space ceiling**, never a
security lower bound. Evidence exceeding that ceiling is refused. The ceiling
is also attached to the coupled ordinary families: recovering s_b exposes the
encapsulation messages containing s. It does not replace analysis of RLWE,
low-Hamming attacks, correlations, or the reciprocal-key cycle. SEAL tc128 is
never substituted for a sparse-family estimate.

## Representative measured run

Ordinary context N=4096, primes 68719206401, 68719230977, 68719403009 (special),
Delta0=2^35. SEAL validates the ordinary context at tc128. The candidate sparse
weight is h=64; **its concrete security is unassessed**. The input message
coefficient envelope is 2^28 and original noise envelope 10000. These are
explicit upstream analytical inputs, not measured errors.

| Quantity | Value |
|---|---:|
| Derived K | 64 |
| Encapsulation canonical bound | 0.010498077142803 |
| Restoration canonical bound | 0.020751957723048 |
| nu_b coefficient bound | 360721183.99379 |
| Observed encapsulation canonical error | 6.9846021985275e-6 |
| Observed restoration canonical error | 1.5205370298202e-5 |
| Observed coefficient lift maximum | 9 (certificate remains 64) |
| Before-CtS rho bound | 0.019531489786849 |
| Signed search-space ceiling | 535.29108146207 bits |

Observed errors vary with freshly generated cryptographic keys. They are test
oracle diagnostics and never inputs to certificate preparation. The before-CtS
rho is (M_m+M_nu_b+Delta0 E_rest)/q_src. Final `rho_cert` is Unknown: the K>1
rejection precedes CtS, whose operator/local error has not been added.

Let q=68719206401, Q=q*68719230977, P=68719403009. The actual public inventory
for that run is:

| Family | Modulus | Samples / components | Secret / relation | Security |
|---|---|---:|---|---|
| ordinary.public | QP | 1 / 2 | iid ternary s, zero message | Unknown |
| ordinary.encryption/0 | QP | 2 / 2 | fresh ternary u, shared pk and e0/e1 | Unknown |
| ordinary.relinearization | QP | 2 / 4 | s; gadget s² | Unknown |
| ordinary.galois | QP | 4 / 8 | s; gadget tau_g(s), including conjugation | Unknown |
| sparse.encapsulation | QP | 2 / 4 | fixed-weight s_b; gadget s | Unknown |
| sparse.restoration | QP | 2 / 4 | s; gadget s_b, reciprocal cycle | Unknown |
| derived.encapsulated | q | 1 / 2 | s_b; deterministic image | Unknown |
| derived.raised | Q | 1 / 2 | s_b; deterministic centered lift | Unknown |
| derived.restored | Q | 1 / 2 | s; deterministic image | Unknown |

For encryption, samples counts scalar RLWE-like equations; the two equations
share one u and are the two ciphertext components. For switching keys, each
sample is a public-key pair. This distinction and all correlations are stated
in the estimator input; the counts do not assert independence.

## Tests and remaining blockers

`test_sparse_bootstrap` covers weights 2/4/8/16, exact alphabet and weight,
zero/random messages, both switching oracles, ModRaise identity, measured lifts
versus analytical K, all known local bounds, Unknown and forged-zero rejection,
key/context/scale/inventory invalidation, missing directional keys, explicit
family inventory, and a real ordinary tc128 context with h=64. Canonical errors
are compared with the shared certificates as well as coefficient errors.
`test_sparse_security` uses clearly labelled mock evidence to test all-family
minimum, missing/stale/duplicate evidence, low-weight ceiling and target failure.
Old PR-0..PR-5 test sources and tolerances are unchanged.

Remaining production blockers:

* First: `UnsupportedEvalRoundDomain`, derived K=64 versus the K=1 executor.
* Independently: `SecurityBudgetExceeded`; concrete evidence is missing for the
  required public families, so lambda_boot is Unknown.
* Final CtS `rho_cert` and end-to-end accuracy are unavailable after the domain
  stop. The reported before-CtS rho is not a substitute. No 1e-10 accuracy claim
  is made; previous deterministic accuracy failures remain unchanged.

PR-7, alternative extraction, tails, tuning and fused switching are not included.

Final validation:

* `cmake --build build-pr0 -j 4`: passed (core without GMP/MPFR dependency).
* `ctest --test-dir build-pr0 --output-on-failure`: **15/15**, 15.98 s.
* `cmake --build build-analysis -j 4`: passed.
* `ctest --test-dir build-analysis --output-on-failure`: **26/26**, 149.70 s.
* Standalone `test_sparse_bootstrap` and `test_sparse_security`: passed.
* `git diff --check`: passed.

Changed implementation/API files: `include/m2424/{bootstrap.hpp,
bootstrap_contract.hpp,m2424.hpp,seal_adapter.hpp,sparse_bootstrap.hpp}`,
`src/ckks/seal_adapter.cpp`, `src/core/{sparse_bootstrap.cpp,
sparse_bootstrap_internal.hpp}`, `src/planning/{bootstrap_planner.cpp,
sparse_bootstrap.cpp}`, and `src/CMakeLists.txt`. New tests:
`tests/{test_sparse_bootstrap.cpp,test_sparse_security.cpp,
sparse_bootstrap_oracle.hpp}`; registration in `tests/CMakeLists.txt`.
Documentation: this report and `docs/seal_backend_bound_audit.md`.
