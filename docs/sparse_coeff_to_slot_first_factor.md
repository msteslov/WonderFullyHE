# PR-6 follow-up: restoration belongs to the first CtS factor

Base revision: 22a640b717e52542dac413692c2055a9c3c1488c.
The mathematical contract is v9 section 3.1–3.2; this follow-up adopts the
implementation specification's first-linear-factor ownership of restoration.
It does not change K=h_b, the identity-switch algorithm, finite-support
formulas, the nine-family security model, or any previous test tolerance.

## Exact execution sequence

1. `encapsulateSparse`: switch s -> s_b at the exact source modulus/scale.
2. `modRaiseSparse`: centered lift to the actual raised modulus. The result is
   still a distinct `SparseRaisedCipher`, decrypting under s_b. Its ideal value
   is y=sigma(u/Delta0), u=m+nu_b+q_src I.
3. `EvalRoundPlusCoeffToSlot::apply(SparseRaisedCipher, certificate)` preflights
   the sparse generation, context, source/raised primes, scale, public inventory,
   restoration key and all linear-factor rotation/conjugation keys.
4. Inside that first-factor executor, one shared `restoreSparseForFirstFactor`
   calls the existing identity-automorphism key switch s_b -> s. This includes
   the audited special-prime ModDown. It preserves the modulus and scale.
5. HP and LP each start their existing BSGS factorization from this restored
   ciphertext under ordinary s, followed by each factor's existing rescale.
   The LP first factor, including the phase-corrected second-half factor,
   still contains exact alpha=Delta0/q_src in its encoded diagonals.
6. Remaining factors and the final conjugation/add projection operate under s.

There is exactly **one** restoration for the four HP/LP half outputs. It is
owned by the CtS first-factor executor, not a separate Bootstrapper stage.
Preparation performs no restoration. The source `SparseRaisedCipher` remains
unchanged; no public conversion to ordinary `RaisedCipher` is used by this path.

The legacy `restoreSparse` and `executeSparseBootstrap` APIs are compiled only
with `M2424_ENABLE_SPARSE_DIAGNOSTICS`, enabled by BUILD_TESTING. Their existing
unit tests remain unchanged. Production builds without testing do not expose
those standalone APIs; the private first-factor implementation remains available
only to its CtS owner. Both production Bootstrapper source files contain no
standalone restoration call.

## Error ledger and recurrence

The existing shared helpers produce the same E_rest as before:

    E_rest = finiteSupportKeyNoise(actual raised primes, P, Delta0)
           + finiteSupportDivideRound(N, H=1, components=2, Delta0).

The sparse certificate now also exposes those two summands separately with
provenance. CtS preparation checks that the total and both summands are known,
positive and at least the configured backend's shared finite-support bounds.
Unknown, missing, zero or understated restoration envelopes are rejected.
There is no new key-switch bound implementation.

For every HP/LP half, the ideal incoming magnitude is the **unrestored**
`raisedMagnitude`, while the incoming semantic error is E_rest:

    M_0 = bound on |sigma(u/Delta0)|,
    E_0 = E_rest.

The existing `prepareRootLinearTransform` machinery then applies, at every factor:

    M_next <= kappa M,
    E_next <= kappa E + delta (M+E) + B_local.

The final projection continues the same recurrence. E_rest is never added to
the final branch error separately. It is not part of B_local. Encoding/scale
terms that depend on the full M+E envelope retain that dependency; this is
operator perturbation of the actual input, not a second additive restoration.

`inputSemanticError`, `firstFactorRestoration`, and each factor's
`incomingSemanticError` expose this distinction. The first factor's
`arithmeticTerms` identifies incoming restoration key noise and ModDown as
**not B_local**. Its existing terms retain first-factor operator perturbation,
BSGS runtime error, rescale rounding and scale-representation provenance.
The observer emits one shared restoration event followed by the unchanged
40 BSGS/rescale/projection events for the selected depth-four profile.

Encapsulation error still belongs to nu_b. Restoration does not change nu_b or I:

    HP = u^(a)/Delta0 + eta_HP,
    LP = I^(a) + (m^(a)+nu_b^(a))/q_src + eta_LP.

The propagated restoration contribution belongs to eta_HP/eta_LP. Thus the
existing domain function computes exactly

    rho_cert = (M_m+M_nu_b)/exact q_src + E_LP.

It does not add E_rest again. `rho.beforeCtS` remains a diagnostic and is not
used as the EvalRound domain certificate or as another incoming error.
No final bootstrap-output restoration term is introduced.

## Full planner and security

The sparse planner now retains a prepared sparse CtS certificate, available
through `BootstrapPlan::coeffToSlot()`. It records E_HP, E_LP and final rho_cert
when the CtS prerequisites are available. It then rejects K>1 with
`UnsupportedEvalRoundDomain`, without constructing/executing an EvalRound DAG.
If CtS keys/levels are unavailable, their diagnostic rejection is retained;
it cannot authorize substituting a smaller K. This preserves earlier PR-6
small-profile rejection behavior.

No new public key material is introduced. The nine family statements, exact
moduli, counts, relations and security evidence rules are unchanged. The test
compares their full statements before and after preparation. Each required
family remains Unknown; lambda_boot remains Unknown/`SecurityBudgetExceeded`.
The sparse weight in the test is a candidate, not a claimed concrete security
estimate. General-K EvalRound, PR-7 and optimizations remain out of scope.

## Selected sparse test profile and results

The new test reuses the existing certified CtS parameter profile:
N=16384, seven 60-bit primes, Delta0=2^59.5, depth four, factorization [6,5,5,3].
The sparse candidate has h_b=64 and the certificate retains **K=64**.
The message is the exactly encoded zero polynomial; the upstream noise envelope
is the same finite-support asymmetric-encryption bound used by existing tests.
No observed value is an input to a certificate.

| Quantity | Certified | Observed in the standalone validation run |
|---|---:|---:|
| Restoration canonical semantic error | 4.1652950178397e-8 | 4.7621153291824e-9 |
| Final HP error | 5.2020794443499e-5 | 1.8511862250863e-11 |
| Final LP error | 5.2015379743736e-5 | 1.8827384658613e-11 |
| Final rho_cert | 5.202038621366e-5 | not used to derive the bound |

Fresh key generation changes observed errors; the deterministic certificates
are unchanged. Every observed stage error is below its certificate. The maximum
observed/certified ratio was 0.114328405 including restoration, and 8.18349e-6
among the 40 subsequent CtS stages. Both branches consume four levels and eight
half-branch rescales, exactly as before; restoration consumes no level/rescale.

The full bootstrap's first blocker remains `UnsupportedEvalRoundDomain` (K=64).
Independently, concrete all-family security evidence is missing, and these
honest deterministic branch errors exceed the requested 1e-10 budget. Nothing
in this follow-up claims a Certified complete bootstrap.

## Validation coverage

`test_sparse_coeff_to_slot` checks sparse type/semantic preservation, exactly
one actual restoration operation, immutable source ciphertext, missing keys
before arithmetic, generation/context invalidation, Unknown/forged-zero bound
rejection, every runtime stage against a test-only oracle, both independent
coefficient-half targets, exact recurrence replay from E_rest, tight outward
rho with no duplicate restoration term, folded normalization level counts,
K=h_b, unchanged security statements and the full K>1 execution rejection.
The earlier PR-0..PR-6 test source files/tolerances remain unchanged; the existing
oracle helper only gains a test-only operation counter.

Final validation:

* Baseline build passed; full baseline ctest: **15/15**, 15.69 s.
* Analysis build passed; full analysis ctest: **27/27**, 219.00 s.
* The new sparse CtS test also passed as a standalone invocation.
* A clean Release core build with `BUILD_TESTING=OFF` and
  `M2424_ENABLE_EVALMOD_ANALYSIS=OFF` passed. A compile-time API check confirmed
  neither legacy standalone restoration API is exposed in that configuration.
* `git diff --check` passed.

Changed files: public headers `bootstrap.hpp`,
`evalround_plus_coeff_to_slot.hpp`, `seal_adapter.hpp`, `slot_to_coeff.hpp`,
`sparse_bootstrap.hpp`; `src/CMakeLists.txt`; `src/ckks/{bootstrap.cpp,
seal_adapter.cpp}`; `src/core/{coeff_to_slot_cert_internal.hpp,
sparse_bootstrap.cpp}`; `src/math/coeff_to_slot_cert.cpp`;
`src/planning/{bootstrap_planner.cpp,coeff_to_slot_preparation.cpp,
sparse_bootstrap.cpp}`; `tests/{CMakeLists.txt,sparse_bootstrap_oracle.hpp,
test_sparse_coeff_to_slot.cpp}`; this report and the historical baseline report
link update.
