# PR-5: deterministic EvalRound+ pipeline

Base revision: `2e2938a13794fee8c128e384fca848a1645405ae`.
Mathematical source: `ckks_bootstrapping_model_v9.pdf`, EvalRound+ decomposition
and the global error contract. Scope stops at the first complete baseline.
No sparse-secret production implementation, probabilistic tails, LCR, AKS,
Grafting, Thrifty evaluation or later optimization is added.

## Execution and certificate boundaries

The new public API is `Bootstrapper::prepare`, `preflight`, and `apply`, with
an immutable `BootstrapPlan`. The executed chain is:

```
input
  -> centered ModRaise
  -> EvalRoundPlusCoeffToSlot::prepareCertified / apply (HP and LP)
  -> EvalRoundExecutionPlan / executeEvalRound (each LP half)
  -> certified scalar arithmetic: HP - (q_src / Delta0) * I_hat
  -> PreparedSlotToCoeffPlan / apply
  -> output
```

Preparation remains in the existing optional `m2424_evalmod_analysis` target.
Core execution has no new GMP/MPFR dependency. Parameter tokens are public
zero encryptions used solely to call the existing prepared-plan APIs at
future levels/scales; they are not output operands or sources of bounds.
Production planning/execution does not call decryption or inspect the secret.

The arithmetic builder and interpreter from PR-3 are shared with combination.
Their original one-input EvalRound behavior is retained. Combination adds
two certified inputs and an explicit downward modulus transition. It has no
metadata scale-rewrite operation. CtS/StC keep the PR-1-cert/PR-4 propagation
and certified root-diagonal preparation.

`executionReadiness()` means all executable stages have valid local proofs
and compatible actual schedules. `trace().result` is the global result. A
locally valid CtS bound whose requested `errorBudget` fails is recorded as a
global failure and is never promoted by a later successful gate. Rejected
plans produce no output; `apply` preserves the first planning failure.

## Required lift evidence

`BootstrapLiftBound` has an optional K, an evidence kind and provenance.
Default evidence is Unknown. There is no implicit K=1 in ordinary ModRaise.
Unknown evidence returns `LiftBoundUnavailable`; any supplied K other than
1 currently returns `UnsupportedEvalRoundDomain`, before EvalRound is
compiled. A caller claiming analytical evidence supplies an explicit
conditional input contract; this implementation does not derive production
K=1 from a normal ciphertext.

`TestFixtureAssumption` is separately marked and always prevents global
production certification. The ordinary sparse-secret production gate also
remains unavailable, and actual SEAL security status is recorded separately.
These checks cannot be bypassed by a successful mathematical EvalRound plan.

## Backward budget and combination

The LP radius is taken directly from the certified CtS domain result:

```
rho_cert = (M_m + M_nu) / exact q_src + E_LP.
```

A missing upstream bound remains Unknown. rho >= 1/2 rejects before EvalRound.
The source modulus is an integer product of actual primes; Delta0 is the
exact dyadic represented by the input scale bits. Gamma multiplication uses
that exact rational throughout constant preparation.

The planner first fixes the output-side modulus/scale schedule and obtains a
real StC certificate. Its actual `gamma` supplies `Gamma_StC`. To reserve
arithmetic that precedes StC, StC is initially prepared with a pair-error
envelope equal to the global target. The actual pair error is checked against
this envelope. The StC intrinsic term is evaluated using the **same** factor
bounds with zero incoming error, retaining local bounds prepared for the full
envelope. Thus the certified affine estimate is

```
E_StC <= Gamma_StC * E_pair + B_StC_intrinsic.
```

Zero incoming error here isolates an affine term; it does not replace an
unknown input error. With all other terms reserved, the remaining budget is

```
remaining = target
          - N*M_nu/Delta0
          - B_StC_intrinsic
          - Gamma_StC*(E_HP + B_comb_reserved)
requiredIntegerError = remaining / (Gamma_StC * exact gamma).
```

The last division is rounded downward. Nonpositive remaining budget prevents
EvalRound compilation/execution. Neither N nor the ideal StC norm is used
instead of the actual prepared gain.

The combination reservation uses the shared finite-support divide-round
helper, configured minimum arithmetic/constant scales, and a conservative
normal-binary64 relative scale ratio bound 2^-50. For each allowed integer
rescale it reserves component rounding and representation error. Constant
encoding reserves `1/(2*minimumConstantScale)` times the input envelope; each
final alignment rescale reserves its actual output-scale rounding. Two
binary64 operations and prime-to-double conversion are covered by 2^-50,
which exceeds their combined unit-roundoff bound for the checked normal
positive scales. No probabilistic argument is used.

After the backward pass, the PR-2 planner and PR-3 compiler select the actual
binary extraction/cleaning schedule from `requiredIntegerError`. The actual
combination DAG then:

1. rescales the integer ciphertext if needed, within the reserved count and
   minimum arithmetic scale;
2. modulus-switches HP and integer to the reserved common prime base;
3. multiplies HP by an encoded exact one and integer by encoded exact gamma,
   using real plaintext scales chosen so their subsequent rescales have
   identical requested scale bits;
4. rescales both operands and subtracts them.

If no such multiply/rescale scale exists, planning rejects. All constant
encoding, finite-support rescale rounding and exact scale representation
errors are recomputed by the shared arithmetic builder. The actual B_comb
must fit the reservation. The resulting contract is

```
E_pair <= E_HP + exact gamma * E_I + B_comb.
```

The final StC preparation uses this actual E_pair. Its gain must equal the
previously prepared gain. The final aggregation is

```
E_boot = actual StC.outputError + N*M_nu/Delta0.
```

StC.outputError already includes amplification of E_pair. The planner does
not add `Gamma_StC * E_pair` a second time. The target is the canonical
embedding of the source message polynomial m/Delta0.

## Unified trace

`BootstrapTrace` retains exact source/raised primes and scale in its plan
metadata, plus ordered gate results, named bounds, selected extractor/cleaning
counts, exact gamma, planned output scale, key requirements, failure-event
status and security status. Unreached quantities are explicit Unknown bounds
with the blocking gate's provenance.

During execution it appends actual active primes, chain index, scale bits,
ideal magnitude, semantic error and local arithmetic error for every CtS,
EvalRound, combination and StC callback. Existing stage executors compare
actual metadata with their compiled schedules. The headroom gate records the
maximum exact centered utilization across compiled output states, alongside
the independently checked internal BSGS QP headroom proofs.

## Explicit test-only integration fixture

The integration test uses N=16, 15 actual 50-bit primes including the special
prime, input scale 2^49, and CtS/StC depths two. This tiny fixture explicitly
uses SEAL security level `none`; normal `SealAdapter::create` and all existing
profiles/security settings are unchanged. Its factory is compiled only with
BUILD_TESTING and declared in `tests/bootstrap_fixture.hpp`.

The fixture ciphertext is synthetic `c=(m,1)`, not production encryption.
Its source noise is the ordinary ternary secret polynomial, with coefficient
support one. The fixture constructor neither reads the secret nor decrypts.
Messages use exact coefficients `(j mod 5 - 2)*2^30`. They remain far from
source wrap, so this fixture has zero actual lifts and satisfies its explicit
external K=1 assumption. It tests both coefficient halves, phase handling,
nonzero LP arguments, cleaning, real gamma arithmetic and reconstruction.
The existing EvalRound tests continue to cover all three integer centers.

The test uses a separate integration budget **1e-6** so the deterministic
pipeline can execute; it does not claim feasibility at 1e-10. The same test
also requests 1e-10 and a tighter budget and checks rejection before EvalRound.
The global successful-execution result remains `TestOnlyAssumption`, never
`Certified`, irrespective of observed accuracy.

Recorded fixture schedule:

| Stage | Chain transition | Scale at stage output |
|---|---|---:|
| Source / ModRaise | 0 -> 13 | 562949953421312 |
| CtS, each HP/LP half | 13 -> 11 | 562949953422159 |
| EvalRound, each LP half | 11 -> 8 | 1.0141204802170829e31 |
| Integer scale reduction | 8 -> 7 | 9007199255068904 |
| HP / integer modulus alignment | 11 / 7 -> 5 | unchanged per operand |
| Real scalar multiply and rescale | 5 -> 4 | 562949953421312 |
| StC | 4 -> 2 | 562949953425279 |

Eleven chain levels are consumed after ModRaise, including two skipped by
integer modulus switching. CtS consumes two rescales per half; EvalRound
consumes three levels on its longest path; combination consumes one integer
conditioning rescale and one final rescale per operand; StC consumes two
levels. HP alignment is real arithmetic, not a metadata change.

Representative deterministic bounds:

| Quantity | Value |
|---|---:|
| rho_cert | 1.9077466371209822e-6 |
| K | 1, external test-only assumption |
| exact gamma, shown rounded | 1.9999999999911342 |
| actual Gamma_StC | 16.000000000000213 |
| required E_I | 3.0633557100258561e-8 |
| achieved certified E_I | 1.5580967743616881e-10 |
| E_HP | 4.3142734234725787e-10 |
| E_LP | 3.9800341184839951e-10 |
| actual B_comb | 5.1869619982602422e-13 |
| reserved B_comb | 7.4400752445976486e-10 |
| E_pair | 7.4356539341804014e-10 |
| E_StC, including input amplification | 1.2816232800124964e-8 |
| E_boot | 1.2816261221834394e-8 |
| source-noise output contribution | 2.8421709430404007e-14 |

The selected extractor is `BinaryQuadraticK1/SEAL-baseline`, radix two, with
one cleaning iteration per digit. The conditional primary-event union has
failure probability zero: no probabilistic events were introduced. This does
not remove the test-only lift assumption or the failed production gates.

One oracle run observed end-to-end error `5.5137861795e-11`. This varies with
random evaluation keys and is not a certificate. Each executed node is
compared with an independent exact-rational/MPFR polynomial oracle or the
canonical embedding reference, and checked against its stage certificate.

## Why 1e-10 is not certified

For an ordinary production input, the first gate is `LiftBoundUnavailable`
unless external analytical lift evidence is supplied. Larger K is explicitly
unsupported by the current ciphertext EvalRound executor. No observed lift
is used to change that decision.

Even this small fixture has E_HP and E_LP above 1e-10. At target 1e-10 the
first global failure is the CtS HP branch error budget; the backward budget
also becomes nonpositive and no EvalRound plan is compiled. The reserved
fixed output-error contribution is about `1.9726172796058213e-8`.

In the completed fixture arithmetic, the largest output contributions are
approximately `Gamma_StC*E_HP = 6.90e-9` and
`Gamma_StC*gamma*E_I = 4.99e-9`. StC intrinsic arithmetic contributes about
`9.19e-10`; actual combination local arithmetic contributes about `8.30e-12`
after StC amplification. Source noise is smaller still. Bounds are retained
without tail improvements or tolerance changes.

The PR-1-cert production-profile bounds remain as previously reported:
E_HP `0.0110686011341`, E_LP `0.0110729811405`; PR-4 standalone StC remains
`0.0340233987753`. They were not replaced by the tiny fixture's numbers.
Those are different profiles and cannot be spliced into a fabricated
end-to-end certificate.

Remaining production blockers are lift evidence/supported K, deterministic
accuracy, an adequate real level/scale schedule, and the unavailable
sparse-secret production/security certificate. The tiny fixture additionally
has explicit security level zero and test-only lift evidence. Unknown gates
remain Unknown, and no later PR is implemented here.
