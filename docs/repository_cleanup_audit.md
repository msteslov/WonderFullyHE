# Repository cleanup audit

Baseline commit: `314cec9e1a02bcb75cca0e9a5f4d516934f9ca71` (`master`).

## Baseline

- Release configuration: `BUILD_TESTING=ON`, `M2424_BUILD_APPS=ON`,
  `M2424_ENABLE_EVALMOD_ANALYSIS=ON`.
- Build: passed.
- CTest inventory: 29 tests. Passed: 26. Timed out: 3.
  `test_evalmod_feasibility` printed `PASS` before the generic 60 s timeout;
  `test_sparse_coeff_to_slot` printed its full certificate trace and `PASS` before
  timeout; `test_evalround_k64_synthesis` timed out without an assertion report.
- Smoke executables passed with exit code 0: `demo_basic`,
  `demo_profile_report`, `demo_security_report`, `bench_ckks`,
  `bench_chain_accuracy`, `bench_bootstrap_candidates`.
- Libraries: `m2424`, optional `m2424_evalmod_analysis`, `seal`.
- Stable app targets: `bench_chain_accuracy`, `bench_bootstrap_candidates`,
  `bench_ckks`, `bench_coeff_to_slot`, `bench_parallel_throughput`, `demo_abft`,
  `demo_basic`, `demo_noise_growth`, `demo_precision_profiles`,
  `demo_profile_report`, `demo_secure_stats`, `demo_security_report`.
- Optional report targets: `report_evalmod_synthesis`,
  `report_evalmod_feasibility`, `report_evalround_k64_synthesis`.
- Public headers: 21 headers directly under `include/m2424/`; experimental API:
  12 headers under `include/m2424/experimental/evalmod_analysis/`.
- Source layout: `src/core` 17 files, `src/ckks` 3 files, `src/math` 9 files,
  `src/planning` 15 files before cleanup.

The baseline target inventory also contained all 29 test executables and the
standard CTest/Doxygen targets. No target was removed.

## Component audit

| Component | Current responsibility | Problem | Action | Behavior changed? | Tests covering action | Follow-up optimization |
|---|---|---|---|---|---|---|
| Public API | Stable contracts/plans/runtime types | Boundary was described only indirectly; umbrella include order was misleading | Documented stable/experimental/internal boundary; normalized umbrella include placement | no | full build, public-header consumers | consider install/export compatibility test |
| Experimental analysis | Exact arithmetic, reference/oracle, synthesis/search, DAG compilation | Production planning and research analysis coexist in one optional target | Kept opt-in GMP/MPFR target; documented file-level roles | no | EvalMod/EvalRound analysis tests | split production plan preparation from research search only in a separate API PR |
| EvalRound runtime | Preflight and execution of certified DAG with liveness | Same filename stem as compiler obscured responsibility | Runtime name retained; internal runtime declaration renamed explicitly | no | `test_evalround_backend`, `test_evalround_certificates`, `test_evalround_polynomial` | runtime allocation/liveness instrumentation |
| EvalRound compiler | Exact compilation and certification of DAG | `planning/evalround_execution.cpp` looked like a second executor | Renamed to `evalround_execution_compiler.cpp`; planning builder renamed | no | planner/backend/certificate/K64 tests | cache immutable compiled constants after profiling |
| Bootstrap planner | Builds the only certified end-to-end plan | Its dependency on optional exact analysis was not explicit | Documented planning path and optional link boundary | no | `test_bootstrap_end_to_end`, sparse tests | unified whole-bootstrap cost model |
| Bootstrap runtime | Executes ModRaise→CtS→EvalRound→combine→StC | Reference/test paths were near production backend code | Documented single path; moved native test fixture/oracle implementation out of production library | no | bootstrap end-to-end, ModRaise, sparse bootstrap | peak memory instrumentation |
| CoeffToSlot | Raised input, factorized BSGS, dual HP/LP certification/execution | Architecture appears asymmetric with StC without rationale | Documented common lifecycle and intentional dual-branch differences | no | CtS contract/cert/backend, double-hoisted BSGS, sparse CtS | prepared diagonal storage/copy audit |
| SlotToCoeff | Factorized inverse transform and two-half combination | Preparation and certificate location differed from CtS | Documented lifecycle and certificate ownership | no | StC contract/backend tests | scratch reuse after ownership design |
| SealAdapter core | Context, encode/decode, arithmetic and RNS ModRaise | One 1234-line implementation mixed unrelated ownership areas | Added private PIMPL header; extracted serialization and sparse-key units | no | primitives, adapter failures, ModRaise, bootstrap/EvalRound serialization paths | deep copy `Cipher`; deep copy `Plain`; copy-before-inplace in `relinearize`, `rescaleToNext`, `modSwitchTo` |
| SealAdapter keys | Ordinary and sparse evaluation keys | Key state shares one adapter PIMPL | Preserved ownership and semantics; isolated sparse implementation | no | sparse security/bootstrap and key failure tests | key ownership and shared evaluation keys |
| Test fixtures/oracles | Unsafe-context fixture and secret-key diagnostics | Implemented inside production `seal_adapter.cpp` under test macros | Moved definitions to `m2424_test_support` | no | fixture-dependent ModRaise/CtS/bootstrap tests | remove public compile-time diagnostic surface in a compatibility PR |
| Apps | Demos, benchmarks and reports | Flat list hid categories | Added explicit CMake groups; paths and target names unchanged | no | Release build and smoke runs | worker-local `SealAdapter` in parallel benchmark |
| Tests | Contract, backend, integration and certification coverage | Flat target list and insufficient host-load timeouts | Grouped stable tests; increased timeouts only for baseline timeout cases | no | complete CTest | shard heavyweight synthesis/certificate tests |
| CMake | Production plus optional analysis | Source ownership was implicit in one list | Added named core/backend/math/planning source groups; optional analysis remains isolated from production GMP/MPFR | no | analysis ON build and production-only build | real installed-target dependency test |

No assertion, numerical tolerance, certificate gate, security condition, CKKS
profile, factorization preset, polynomial coefficient or operation order changed.

## Post-cleanup verification

- Release build with tests, apps and experimental analysis: passed.
- Production-only Release build with tests and experimental analysis disabled:
  passed; this verifies that `m2424::m2424` does not require GMP/MPFR or
  `m2424_test_support`.
- Full CTest: 29/29 passed in 372.90 s. The three baseline timeout cases passed
  unchanged with host-load-appropriate limits.
- Focused adapter/ModRaise/sparse/bootstrap subset: 6/6 passed.
- Post-change smoke executables: the same six baseline commands all returned 0.
- No CMake target or public header was removed. Public API signatures are
  unchanged.

## Removed or retained code

No runtime path or public API was deleted. Test fixture/oracle definitions were
moved, not removed. The two EvalRound execution files and the two certified
arithmetic internal headers were retained after dependency inspection proved
that they have distinct compiler/runtime responsibilities. Historical archive
documents retain old filenames as historical records.

No uncertain dead code was removed. The diagnostic `restoreSparse` path remains
compile-time test-only because current sparse tests depend on it. The stable
headers whose `prepare` definitions live in the optional analysis library remain
unchanged to avoid a public API break.

## Technical debt not changed in this cleanup

Priority order for follow-up PRs:

1. Define ownership and benchmark deep copy `Cipher` and deep copy `Plain`;
   then address copy-before-inplace in `relinearize`, `rescaleToNext` and
   `modSwitchTo` without changing the public API accidentally.
2. Separate key ownership from worker state and evaluate shared evaluation keys;
   replace worker-local `SealAdapter` in the parallel benchmark only after that
   ownership contract exists.
3. Define memory-pool/scratch ownership and add peak-memory instrumentation.
4. Add a unified whole-bootstrap cost model using the existing certificates and
   metrics; do not make it a certification authority.
5. Specify the CPU/GPU backend boundary below runtime execution.
6. Design a Boolean DAG layer above the existing arithmetic DAG without changing
   EvalRound certification semantics.

Also retained: optional-preparation methods declared in stable headers but
implemented in `m2424_evalmod_analysis`; compile-time diagnostic members exposed
when `BUILD_TESTING=ON`; flat physical `apps/` and `tests/` layouts to avoid a
large path-only diff.
