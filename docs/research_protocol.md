# Research protocol

## Research question

How do CKKS parameter selection, multiplicative depth, parallel execution, and
CoeffToSlot factorization affect the joint trade-off between latency, memory cost,
throughput, and numerical error in WonderFullyHE?

The practical goal is to identify the least expensive configuration that satisfies
an explicit error threshold. This is a systems and experimental contribution; the
current repository does not claim a new cryptographic construction or a fully
certified production bootstrap.

## Hypotheses

1. Numerical error grows with multiplicative depth, while adding unused modulus
   levels does not materially improve a fixed-depth computation.
2. Matching scale and working-prime size improves depth-two accuracy, subject to
   the modulus budget and security constraints.
3. Reusing a prepared CoeffToSlot plan reduces warm-application latency relative
   to the first application.
4. Parallel throughput improves with worker count but remains sublinear because
   every worker owns a context and evaluation keys.
5. Bootstrap key setup and serialized key volume increase with polynomial degree
   and the supported rotation set.

## Experimental design

The canonical experiment is `python/configs/research.json`. It covers:

- five built-in CKKS profiles;
- chain length, scale, working-prime size, and multiplicative depth;
- one, two, four, and eight worker threads;
- all registered bootstrap candidates;
- all five ranked CoeffToSlot factorizations.
- a complete test-only ModRaise → CoeffToSlot → EvalRound → SlotToCoeff fixture.

Each external configuration is repeated at least three times. Some C++ benchmarks
also perform internal trials. Raw stdout and stderr are retained. The manifest
records the Git commit, dirty-worktree state, exact commands, platform, Python
version, return codes, and wall-clock duration.

## Measurements

- operation and end-to-end latency in milliseconds;
- throughput in tasks per second;
- maximum and mean absolute error;
- serialized public, relinearization, and Galois key sizes;
- CoeffToSlot operation counts and prepared-plaintext bytes;
- mean, median, minimum, maximum, standard deviation, and 95% confidence interval.

The 95% interval is the normal approximation `1.96 * s / sqrt(n)`. It describes
run-to-run variability; it is not a proof of cryptographic correctness.

## Acceptance and interpretation

A configuration is admissible only if the underlying benchmark reports success
and its observed error is below the relevant contract threshold. Performance is
compared only among admissible configurations. Results from a dirty worktree must
be labelled as exploratory and rerun from a clean commit before publication.

The end-to-end fixture verifies integration and numerical behavior of every
arithmetic stage at insecure `N=16`. Its `TestOnlyAssumption(K=1)` is recorded in
every output row. The production bootstrap remains a limitation until its lift and
EvalRound contracts are certified without test-only assumptions; fixture success
must not be presented as production security evidence.
