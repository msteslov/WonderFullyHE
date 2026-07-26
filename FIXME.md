# FIXME

## Bootstrap research blockers

These items are deliberately retained as research work, not exposed as a
successful refresh implementation.

- **EvalMod P3 scale alignment.** In the STC-first FFT-like path, `u^3` can
  reach a scale larger than the desired output scale. The coefficient-weighting
  step then requires a non-positive plaintext scale and aborts with
  `cannot align EvalMod term scale`. A correct fix must define a level-aware
  scale-squash policy for the cubic branch and align the linear branch to it.

- **SmallSlots4Butterfly level budget.** The backend consumes more transform
  levels than the current STC-first default target chain index preserves. It
  reaches a final modulus where the next plaintext multiplication is out of
  bounds. The required input chain index must be derived from the selected DFT
  plan instead of using a shared fixed value.

- **Encrypted ScaleDown/ModUp semantic invariant.** Structural ModUp does not
  yet prove the pre-EvalMod lattice form required for a value-preserving
  bootstrap. Keep using the reference and encrypted-vs-reference diagnostics
  before promoting a refresh configuration.

## Removed prototype

`ModUpThenEncode` was removed from the public research API because it had no
validated end-to-end scenario and duplicated an unproven ScaleDown/ModUp path.
It must not be restored without a checked encrypted test that verifies value
preservation, level restoration, and the lattice invariant.
