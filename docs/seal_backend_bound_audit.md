# Deterministic backend audit (PR-6, part A)

Audited baseline: bb04c5c05120c164de85698418480b1536a07cf6;
vendored SEAL: a78c3afb01cd2bef1ffe22ef730812738a056795.
Both build-pr0 and build-analysis have SEAL_USE_GAUSSIAN_NOISE undefined.
Result: PASS for the configured CBD backend and the PR-3/4/5 two-component
paths. No shared mathematical bound was changed. This is a finite-support
correctness audit, not a concrete security estimate or an independence claim.

## Source-to-formula mapping

Paths below are relative to extern/seal/native/src/seal/.

* keygenerator.cpp, KeyGenerator::generate_sk and util/rlwe.cpp,
  sample_poly_ternary: independent uniform {-1,0,1} coefficients, represented
  consistently in every RNS limb before NTT. This is not a fixed-weight secret.
  Its coefficient support is 1; its canonical norm is at most N.
* util/rlwe.cpp, sample_poly_cbd: two Hamming weights of 21 bits (8+8+5),
  hence error coefficient support [-21,21]. encrypt_zero_symmetric constructs
  (-a*s-e,a). The same integer error is reduced into all key-modulus limbs.
  The Gaussian configuration is explicitly Unknown in our helper; the CBD
  certificate must not be reused for it.
* keygenerator.cpp, generate_one_kswitch_key: for each data prime q_j, encrypt
  zero at the full key modulus Q_top*P, then add (P mod q_j)*s_old only in
  limb j of component zero. Relinearization uses s_old=s^2 (one switch for
  size 3); Galois keys use s_old=tau_g(s). Multiple removed components would
  require a separate switching charge for each component.
* evaluator.cpp, switch_key_inplace: inverse NTT produces digit d_j with
  coefficients in [0,q_j-1], then lifts that digit to every active q_i and P.
  util/ntt.cpp, inverse_ntt_negacyclic_harvey reduces the lazy output to [0,q).
  Exact modular products accumulate d_j times the evaluation-key pair.
  CRT reconstructs the intended message P*c*s_old plus sum_j d_j*e_j.
  Each convolution has coefficient norm <= N*(q_j-1)*B_e and canonical
  norm <= N^2*(q_j-1)*B_e. After special-prime division and semantic scale S:

      K = N^2 B_e sum_j(q_j-1)/(P S).

  This is exactly finiteSupportKeyNoise. It is noise alone: it does not
  include ModDown rounding. No probabilistic independence is used.
* evaluator.cpp, switch_key_inplace (CKKS special-prime branch): inverse NTT
  the P limb, add floor(P/2), reduce, subtract the resulting centered residue
  in each q_i and multiply by P^{-1} mod q_i. This is nearest integer division
  componentwise. The residual is at most 1/2 per coefficient (P is odd).
* util/rns.cpp, RNSTool::divide_and_round_q_last_ntt_inplace, called by
  evaluator.cpp, mod_switch_scale_to_next: the same centered operation on
  the dropped data prime. With c ciphertext components and secret coefficient
  support H, each residual has canonical norm <= N/2 and |sigma(s)|<=NH:

      R = N/(2 S) sum_{i=0}^{c-1}(N H)^i.

  This is finiteSupportDivideRound. S is the scale AFTER rescale, but the
  unchanged ciphertext scale for special-prime ModDown. No division by P
  belongs on this rounding term. Scale representation is a separate term.
* evaluator.cpp, relinearize_internal and apply_galois_inplace: both call
  switch_key_inplace. Conjugation (g=2N-1) therefore incurs K+R, not zero.
  Automorphisms themselves are exact signed coefficient permutations.

## Hoisting and linear factors

The vendored rotate_vector_many_hoisted lifts before the automorphism. Its
signed digits still have absolute coefficient bound q_j-1.
apply_bsgs_double_hoisted keeps baby products in QP; identity babies lift P*c
without key noise. Nonidentity babies incur K at input scale. Each group
multiplies by its full-key-level encoded diagonals and performs one inner
ModDown. Each nonidentity giant incurs K at product scale; the final sum
performs one outer ModDown. Thus slot_to_coeff_preparation.cpp charges:

    weighted baby K + group_count*R_product + nonidentity_giants*K_product
    + R_product + R_rescale + product/rescale scale-representation errors.

The baby weight includes encoded-operator perturbation, kappa_baby+delta_baby.
There is no missing per-baby ModDown: it is deferred until the group sum.
The second-half phase and final conjugation/add projection are covered by the
same machinery. PR-3 arithmetic adds K+R for every actual relinearization or
conjugation and R for each rescale. Exact prime products/dyadic scales, constant
encoding errors and centered headroom checks remain separate obligations.
The identities are modulo the active modulus: semantic interpretation also
requires the existing no-wrap checks; finite support alone does not prove it.

## Limits

This audit establishes these formulas only for the operations and sampler
above. It does not certify a sparse security level, arbitrary externally
loaded key distributions, or a different SEAL backend. Ordinary secret H=1
is a deterministic support assumption, not a claim that all coefficients are
nonzero. Bounds are deliberately conservative; observed errors are not used.

## PR-6 identity switching extension

`src/ckks/seal_adapter.cpp::identitySwitchKey` follows the audited
`generate_one_kswitch_key` limb construction under an explicitly different
destination secret. SEAL's public `apply_galois_inplace` accepts g=1, applies
an identity permutation (it does not return early), then uses precisely the
same `switch_key_inplace`. This differs from `rotate_internal(steps=0)`, which
would return early and is deliberately not used. Both directional switches
are charged separately. The source/destination test-only oracles exercise the
actual operation, including its special-prime ModDown.
