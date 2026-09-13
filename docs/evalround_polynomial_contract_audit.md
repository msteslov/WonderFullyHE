# PR-3 polynomial follow-up audit (before implementation)

Base dd2ee4721045053560214e6789d7ce1d89d78277.
Reviewed evalround.hpp, evalround_execution.hpp, the reference planner/compiler,
evalround_reference.cpp, evalmod_synthesis.cpp's compileEvalModPolynomial and
MPFR evaluator, approximation_lab.cpp's interval certificate, and v9 section 4.

No concrete whole-domain-certified general-K extraction polynomial/circuit is
present. v9 supplies the disjoint interval domain, digit semantics, candidate
routes and cleaning/reconstruction inequalities, not instantiated general-K
coefficients and their error certificate. BinaryQuadraticK1 is the existing
analytic polynomial exception. PiecewiseReference and exact ternary phase are
reference targets, not ciphertext implementations. K=64 remains blocked on an
actual polynomial list and an error certificate for every offset-binary digit.

EvalModPolynomial already stores basis and exact decimal coefficients, parsed
as rationals by parseExactDecimal. Its storage/parser can be shared without
making GMP a core dependency. The old experimental evaluator's arithmetic
certificate is not transferable: the new compiler must use the current Builder.

The MPFR interval multiplication/Horner enclosure mathematics applies directly
to a fixed polynomial minus a constant digit on each closed interval. This
needs no new extraction assumption: the integer offset digit is already defined
by v9. The EvalMod target subtraction (p(x)-(x-I)) and derivative correction
must not be reused as a digit bound. A separate target-specific verifier may
reuse the same outward interval primitives over complete cells, subtracting
bit_j(I+K), without synthesizing a polynomial or accepting grid diagnostics.
