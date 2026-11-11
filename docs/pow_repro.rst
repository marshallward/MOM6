Reproducible real powers
========================

``MOM_intrinsic_functions`` exports the elemental default-real function
``pow_repro(a, b)``, implemented by the ``MOM_pow`` submodule. Both arguments
are real. Like the existing reproducible exponential and logarithm, its
intended arithmetic is IEEE binary64 (MOM's default-real-8 builds).

Algorithm and accuracy
----------------------

For a positive finite base, the general algorithm uses the mathematical
identity ``a**b = exp(b*log(a))``, but carries high/low pairs through the
intermediate calculations. Simply composing the public ``log_repro`` and
``exp_repro`` functions loses the logarithm and product rounding residuals;
exponentiation amplifies those errors to hundreds of output ULPs.

The dedicated power kernel instead:

#. Reduces ``a = 2**k * z`` and uses the existing reciprocal logarithm table
   to form ``1 + r = z * inv_c``. The multiplication residual is retained.
#. Reconstructs ``log(a)`` as a high/low pair from ``k*log(2)``, the table
   logarithm and its tail, and a compensated degree-10 ``log1p`` polynomial.
   Inputs close to one use a direct polynomial to avoid cancellation of
   table constants when a large exponent will amplify the residual error.
#. Multiplies this pair by ``b`` with compensated, bit-split multiplication.
   No hardware FMA or runtime extended-precision library is required.
#. Passes both product parts to a private exponential interface. The low
   part is added after exponential range reduction, before polynomial
   evaluation and final scaling. The public exponential's arithmetic is
   unchanged when the optional residual is absent.

``MOM_pow_log_data`` supplies the extra logarithm table tails. Its Python
generator uses 100-digit Decimal arithmetic and the exact binary64 values
of the existing logarithm table. Regenerate it if that table changes::

  cd src/framework/intrinsics
  python3 generators/generate_MOM_pow_log_data.py > MOM_pow_log_data.F90

For negative finite bases, an integral exponent is required; the magnitude
is evaluated with ``abs(a)`` and the sign is restored for odd exponents.
Integrality and parity checks do not convert large real exponents to integers.

Results are near roundoff accuracy, but are not guaranteed correctly rounded.
Quad-reference tests enforce a 2-ULP bound on general powers, including
large exponents near one and logarithm table boundaries, and a 3-ULP bound
on a broader grid that also exercises the nested-root shortcuts. These
tested bounds are regression thresholds, not a proof for every input.
Additional tests cover subnormal inputs and results and range limits.

Fast paths avoid unnecessary transcendental evaluations:

* Exponents zero and one preserve the corresponding exact identities.
* Exponents two and minus one use direct multiplication and division.
  Other integral exponents use the compensated general path: repeated
  multiplication can accumulate substantially more rounding error.
* Exponents with magnitude 1/2, 1/4, or 1/8 use one, two, or three square
  roots, followed by a reciprocal for negative exponents. Each square
  root is correctly rounded, but the composed operation need not be.

The bounded set of root shortcuts avoids an arbitrary-length root loop
on the general path. The timing driver compares runtime intrinsic powers,
``pow_repro``, the direct exp/log composition, and explicit nested roots.
Compiler optimization of constant powers can differ from runtime powers;
the tests do not assume a particular Intel or GNU library implementation.

Special values
--------------

* ``a**0`` is one, including ``0**0``, infinity, and quiet NaN bases.
* ``1**b`` is one, including quiet NaN exponents.
* Other quiet NaNs propagate. Signaling NaNs raise invalid and propagate
  even through the preceding identities.
* Negative finite nonzero bases with nonintegral finite exponents return
  NaN and raise invalid.
* Zero to a negative finite exponent returns infinity and raises divide
  by zero. Zero to a positive exponent returns zero.
* Infinite bases to positive finite exponents return infinity; negative
  finite exponents return zero.
* Zero and infinity results have a negative sign only for negative bases
  and odd integral exponents. In particular, ``(-0)**0.5`` is positive zero.
* For infinite exponents, the result depends on whether ``abs(a)`` is
  below, equal to, or above one. ``(-1)**(+/-Inf)`` is one.

The general path bounds the exponential argument before multiplying very
large exponents, avoiding a spurious multiplication overflow when the
requested power should underflow instead.

Reproducibility and testing
---------------------------

Use round-to-nearest arithmetic, gradual underflow, preserved parentheses,
and disabled floating-point contraction. Unsafe reassociation, approximate
square roots, and flush-to-zero modes are outside the reproducibility
contract. NaN payload bits are not a portable result guarantee.

The existing intrinsic unit-test driver registers ``MOM_pow_repro_tests``.
The timing driver is ``time_MOM_intrinsic_functions``. With configured
testing builds, focused checks are::

  make -C .testing/build/unit test_MOM_intrinsic_functions
  .testing/build/unit/test_MOM_intrinsic_functions
  make -C .testing/build/timing time_MOM_intrinsic_functions
  .testing/build/timing/time_MOM_intrinsic_functions

This addition introduces a new interface without changing existing model
call sites or numerical answers. Replacing existing powers in model code
is a separate, potentially answer-changing operation.
