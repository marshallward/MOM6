! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Bitwise reproducible complementary error function.
!! Rational approximations and coefficients are from W. J. Cody, "Rational
!! Chebyshev approximations for the error function", Math. Comp. 23 (1969),
!! 631-637, as tabulated in https://netlib.org/specfun/erf (CALERF, 1990).
!! Evaluation uses compensated rational functions and exp_repro(), with a split
!! square and delayed underflow.  This is a roundoff-accuracy approximation,
!! not a correctly rounded implementation.  The accuracy target is binary64
!! default real; exp_repro() range reduction is also tuned for that precision.
!! See the interface for FP requirements.
submodule (MOM_intrinsic_functions) MOM_erfc

use, intrinsic :: iso_fortran_env, only : int32

implicit none

integer, parameter :: erfc_int_kind = merge(int64, int32, storage_size(real_mold) > 32)
  !< Integer kind matching default real
integer(erfc_int_kind), parameter :: erfc_int_mold = 0
  !< Integer mold for bit transfers
integer(erfc_int_kind), parameter :: infinity_bits = ishft(2_erfc_int_kind**expwidth - 1, expbit)
  !< Positive infinity bit pattern
!> Compensated working value; all components have the same nondimensional units.
type :: erfc_pair
  real :: hi !< Leading part [nondim]
  real :: lo !< Rounding residual [nondim]
end type erfc_pair

contains

!> Evaluate erfc without cancellation in the positive tail.
module procedure erfc_repro
  real, parameter :: inv_sqrt_pi = 0.564189583547756286948079451560772586
    !< Reciprocal square root of pi [nondim]
  real, parameter :: ln2 = 0.693147180559945309417232121458176568
    !< Natural logarithm of two [nondim]
  real, parameter :: ln2_hi = 0.69314718036912381649017333984375
    !< Leading part of ln(2), as in exp_repro() [nondim]
  real, parameter :: ln2_lo = 1.90821492927058770002e-10
    !< Trailing part of ln(2) [nondim]
  integer, parameter :: tail_scale = maxexponent(real_mold) / 2
    !< Power-of-two bias keeping tail exponentials and residuals normal
  real, parameter :: zero_bound = sqrt(real(digits(real_mold) - minexponent(real_mold) + 2) * ln2)
    !< Conservative positive bound beyond which erfc rounds to zero [nondim]
  real, parameter :: two_bound = sqrt(real(digits(real_mold)) * ln2)
    !< Conservative magnitude beyond which negative inputs round to two [nondim]
  real, parameter :: small_num(5) = [ &
      3.16112374387056560, 1.13864154151050156e2, 3.77485237685302021e2, &
      3.20937758913846947e3, 1.85777706184603153e-1 ]
    !< Numerator coefficients for erf(x)/x near zero [nondim]
  real, parameter :: small_den(4) = [ &
      2.36012909523441209e1, 2.44024637934444173e2, &
      1.28261652607737228e3, 2.84423683343917062e3 ]
    !< Denominator coefficients for erf(x)/x near zero [nondim]
  real, parameter :: mid_num(9) = [ &
      5.64188496988670089e-1, 8.88314979438837594, 6.61191906371416295e1, &
      2.98635138197400131e2, 8.81952221241769090e2, 1.71204761263407058e3, &
      2.05107837782607147e3, 1.23033935479799725e3, 2.15311535474403846e-8 ]
    !< Numerator coefficients for exp(x*x)*erfc(x), 0.46875 < x <= 4 [nondim]
  real, parameter :: mid_den(8) = [ &
      1.57449261107098347e1, 1.17693950891312499e2, 5.37181101862009858e2, &
      1.62138957456669019e3, 3.29079923573345963e3, 4.36261909014324716e3, &
      3.43936767414372164e3, 1.23033935480374942e3 ]
    !< Denominator coefficients for the middle interval [nondim]
  real, parameter :: tail_num(6) = [ &
      3.05326634961232344e-1, 3.60344899949804439e-1, 1.25781726111229246e-1, &
      1.60837851487422766e-2, 6.58749161529837803e-4, 1.63153871373020978e-2 ]
    !< Numerator coefficients for the asymptotic correction, x > 4 [nondim]
  real, parameter :: tail_den(5) = [ &
      2.56852019228982242, 1.87295284992346047, 5.27905102951428412e-1, &
      6.05183413124413191e-2, 2.33520497626869185e-3 ]
    !< Denominator coefficients for the asymptotic correction [nondim]
  real :: y, delta, exponential, correction
    ! Absolute input, exponential argument residual, exponential and its correction [nondim]
  type(erfc_pair) :: z, ratio, scaled_erfc, square, argument, factor, value
    ! Compensated polynomial argument, rational value, scaled complement and reconstruction
  integer(erfc_int_kind) :: xb, magnitude
    ! Input bits and unsigned magnitude bits
  integer :: rescale ! Binary exponent adjustment for the final result

  ! Classify by bits before comparisons or squaring, including signaling NaNs.
  xb = transfer(x, erfc_int_mold)
  magnitude = ibclr(xb, signbit)
  if (magnitude >= infinity_bits) then
    if (magnitude == infinity_bits) then
      a = merge(2., 0., btest(xb, signbit))
    else
      a = x + x
    endif
    return
  endif

  if (x <= -two_bound) then
    a = 2.
    return
  else if (x >= zero_bound) then
    a = 0.
    return
  endif

  y = abs(x)
  if (y < epsilon(real_mold) * 0.125) then
    ! Both signed zeros and tiny inputs round to one; avoid squaring subnormals.
    a = 1.
    return
  else if (y <= 0.46875) then
    z = erfc_multiply(erfc_pair(y, 0.), erfc_pair(y, 0.))
    ratio = erfc_rational(z, small_num, small_den)
    value = erfc_multiply(erfc_pair(-x, 0.), ratio)
    value = erfc_add(erfc_pair(1., 0.), value)
    a = value%hi + value%lo
    return
  else if (y <= 4.) then
    scaled_erfc = erfc_rational(erfc_pair(y, 0.), mid_num, mid_den)
  else
    square = erfc_multiply(erfc_pair(y, 0.), erfc_pair(y, 0.))
    z = erfc_divide(erfc_pair(1., 0.), square)
    ratio = erfc_rational(z, tail_num, tail_den)
    value = erfc_multiply(z, ratio)
    value = erfc_add(erfc_pair(inv_sqrt_pi, 0.), erfc_pair(-value%hi, -value%lo))
    scaled_erfc = erfc_divide(value, erfc_pair(y, 0.))
  endif

  ! Keep the rounded square and its residual.  The low part of the exponential
  ! argument must not be added back before calling exp_repro(): its loss would
  ! amplify relative error by O(x*x).  In the far tail, bias the exponent by an
  ! accurately split multiple of ln(2) so all compensated products stay normal.
  square = erfc_multiply(erfc_pair(y, 0.), erfc_pair(y, 0.))
  argument = erfc_pair(-square%hi, -square%lo)
  rescale = 0
  if (square%hi > -0.5 * real(minexponent(real_mold)) * ln2) then
    argument = erfc_add(argument, erfc_pair(real(tail_scale) * ln2_hi, real(tail_scale) * ln2_lo))
    rescale = -tail_scale
  endif

  ! exp(argument%lo) is almost one.  Preserve expm1(argument%lo) as a low
  ! component instead of rounding it away in a second exponential call.
  ! |delta| < 2^-43 in binary64; the omitted cubic is below 2^-129.
  delta = argument%lo
  correction = 0.5 * (delta * delta)
  factor = erfc_add(erfc_pair(1., 0.), erfc_pair(delta, correction))
  value = erfc_multiply(scaled_erfc, factor)
  exponential = exp_repro(argument%hi)
  value = erfc_multiply(value, erfc_pair(exponential, 0.))
  ! Scaling is exact for normal outputs.  Subnormals round only here, although
  ! rounding hi+lo first can introduce double rounding at a subnormal midpoint.
  a = scale(value%hi + value%lo, rescale)
  if (x < 0.) a = 2. - a
end procedure erfc_repro

!> Add two compensated values, retaining the leading addition's roundoff.
pure function erfc_add(a, b) result(c)
  type(erfc_pair), intent(in) :: a, b !< Operands
  type(erfc_pair) :: c !< Compensated sum
  real :: s, t, v ! Leading sum, residual, and recovered addend [nondim]

  s = a%hi + b%hi
  v = s - a%hi
  t = ((a%hi - (s - v)) + (b%hi - v)) + (a%lo + b%lo)
  c%hi = s + t
  c%lo = (s - c%hi) + t
end function erfc_add

!> Multiply compensated values using Dekker splitting, without an FMA.
!! Operands on the erfc paths are bounded away from overflow; the exponential
!! is scaled to keep these products and their residuals normal in binary64.
pure function erfc_multiply(a, b) result(c)
  type(erfc_pair), intent(in) :: a, b !< Operands
  type(erfc_pair) :: c !< Compensated product
  real, parameter :: splitter = scale(1., (digits(real_mold)+1)/2) + 1.
    !< Exact splitter for Dekker's product residual [nondim]
  real :: ah, al, bh, bl, t, p, e ! Split operands, product, and residual [nondim]

  t = splitter * a%hi
  ah = t - (t - a%hi)
  al = a%hi - ah
  t = splitter * b%hi
  bh = t - (t - b%hi)
  bl = b%hi - bh
  p = a%hi * b%hi
  e = (((ah * bh - p) + ah * bl) + al * bh) + al * bl
  e = e + (a%hi * b%lo + a%lo * b%hi)
  c%hi = p + e
  c%lo = (p - c%hi) + e
end function erfc_multiply

!> Correct a quotient using the compensated division residual.
pure function erfc_divide(a, b) result(c)
  type(erfc_pair), intent(in) :: a, b !< Numerator and nonzero denominator
  type(erfc_pair) :: c, product, residual !< Quotient and intermediate values
  real :: q, r ! Leading quotient and correction [nondim]

  q = a%hi / b%hi
  product = erfc_multiply(erfc_pair(q, 0.), b)
  residual = erfc_add(a, erfc_pair(-product%hi, -product%lo))
  r = (residual%hi + residual%lo) / b%hi
  c = erfc_add(erfc_pair(q, 0.), erfc_pair(r, 0.))
end function erfc_divide

!> Evaluate a Cody rational function with compensated Horner steps and division.
!! The last numerator coefficient is leading; the others follow in descending
!! order.  The denominator is monic with its remaining coefficients descending.
pure function erfc_rational(z, cn, cd) result(r)
  type(erfc_pair), intent(in) :: z !< Polynomial argument
  real, intent(in) :: cn(:), cd(:) !< Numerator and denominator coefficients [nondim]
  type(erfc_pair) :: r, num, den !< Rational approximation and polynomial values
  integer :: i ! Coefficient index

  num = erfc_pair(cn(size(cn)), 0.)
  den = erfc_pair(1., 0.)
  do i=1,size(cd)
    num = erfc_add(erfc_multiply(num, z), erfc_pair(cn(i), 0.))
    den = erfc_add(erfc_multiply(den, z), erfc_pair(cd(i), 0.))
  enddo
  r = erfc_divide(num, den)
end function erfc_rational

end submodule MOM_erfc
