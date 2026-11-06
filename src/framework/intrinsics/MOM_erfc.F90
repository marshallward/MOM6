! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Bitwise reproducible complementary error function.
!!
!! Near zero, evaluate erf(x)/x as a mixed-precision Taylor polynomial.
!! Elsewhere, approximate exp(x*x)*erfc(x) with piecewise Chebyshev
!! polynomials, and retain hi/lo components through the exponential and final
!! rounding.  Results near rounding midpoints use a three-component fallback
!! with longer polynomials.  Small high-order terms use ordinary arithmetic;
!! significant stages retain compensation. The exponential shares exp_repro()'s
!! range-reduction table, but retains extra precision.
!!
!! Coefficients are generated with MPFR; see .testing/tools/erfc/generate.c.
!! The accuracy target is binary64 default real.  Correct rounding is tested
!! against MPFR, but has not been proved for every binary64 input.
!! See the interface for FP requirements.
submodule (MOM_intrinsic_functions) MOM_erfc

use, intrinsic :: iso_fortran_env, only : int32
use MOM_erfc_data, only : erfcx_hi, erfcx_lo, erf_hi, erf_lo, exp_hi, exp_lo
use MOM_erfc_data, only : erfcx_tail, erf_tail, exp_tail, ln2_parts
use MOM_exp_data_n128, only : ndiv, idiv_scale_lookup, idiv_residual_lookup

implicit none

integer, parameter :: erfc_int_kind = merge(int64, int32, storage_size(real_mold) > 32)
  !< Integer kind matching default real
integer(erfc_int_kind), parameter :: erfc_int_mold = 0
  !< Integer mold for bit transfers
integer(erfc_int_kind), parameter :: infinity_bits = ishft(2_erfc_int_kind**expwidth - 1, expbit)
  !< Positive infinity bit pattern

!> Compensated working value; all components have the same nondimensional
!! units.
type :: erfc_pair
  real :: hi
    !< Leading part [nondim]
  real :: lo
    !< Rounding residual [nondim]
end type erfc_pair

!> Three-component expansion used only when the fast result is near a midpoint.
type :: erfc_wide
  real :: v(3)
    !< Components in decreasing magnitude [nondim]
end type erfc_wide

contains

!> Evaluate erfc without cancellation or premature rounding in the positive
!! tail.
module procedure erfc_repro
  real, parameter :: ln2 = 0.693147180559945309417232121458176568
    !< Natural logarithm of two [nondim]
  real, parameter :: zero_bound = sqrt(real(digits(real_mold) - minexponent(real_mold) + 2) * ln2)
    !< Conservative positive bound beyond which erfc rounds to zero [nondim]
  real, parameter :: two_bound = sqrt(real(digits(real_mold)) * ln2)
    !< Conservative magnitude beyond which negative inputs round to two [nondim]
  real :: y, t, error, lower, upper ! Input, reduced argument, error budget, and rounded bounds [nondim]
  type(erfc_pair) :: square, scaled_erfc, exponential, value
    ! Compensated square, scaled complement, exponential mantissa, and result
  integer(erfc_int_kind) :: xb, magnitude
    ! Input bits and unsigned magnitude bits
  integer :: rescale, bin ! Binary exponent adjustment and Chebyshev interval

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
    square = erfc_multiply(erfc_pair(y, 0.), erfc_pair(y, 0.))
    value = erfc_polynomial(square, erf_hi(0:15), erf_lo(0:15), 6)
    value = erfc_multiply(erfc_pair(-x, 0.), value)
    ! Working error budget includes mixed-precision evaluation and truncation.
    ! Compare rounded interval endpoints before accepting a near-midpoint result.
    error = scale(abs(value%hi), -72)
    lower = erfc_round_sum(1., erfc_pair(value%hi, value%lo-error))
    upper = erfc_round_sum(1., erfc_pair(value%hi, value%lo+error))
    a = lower
    if (lower /= upper) a = erfc_accurate(x)
    return
  endif

  ! Map each interval exactly onto [-1,1].  The interval boundaries, centers,
  ! and radii are binary fractions, and the subtraction is exact by Sterbenz.
  if (y < 1.) then
    bin = 0
    t = 2. * y - 1.
  else
    bin = exponent(y)
    t = scale(y, 2-bin) - 3.
  endif
  scaled_erfc = erfc_chebyshev(t, bin)
  square = erfc_multiply(erfc_pair(y, 0.), erfc_pair(y, 0.))
  call erfc_exp(erfc_pair(-square%hi, -square%lo), exponential, rescale)
  value = erfc_multiply(scaled_erfc, exponential)
  if (x < 0.) then
    ! Preserve the low part until after reflection, avoiding double rounding.
    value = erfc_pair(-scale(value%hi, rescale), -scale(value%lo, rescale))
    error = scale(abs(value%hi), -72)
    lower = erfc_round_sum(2., erfc_pair(value%hi, value%lo-error))
    upper = erfc_round_sum(2., erfc_pair(value%hi, value%lo+error))
  else
    error = scale(abs(value%hi), -72)
    lower = erfc_round(erfc_pair(value%hi, value%lo-error), rescale)
    upper = erfc_round(erfc_pair(value%hi, value%lo+error), rescale)
  endif
  a = lower
  if (lower /= upper) a = erfc_accurate(x)
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


!> Multiply a scalar and a compensated value, deferring normalization to
!! erfc_add.  Used only for bounded Clenshaw operands followed immediately by
!! an addition.
pure function erfc_scalar_product(a, b) result(c)
  real, intent(in) :: a !< Scalar operand [nondim]
  type(erfc_pair), intent(in) :: b !< Compensated operand
  type(erfc_pair) :: c !< Unnormalized product
  real, parameter :: splitter = scale(1., (digits(real_mold)+1)/2) + 1.
    !< Exact splitter for Dekker's product residual [nondim]
  real :: ah, al, bh, bl, t, p, e ! Split operands, product, and residual [nondim]

  t = splitter * a
  ah = t - (t - a)
  al = a - ah
  t = splitter * b%hi
  bh = t - (t - b%hi)
  bl = b%hi - bh
  p = a * b%hi
  e = (((ah * bh - p) + ah * bl) + al * bh) + al * bl
  c%hi = p
  c%lo = e + a * b%lo
end function erfc_scalar_product


!> Add compensated values using FastTwoSum when abs(a%hi) >= abs(b%hi).
pure function erfc_dominant_add(a, b) result(c)
  type(erfc_pair), intent(in) :: a, b !< Magnitude-ordered operands
  type(erfc_pair) :: c !< Compensated sum
  real :: s, e ! Leading sum and residual [nondim]

  s = a%hi + b%hi
  e = ((a%hi - s) + b%hi) + (a%lo + b%lo)
  c%hi = s + e
  c%lo = (s - c%hi) + e
end function erfc_dominant_add


!> Evaluate small high-order terms ordinarily, then use compensated Horner
!! steps.
pure function erfc_polynomial(z, ch, cl, cutoff) result(r)
  type(erfc_pair), intent(in) :: z
    !< Polynomial argument
  real, intent(in) :: ch(0:), cl(0:)
    !< Leading coefficients and residuals [nondim]
  integer, intent(in) :: cutoff
    !< First ordinary stage; lower indices use compensation

  type(erfc_pair) :: r
    !< Polynomial value
  real :: high
    ! Ordinary high-order polynomial [nondim]
  integer :: i
    ! Coefficient index

  high = ch(size(ch)-1)
  do i=size(ch)-2,cutoff,-1
    high = ch(i) + z%hi * high
  enddo
  r = erfc_pair(high, 0.)
  do i=cutoff-1,0,-1
    r = erfc_add(erfc_multiply(r, z), erfc_pair(ch(i), cl(i)))
  enddo
end function erfc_polynomial


!> Mixed-precision Clenshaw evaluation of the scaled complementary error function.
pure function erfc_chebyshev(t, bin) result(r)
  real, intent(in) :: t
    !< Reduced input in [-1,1] [nondim]
  integer, intent(in) :: bin
    !< Table interval
  type(erfc_pair) :: r, b1, b2, next
    !< Polynomial value and Clenshaw recurrence

  integer, parameter :: degree(0:5) = [24, 22, 26, 29, 31, 31]
    !< Degrees placing discarded coefficient tails below 2^-80 relative
  real :: h1, h2, hnext
    ! Ordinary high-order recurrence [nondim]
  integer :: i
    ! Coefficient index

  h1 = 0.
  h2 = 0.

  do i=degree(bin),12,-1
    hnext = ((2.*t)*h1 - h2) + erfcx_hi(i,bin)
    h2 = h1
    h1 = hnext
  enddo

  b1 = erfc_pair(h1, 0.)
  b2 = erfc_pair(h2, 0.)

  do i=11,1,-1
    next = erfc_scalar_product(2.*t, b1)
    next = erfc_add(next, erfc_pair(-b2%hi, -b2%lo))
    ! Coefficient-tail bounds give abs(next%hi)/abs(coefficient) < 0.541.
    next = erfc_dominant_add(erfc_pair(erfcx_hi(i,bin), erfcx_lo(i,bin)), next)
    b2 = b1
    b1 = next
  enddo

  r = erfc_scalar_product(t, b1)
  r = erfc_add(r, erfc_pair(-b2%hi, -b2%lo))
  r = erfc_add(r, erfc_pair(erfcx_hi(0,bin), erfcx_lo(0,bin)))
end function erfc_chebyshev


!> Compute exp(x) = (a%hi+a%lo)*2**k without prematurely rounding to one double.
!! The input is the exact compensated negative square, bounded by the erfc cutoff.
pure subroutine erfc_exp(x, a, k)
  type(erfc_pair), intent(in) :: x
    !< Compensated negative square
  type(erfc_pair), intent(out) :: a
    !< Exponential mantissa
  integer, intent(out) :: k
    !< Binary exponent

  real, parameter :: inv_ln2 = 1.44269504088896340735992468100189214
    !< Reciprocal ln(2) [nondim]
  real, parameter :: ln2_hi = 0.69314718036912381649017333984375
    !< Leading part of ln(2), whose product with the reduction integer is exact [nondim]
  real, parameter :: ln2_lo = 1.90821492927058770002e-10
    !< Second part of ln(2) [nondim]
  real, parameter :: ln2_tail = 1.1612227229362532e-26
    !< Residual of the first two parts of ln(2) [nondim]

  type(erfc_pair) :: r, correction, table
    !< Reduced argument, reduction residual, and table value
  real :: z
    ! Reduction integer divided by the table size [nondim]
  integer :: n, j
    ! Reduction integer and table index

  n = nint(x%hi * (real(ndiv) * inv_ln2))
  j = modulo(n, ndiv)
  k = (n-j) / ndiv
  z = real(n) / real(ndiv)
  r = erfc_add(x, erfc_pair(-z * ln2_hi, 0.))
  correction = erfc_multiply(erfc_pair(-z, 0.), erfc_pair(ln2_lo, ln2_tail))
  r = erfc_add(r, correction)
  a = erfc_polynomial(r, exp_hi(0:7), exp_lo(0:7), 3)
  table = erfc_multiply(erfc_pair(idiv_scale_lookup(j), 0.), &
      erfc_pair(1., idiv_residual_lookup(j)))
  a = erfc_multiply(a, table)
end subroutine erfc_exp


!> Round base + a hi/lo correction without losing a tiny midpoint displacement.
!! Returning a normalized hi/lo sum first would discard the third component
!! when base is much larger than the correction (especially for tiny inputs).
pure function erfc_round_sum(base, correction, tail) result(a)
  real, intent(in) :: base
    !< Exactly representable offset [nondim]
  type(erfc_pair), intent(in) :: correction
    !< Signed correction
  real, optional, intent(in) :: tail
    !< Third component of the correction [nondim]

  real :: a
    !< Rounded sum [nondim]
  real :: h, v, residual, direction, neighbor, half_gap, last
    ! Summation and rounding terms [nondim]
  type(erfc_pair) :: low, distance
    !< Residual and distance past the rounding midpoint

  h = base + correction%hi
  v = h - base
  residual = (base - (h - v)) + (correction%hi - v)
  last = 0.
  if (present(tail)) last = tail
  low = erfc_add(erfc_pair(residual, 0.), erfc_pair(correction%lo, last))
  a = h
  if (low%hi == 0.) return
  direction = sign(1., low%hi)
  neighbor = nearest(h, direction)
  half_gap = 0.5 * abs(neighbor - h)
  distance = erfc_add(erfc_pair(direction * residual - half_gap, 0.), &
      erfc_pair(direction * correction%lo, direction * last))
  if (distance%hi > 0.) then
    a = neighbor
  elseif (distance%hi == 0.) then
    if (distance%lo > 0. .or. (distance%lo == 0. .and. btest(transfer(h, erfc_int_mold),0))) a = neighbor
  endif
end function erfc_round_sum


!> Round a positive compensated value after binary scaling, including subnormals.
pure function erfc_round(value, k) result(a)
  type(erfc_pair), intent(in) :: value
    !< Positive mantissa
  integer, intent(in) :: k
    !< Binary exponent
  real :: a
    !< Rounded and scaled result [nondim]
  integer, parameter :: quantum_exp = minexponent(real_mold) - digits(real_mold)
    !< Exponent of the least positive subnormal
  integer(erfc_int_kind) :: n
    ! Integer number of subnormal quanta
  type(erfc_pair) :: units, fraction
    !< Result in subnormal quanta and fractional remainder

  if (exponent(value%hi) + k > minexponent(real_mold)) then
    a = scale(value%hi + value%lo, k)
  else
    ! Round to integer multiples of the subnormal quantum BEFORE scaling.
    ! Rounding hi+lo first would sometimes double-round at subnormal midpoints.
    units = erfc_pair(scale(value%hi, k-quantum_exp), scale(value%lo, k-quantum_exp))
    n = int(units%hi, erfc_int_kind)
    fraction = erfc_add(units, erfc_pair(-real(n), 0.))
    if (fraction%hi > 0.5) then
      n = n + 1
    else if (fraction%hi == 0.5) then
      if (fraction%lo > 0. .or. (fraction%lo == 0. .and. btest(n,0))) n = n + 1
    endif
    a = scale(real(n), quantum_exp)
  endif
end function erfc_round


!> Error-free sum of two ordinary reals, with no multiply-add contraction.
pure subroutine erfc_two_sum(a, b, s, e)
  real, intent(in) :: a, b
    !< Addends [nondim]
  real, intent(out) :: s, e
    !< Rounded sum and exact residual [nondim]

  real :: v
    ! Recovered addend [nondim]

  s = a + b
  v = s - a
  e = (a - (s - v)) + (b - v)
end subroutine erfc_two_sum


!> Accumulate terms into three components, keeping two levels of sum residuals.
pure function erfc_accumulate(terms) result(r)
  real, intent(in) :: terms(:)
    !< Terms to accumulate [nondim]

  type(erfc_wide) :: r
    !< Three-component sum
  real :: h, l, t, s, e, f, u, v
    ! Accumulators and error-free sum outputs [nondim]
  integer :: i
    ! Term or normalization index

  h = 0.
  l = 0.
  t = 0.

  do i=1,size(terms)
    call erfc_two_sum(h, terms(i), s, e)
    h = s
    call erfc_two_sum(l, e, s, f)
    l = s
    t = t + f
  enddo

  do i=1,2
    call erfc_two_sum(l, t, s, e)
    call erfc_two_sum(h, s, u, v)
    h = u
    l = v
    t = e
  enddo

  r%v(:) = [h, l, t]
end function erfc_accumulate


!> Three-component addition for the rare accurate path.
pure function erfc_wide_add(a, b) result(r)
  type(erfc_wide), intent(in) :: a, b
    !< Operands
  type(erfc_wide) :: r
    !< Sum

  r = erfc_accumulate([a%v(1), b%v(1), a%v(2), b%v(2), a%v(3), b%v(3)])
end function erfc_wide_add


!> Three-component multiplication using exact two-component products.
pure function erfc_wide_multiply(a, b) result(r)
  type(erfc_wide), intent(in) :: a, b
    !< Operands
  type(erfc_wide) :: r
    !< Product

  type(erfc_pair) :: p
    ! Exact product of two ordinary components
  real :: terms(18)
    ! Products and their residuals [nondim]
  integer :: i, j, n
    ! Component and term indices

  n = 0

  do i=1,3
    do j=1,3
      p = erfc_multiply(erfc_pair(a%v(i), 0.), erfc_pair(b%v(j), 0.))
      terms(n+1:n+2) = [p%hi, p%lo]
      n = n + 2
    enddo
  enddo

  r = erfc_accumulate(terms(:))
end function erfc_wide_multiply


!> Three-component Horner evaluation including coefficient residuals.
pure function erfc_wide_polynomial(z, ch, cl, ct) result(r)
  type(erfc_wide), intent(in) :: z
    !< Argument
  real, intent(in) :: ch(0:), cl(0:), ct(0:)
    !< Three coefficient components [nondim]
  type(erfc_wide) :: r
    !< Polynomial value
  integer :: i
    ! Coefficient index

  i = size(ch)-1
  r = erfc_wide([ch(i), cl(i), ct(i)])

  do i=size(ch)-2,0,-1
    r = erfc_wide_add(erfc_wide_multiply(r, z), erfc_wide([ch(i), cl(i), ct(i)]))
  enddo
end function erfc_wide_polynomial


!> Higher-precision fallback, using three-component arithmetic and longer polynomials.
!! No hardware extended precision or platform transcendental functions are used.
pure function erfc_accurate(x) result(a)
  real, intent(in) :: x
    !< Finite input inside the saturation bounds [nondim]

  real :: a
    !< Rounded complementary error function [nondim]
  real :: y, t
    ! Absolute input and reduced Chebyshev argument [nondim]
  type(erfc_wide) :: square, value, b1, b2, next, r, exponential, units
    ! Compensated square, polynomial recurrence, exponential and rounding values
  integer :: bin, i, k
    ! Interval, coefficient index, and binary exponent
  integer(erfc_int_kind) :: n
    ! Integer count of subnormal quanta
  integer, parameter :: quantum_exp = minexponent(real_mold) - digits(real_mold)
    !< Exponent of the least positive subnormal
  type(erfc_pair) :: fraction
    ! Fractional distance beyond a subnormal midpoint

  y = abs(x)
  square = erfc_wide_multiply(erfc_wide([y, 0., 0.]), erfc_wide([y, 0., 0.]))

  if (y <= 0.46875) then
    value = erfc_wide_polynomial(square, erf_hi, erf_lo, erf_tail)
    value = erfc_wide_multiply(erfc_wide([-x, 0., 0.]), value)
    a = erfc_round_sum(1., erfc_pair(value%v(1), value%v(2)), value%v(3))
    return
  endif

  if (y < 1.) then
    bin = 0
    t = 2.*y - 1.
  else
    bin = exponent(y)
    t = scale(y, 2-bin) - 3.
  endif

  b1 = erfc_wide([0., 0., 0.])
  b2 = b1

  do i=ubound(erfcx_hi,1),1,-1
    next = erfc_wide_multiply(erfc_wide([2.*t, 0., 0.]), b1)
    next = erfc_wide_add(next, erfc_wide([-b2%v(1), -b2%v(2), -b2%v(3)]))
    next = erfc_wide_add(next, erfc_wide([erfcx_hi(i,bin), erfcx_lo(i,bin), erfcx_tail(i,bin)]))
    b2 = b1
    b1 = next
  enddo

  value = erfc_wide_multiply(erfc_wide([t, 0., 0.]), b1)
  value = erfc_wide_add(value, erfc_wide([-b2%v(1), -b2%v(2), -b2%v(3)]))
  value = erfc_wide_add(value, erfc_wide([erfcx_hi(0,bin), erfcx_lo(0,bin), erfcx_tail(0,bin)]))

  ! Independent full-ln(2) reduction avoids the fast path's two-component table.
  k = nint(-square%v(1) / ln2_parts(1))
  r = erfc_wide([-square%v(1), -square%v(2), -square%v(3)])
  next = erfc_wide_multiply(erfc_wide([-real(k), 0., 0.]), erfc_wide(ln2_parts(:)))
  r = erfc_wide_add(r, next)
  exponential = erfc_wide_polynomial(r, exp_hi, exp_lo, exp_tail)
  value = erfc_wide_multiply(value, exponential)

  if (x < 0.) then
    a = erfc_round_sum(2., erfc_pair(-scale(value%v(1), k), -scale(value%v(2), k)), -scale(value%v(3), k))
  else if (exponent(value%v(1)) + k > minexponent(real_mold)) then
    a = scale(erfc_round_sum(value%v(1), erfc_pair(value%v(2), value%v(3))), k)
  else
    units = erfc_wide([scale(value%v(1), k-quantum_exp), &
        scale(value%v(2), k-quantum_exp), scale(value%v(3), k-quantum_exp)])
    n = int(units%v(1), erfc_int_kind)
    fraction = erfc_add(erfc_pair((units%v(1)-real(n)) - 0.5, 0.), erfc_pair(units%v(2), units%v(3)))

    if (fraction%hi > 0.) then
      n = n + 1
    else if (fraction%hi == 0.) then
      if (fraction%lo > 0. .or. (fraction%lo == 0. .and. btest(n,0))) n = n + 1
    endif

    a = scale(real(n), quantum_exp)
  endif
end function erfc_accurate

end submodule MOM_erfc
