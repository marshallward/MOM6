! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Reproducible real powers with IEEE-style special values.
submodule (MOM_intrinsic_functions) MOM_pow

use MOM_log_data_n128, only : inv_c, log_c, binbits
use MOM_pow_log_data, only : log_c_tail

implicit none

contains

!> Evaluate a real power in a fixed arithmetic order.
!! Zero exponents (including 0**0) and unit bases return one, even with quiet
!! NaNs. Signaling NaNs propagate through arithmetic. Negative finite bases
!! require integral exponents; otherwise invalid is raised and NaN returned.
!! Signed zero and infinity results are negative only for negative bases and
!! odd integral exponents. Zero to a negative finite power raises divide by zero.
!! Squares use multiplication; +/-2**(-n), n=1,2,3, use nested sqrt.
!! Repeated roots need not be correctly rounded as a single power operation.
module procedure pow_repro
  integer(kind=int_rsize) :: ab, bb, ax, bx
    ! Input bit patterns and their unsigned magnitudes
  logical :: odd_power
    ! Whether b is an odd integer
  real :: base, power, log_base, log_tail, product, product_tail
    ! Positive base, exponent magnitude, logarithm and product high/low parts [nondim]
  real, parameter :: exp_limit = 2. * maxexponent(real_mold) * ln2
    ! Exponential argument safely beyond overflow and underflow limits [nondim]

  ab = transfer(a, int_mold)
  bb = transfer(b, int_mold)
  ax = iand(ab, not(signmask))
  bx = iand(bb, not(signmask))

  ! Quiet signaling NaNs before identities, without comparing NaNs as reals.
  if (ax > expmask) then
    if (.not. btest(ax, expbit - 1)) then
      p = a + b
      return
    endif
  endif
  if (bx > expmask) then
    if (.not. btest(bx, expbit - 1)) then
      p = a + b
      return
    endif
  endif
  if (bx == 0_int_rsize .or. ab == transfer(1., int_mold)) then
    p = 1.
    return
  endif
  if (ax > expmask .or. bx > expmask) then
    p = a + b
    return
  endif

  base = abs(a)
  if (bx == expmask) then
    p = 1.
    if (base /= 1.) then
      p = 0.
      if ((base > 1.) .eqv. (b > 0.)) p = transfer(expmask, real_mold)
    endif
    return
  endif

  ! No integer conversion: even the largest finite b is safe here. Beyond
  ! 2**digits(real_mold), all representable exponents are even integers.
  odd_power = .false.
  if (btest(ab, signbit)) odd_power = mod(abs(b), 2.) == 1.
  if (ax == 0_int_rsize .or. ax == expmask) then
    p = base
    if (b < 0.) p = 1. / base
    if (odd_power) p = -p
    return
  endif
  if (a < 0.) then
    if (b /= aint(b)) then
      p = (a - a) / (a - a)
      return
    endif
  endif

  power = abs(b)
  if (power == 1.) then
    p = base
    if (b < 0.) p = 1. / base
  else if (b == 2.) then
    p = base * base
  else if (power == 0.5 .or. power == 0.25 .or. power == 0.125) then
    p = sqrt(base)
    if (power <= 0.25) p = sqrt(p)
    if (power == 0.125) p = sqrt(p)
    if (b < 0.) p = 1. / p
  else
    call pow_log_extended(base, log_base, log_tail)
    ! Avoid overflow in b*log(base), which would otherwise give spurious
    ! overflow for a power that should instead underflow.
    if (power > 1.) then
      if (abs(log_base) > exp_limit / power) then
        p = exp_repro(sign(exp_limit, b) * sign(1., log_base))
        if (odd_power) p = -p
        return
      endif
    endif
    call pow_two_product(b, log_base, product, product_tail)
    product_tail = product_tail + b * log_tail
    p = exp_repro_extended(product, product_tail)
  endif
  if (odd_power) p = -p
end procedure pow_repro

!> Positive finite logarithm with a low-order residual for power evaluation.
!! Table reduction follows log_repro, but retains the product residual and
!! log(c) table tail. A separate near-one path avoids cancellation of table
!! constants when b is large and log(a) is small.
pure subroutine pow_log_extended(x, hi, lo)
  real, intent(in) :: x !< Positive finite input [nondim]
  real, intent(out) :: hi, lo !< High/low logarithm [nondim]
  real, parameter :: log2_hi = 0.693147180559945309417232121458176568
    !< Rounded natural logarithm of two [nondim]
  real, parameter :: log2_lo = 2.319046813846299558417771099e-17
    !< Residual of the rounded natural logarithm of two [nondim]
  integer(kind=int_rsize) :: xb, frac
    ! Normalized input and fraction bit patterns
  integer :: k, idiv
    ! Binary exponent and table index
  real :: z, r, rlo, r2, r2lo, tail, t, tlo, w, wlo
    ! Reduced input, powers, polynomial tail and compensated sums [nondim]

  if (abs(x - 1.) < 0.00390625) then
    r = x - 1. ! Exact subtraction near one.
    rlo = 0.
    t = 0.
    tlo = 0.
  else
    xb = transfer(x, int_mold)
    if (xb < expstep) then
      ! Keep intrinsic normalization for the uncommon subnormal case.
      k = exponent(x) - 1
      z = 2. * fraction(x)
      xb = transfer(z, int_mold)
    else
      ! Extract normal fields directly: exponent() and fraction() can each
      ! become a separate frexp library call even in optimized builds.
      k = int(shiftr(xb, expbit)) - expbias
      z = transfer(ior(iand(xb, expstep - 1), transfer(1., int_mold)), real_mold)
    endif
    frac = iand(xb, expstep - 1)
    idiv = int(shiftr(frac, expbit - binbits))
    call pow_two_product(z, inv_c(idiv), w, rlo)
    r = w - 1. ! Exact by Sterbenz's lemma; rlo retains product rounding.
    call pow_two_product(real(k), log2_hi, t, tlo)
    tlo = tlo + (real(k) * log2_lo + log_c_tail(idiv))
    call pow_two_sum(t, log_c(idiv), w, wlo)
    t = w
    tlo = tlo + wlo
  endif

  ! log1p(r+rlo) = r - r*r/2 + higher terms + rlo/(1+r).
  ! |r| < 1/256; the degree-10 Taylor remainder is below 3e-28.
  ! Compensating the leading two terms makes ordinary binary64 evaluation
  ! sufficient for the small cubic-and-higher tail.
  call pow_two_product(r, r, r2, r2lo)
  tail = (r * r2) * ((1./3.) + r * ((-1./4.) + r * ((1./5.) + &
      r * ((-1./6.) + r * ((1./7.) + r * ((-1./8.) + r * ((1./9.) - r / 10.)))))))
  call pow_two_sum(r, -0.5 * r2, w, wlo)
  wlo = wlo + ((tail - 0.5 * r2lo) + rlo / (1. + r))
  call pow_two_sum(t, w, hi, lo)
  lo = lo + (tlo + wlo)
  ! Normalize the pair before multiplying by a potentially large exponent.
  call pow_two_sum(hi, lo, w, wlo)
  hi = w
  lo = wlo
end subroutine pow_log_extended

!> Error-free sum of two finite values, with no magnitude-order assumption.
pure subroutine pow_two_sum(x, y, hi, lo)
  real, intent(in) :: x, y !< Summands [nondim]
  real, intent(out) :: hi, lo !< Sum and rounding residual [nondim]
  real :: v !< Portion of the rounded sum attributed to y [nondim]

  hi = x + y
  v = hi - x
  lo = (x - (hi - v)) + (y - v)
end subroutine pow_two_sum

!> Compensated product using bit-split operands instead of hardware FMA.
!! Truncating the low 27 fraction bits avoids splitter multiplication overflow.
!! The residual is accurate to second order in epsilon for normal products.
pure subroutine pow_two_product(x, y, hi, lo)
  real, intent(in) :: x, y !< Factors with a finite product [nondim]
  real, intent(out) :: hi, lo !< Product and rounding residual [nondim]
  integer(kind=int_rsize), parameter :: split_mask = -(2_int_rsize**27)
    !< Mask retaining the leading 26 significand bits
  real :: xhi, xlo, yhi, ylo !< High and low factors [nondim]

  xhi = transfer(iand(transfer(x, int_mold), split_mask), real_mold)
  yhi = transfer(iand(transfer(y, int_mold), split_mask), real_mold)
  xlo = x - xhi
  ylo = y - yhi
  hi = x * y
  lo = (((xhi * yhi - hi) + xhi * ylo) + xlo * yhi) + xlo * ylo
end subroutine pow_two_product

end submodule MOM_pow
