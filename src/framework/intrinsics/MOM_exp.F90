! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!! This submodule provides a bitwise reproducible implementation of exp().

#include "MOM_exp.h"

submodule (MOM_intrinsic_functions) MOM_exp

use, intrinsic :: iso_fortran_env, only : int32, int64, real64
use, intrinsic :: ieee_arithmetic, only : ieee_rint

implicit none

! Mold for transfers and numerical format queries
integer(kind=int64), parameter :: int64_mold = 0
real(kind=real64), parameter :: real64_mold = 0.

! Floating point layout
integer, parameter :: expbit = digits(real64_mold) - 1
  !< Position of lowest exponent bit (52 in double)
integer, parameter :: expwidth = storage_size(real64_mold) - expbit - 1
  !< Number of exponent bits (11 in double)
integer, parameter :: expbias = maxexponent(real64_mold) - 1
  !< Exponent bias (1023 in double)

! IEEE 754 special values
integer(int64), parameter :: pos_inf_bits = ishft(2_int64**expwidth - 1, expbit)
  !< IEEE 64-bit +Inf
integer(int64), parameter :: neg_inf_bits = ior(pos_inf_bits, ishft(-1_int64, 63))
  !< IEEE 64-bit -Inf

contains

!> Reproducible exponential function
!!
!! Compute exp(x) with bitwise reproducibility across platforms.
module procedure exp_repro
  ! ln2 estimates
  real(kind=real64), parameter :: ln2 &
      = transfer(int(z'3FE62E42FEFA39EF', int64), 0.0_real64)
    ! 0.693147180559945309417232...
  real(kind=real64), parameter :: I_ln2 &
      = transfer(int(z'3FF71547652B82FE', int64), real64_mold)
    ! 1 / ln2: 1.4426950408889634073599...

  ! The max and min x values between which x remains finite.
  real(kind=real64), parameter :: xmax &
      = real(maxexponent(real64_mold) + 1, real64) * ln2
  real(kind=real64), parameter :: xmin &
      = real(minexponent(real64_mold) - digits(real64_mold) - 1, real64) * ln2

  ! Double-double precision of ln2 used in range reduction (Cody-Waite)
  real(kind=real64), parameter :: ln2_hi &
      = transfer(int(z'3FE62E42FEE00000', int64), real64_mold)
    ! Upper 32 bits of ln2: 6.93147180369123816490e-01
  real(kind=real64), parameter :: ln2_lo &
      = transfer(int(z'3DEA39EF35793C76', int64), real64_mold)
    ! Lower precision bits of ln2: 1.90821492927058770002e-10

  real(kind=real64) :: xc
    ! x after being bounded by xmin and xmax
  real(kind=real64) :: x_ln2
    ! Cached value of x / ln2
  real(kind=real64) :: K
    ! nint(x / ln2), used to minimize r = x - K ln2
    ! NOTE: Stored as a real to avoid costly type conversions
  real(kind=real64) :: r
    ! Range-reduced input, r = x - K ln2
  real(kind=real64) :: e
    ! Exponent of range-reduced input, e = exp(r)

  integer(kind=int64) :: xb, Kb, eb
    ! Bit representations of x, K, e

  integer(int64) :: j, fb
    ! j-scaling for subnormal handling; fb is the scaling factor bits

  logical :: nonfinite
    ! True if input is a nonfinite float (+/-Inf, NaN)

  ! 1. Nonfinite handling
  ! ---------------------
  ! Nonfinites must be handled first to prevent their appearance in
  ! calculations, which may raise unwanted floating point signals.

  xb = transfer(x, int64_mold)
  nonfinite = iand(xb, pos_inf_bits) == pos_inf_bits

  if (nonfinite) then
    ! exp(-Inf) = 0, otherwise pass-through for +Inf and +/-NaN.
    ! Use nonfinite result `x + x = x` to trigger Invalid for signaled NaNs.
    a = merge(0._real64, x + x, xb == neg_inf_bits)
    return
  endif

  ! 2. Range Reduction
  ! ------------------
  ! Apply a range reduction of r = x - K ln2, so that exp(x) = 2**K exp(r).
  ! If K = nint(x / ln2) then r in [-1/2 ln2, 1/2 ln2] and exp(r) can be
  ! estimated by a highly accurate polynomial.

  ! First limit x to the precision-based numerical range of exp().
  ! This allows for safe detection of over/underflow in the subnormal handler.
  xc = min(max(x, xmin), xmax)

  ! Compute K = nint(x / ln2)
  x_ln2 = xc * I_ln2
  K = NEAREST_INT(x_ln2)

  ! NOTE: Performance may degrade if NEAREST_INT contains a function call.

  ! Since K ~ x/ln2, the terms in r = x - K ln2 will nearly cancel and there is
  ! some expected loss of precision.  For 64-bit reals, |K| can reach values of
  ! about 1000 before overflow, causing a loss of about 10 bits of precision.

  ! To compensate, we apply a Cody-Waite correction which separates ln2 into
  ! its upper 32 bits and a lower residual.
  r = (xc - K * ln2_hi) - K * ln2_lo

  ! NOTE: Aggressive optimizers may reduce this to x - K * (ln2_hi + ln2_lo)
  ! which is no better than x - K ln2.

  ! 3. Polynomial approximation
  ! ---------------------------
  ! Estimate exp(r) = 1 + r * P(r) where P(x) is a 10th order Remez "minimax"
  ! polynomial of (exp(x) - 1)/x.  This form is used to ensure that exp(0) = 1.

  e = exp_remez_horner_10(r)

  ! 4. Unscaling
  ! ------------
  ! Compute exp(x) = 2**K exp(r), an exact power-of-2 calculation.
  ! Adjust scaling to compensate for subnormal output.

  ! Determine if a subnormal j-scaling may be required.  The range is reduced
  ! from [-1023,+1022] to [-1020, +1020] to account for the range of exp(r).
  j = merge(1022_int64, 0_int64, K < -1020.0_real64) &
      + merge(-1022_int64, 0_int64, K > 1020.0_real64)

  ! Get the bit representation of exp(r).
  eb = transfer(e, int64_mold)

  ! Convert exponent K to integer form and apply both K and j-scaling to e.
  Kb = int(K, int64)
  eb = eb + ishft(Kb + j, expbit)

  ! Convert 2**(K+j) exp(r) to float, stored to `a`.
  ! For moderate values, j is zero and this completes scaling.
  a = transfer(eb, real64_mold)

  ! Undo the j-scaling correction as a multiplication.
  ! For subnormals, this reverts the j-scaling.
  ! For extreme values of K, this triggers floating point over/underflow.
  fb = ishft(expbias - j, expbit)
  a = a * transfer(fb, real64_mold)
end procedure exp_repro


!> Remez minimax polynomial for exp(x) over [-1/2 ln2, 1/2 ln2] in Horner form
!!
!! Approximate (exp(x) - 1) / x, then compute 1 + x * P(x).
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_horner_10(x) result(e)
  real(kind=real64), intent(in) :: x
    !< Input value; expected range is [-1/2 ln2, 1/2 ln2]
  real(kind=real64) :: e
    !< Approximation of exp(x)

  real(kind=real64) :: p

  ! fpminimax coefficients for (exp(x) - 1) / x on [-1/2 ln2, 1/2 ln2]
  real(kind=real64), parameter :: c(0:10) = [ &
      1.0_real64, &
      5.00000000000000555e-1_real64, &
      1.66666666666666074e-1_real64, &
      4.16666666665738844e-2_real64, &
      8.33333333337716447e-3_real64, &
      1.38888889322647565e-3_real64, &
      1.98412697469850182e-4_real64, &
      2.48015045964426143e-5_real64, &
      2.75573817985163132e-6_real64, &
      2.76262647076892519e-7_real64, &
      2.50621020021886396e-8_real64 &
  ]

  ! Horner evaluation: p(x) = c0 + x*(c1 + x*(c2 + x*(c3 + ...)))
  p = c(10)
  p = c(9) + x * p
  p = c(8) + x * p
  p = c(7) + x * p
  p = c(6) + x * p
  p = c(5) + x * p
  p = c(4) + x * p
  p = c(3) + x * p
  p = c(2) + x * p
  p = c(1) + x * p
  p = c(0) + x * p

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_horner_10


!> An optimized nearest-integer function for floating point reals.
!!
!! The value x is shifted from 2**K (1 + a) to 2**(52+K) (1.5 + a'), causing
!! any fractional bits to be rounded according to current IEEE setting.  In
!! almost all cases, this is nearest ties-to-even.
!!
!! The behavior of this function does not match nint() or anint().  The nint()
!! function always ties away from zero, e.g. nint(2.5) = 3.
!!
!! The +0.5 ensures that biased 2**52 exponent of negative numbers does not
!! drop to 2**51, which can cause half-value rounding.
!!
!! It is essential that compilers not reduce (x+b)-b to x.  This can typically
!! be ensured as long as parentheses are respected.  This is managed by the
!! ENABLE_FAST_RINT macro in exp_repro.h and assigned to NEAREST_INT().  If
!! unset, then ieee_rint() is used.
pure function fast_rint(x) result(n)
  real(real64), intent(in) :: x
    !< Real value to be rounded to the nearest integer
  real(real64) :: n
    !< Nearest integer to x, stored as a real

  real(kind=real64), parameter :: round_bias = 1.5_real64 * 2_int64**52
    !< Binary offset used to trigger rounding of fractional values

  n = (x + round_bias) - round_bias
end function fast_rint

end submodule MOM_exp
