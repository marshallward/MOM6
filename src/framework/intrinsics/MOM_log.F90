! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This submodule provides a bitwise reproducible implementation of log().

submodule (MOM_intrinsic_functions) MOM_log

use, intrinsic :: iso_fortran_env, only : int32, int64

implicit none

integer, parameter :: int_kind &
    = merge(int64, int32, storage_size(real_mold) > storage_size(0_int32))
  !< Integer kind with the same storage size as default real
integer(kind=int_kind), parameter :: int_mold = 0
  !< Integer mold value

! IEEE 754 masks
integer(kind=int_kind), parameter :: exp_mask = 2_int_kind**expwidth - 1_int_kind
  !< Mask for the biased exponent value
integer(kind=int_kind), parameter :: mant_mask = 2_int_kind**expbit - 1_int_kind
  !< Mask for the significand bits
integer(kind=int_kind), parameter :: sign_mask = ishft(-1_int_kind, signbit)
  !< Mask for the sign bit

contains

!> Reproducible natural logarithm function
!!
!! Compute log(x) with bitwise reproducibility across platforms.
module procedure log_repro
  real, parameter :: ln2_hi = 0.69314718036912381649017333984375
    !< Upper 32 bits of ln2: 6.93147180369123816490e-01 [nondim]
  real, parameter :: ln2_lo = 1.90821492927058770002e-10
    !< Lower precision bits of ln2: 1.90821492927058770002e-10 [nondim]
  real, parameter :: sqrt2 = 1.41421356237309504880168872420969808
    !< sqrt(2) [nondim]
  integer(kind=int_kind), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used to normalize subnormal inputs
  real, parameter :: scale_up = transfer(ishft(int(expbias, int_kind) + Kbias, expbit), real_mold)
    !< Exact power-of-two scale factor for subnormal inputs [nondim]

  integer(kind=int_kind) :: xb, mb
    ! Bit representations of x and its normalized significand
  integer(kind=int_kind) :: raw_exp
    ! Biased IEEE exponent field
  integer(kind=int_kind) :: K
    ! Binary exponent in x = 2**K m [nondim]
  real :: xs
    ! Input value, possibly scaled to normalize subnormal numbers [nondim]
  real :: m
    ! Significand of x, adjusted into [1/sqrt(2),sqrt(2)] [nondim]
  real :: y
    ! Reduced argument, y = (m - 1) / (m + 1) [nondim]
  real :: log_m
    ! Approximation to log(m) [nondim]

  xb = transfer(x, int_mold)

  ! Handle exceptional values before arithmetic range reduction.  This gives
  ! log(0) = -Inf with divide-by-zero, log(negative) = NaN with invalid,
  ! log(+Inf) = +Inf, and NaNs pass through with the usual signaling behavior.
  if (x == 0.) then
    a = -1. / abs(x)
    return
  endif

  if (iand(xb, sign_mask) /= 0_int_kind) then
    a = (x - x) / (x - x)
    return
  endif

  raw_exp = iand(ishft(xb, -expbit), exp_mask)

  if (raw_exp == exp_mask) then
    a = x + x
    return
  endif

  ! Range reduction: decompose x = 2**K m with m in [1,2).  Subnormals are
  ! first multiplied by an exact power of two, then compensated in K.
  if (raw_exp == 0_int_kind) then
    xs = x * scale_up
    xb = transfer(xs, int_mold)
    raw_exp = iand(ishft(xb, -expbit), exp_mask)
    K = raw_exp - int(expbias, int_kind) - Kbias
  else
    K = raw_exp - int(expbias, int_kind)
  endif

  mb = ior(iand(xb, mant_mask), ishft(int(expbias, int_kind), expbit))
  m = transfer(mb, real_mold)

  ! Keep m close to 1 so the polynomial only sees a small interval.
  ! Then log(x) = K ln2 + log(m), with m in [1/sqrt(2),sqrt(2)].
  if (m > sqrt2) then
    m = 0.5 * m
    K = K + 1_int_kind
  endif

  ! log(m) = 2 * atanh(y), y = (m - 1) / (m + 1).  This gives a symmetric
  ! reduced range, |y| <= sqrt(2)-1 over sqrt(2)+1, before the polynomial.
  y = (m - 1.) / (m + 1.)
  log_m = log_remez_atanh_horner_17(y)

  a = (real(K) * ln2_hi + log_m) + real(K) * ln2_lo
end procedure log_repro


!> Polynomial estimate of 2 * atanh(x) over the log_repro() reduced range.
!!
!! The log_repro() range reduction keeps |x| <= 0.1716.  This first-pass
!! polynomial uses the odd atanh series; a true Remez fit can replace these
!! coefficients without changing the range-reduction structure.
pure function log_remez_atanh_horner_17(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-0.172, 0.172] [nondim]
  real :: a
    !< Approximation of 2 * atanh(x) [nondim]
  real :: x2
    !< x squared [nondim]

  x2 = x * x
  a = 2. * x * (1. + x2 * (0.333333333333333333333333333333333333 + &
      x2 * (0.2 + x2 * (0.142857142857142857142857142857142857 + &
      x2 * (0.111111111111111111111111111111111111 + &
      x2 * (0.0909090909090909090909090909090909091 + &
      x2 * (0.0769230769230769230769230769230769231 + &
      x2 * (0.0666666666666666666666666666666666667 + &
      x2 * (0.0588235294117647058823529411764705882)))))))))
end function log_remez_atanh_horner_17

end submodule MOM_log
