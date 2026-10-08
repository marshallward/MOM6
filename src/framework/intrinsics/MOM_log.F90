! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This submodule provides a bitwise reproducible implementation of log().

submodule (MOM_intrinsic_functions) MOM_log

use MOM_log_data_n128, only : log_ndiv, log_invc_lookup, logc_lookup
use MOM_log_data_n128, only : log_chi_lookup, log_clo_lookup

implicit none

! IEEE 754 masks
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
  integer(kind=int_kind), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used to normalize subnormal inputs
  real, parameter :: scale_up = transfer(ishft(int(expbias, int_kind) + Kbias, expbit), real_mold)
    !< Exact power-of-two scale factor for subnormal inputs [nondim]
  integer(kind=int_kind), parameter :: two_to_expbit = 2_int_kind**expbit
    !< Integer value 2**expbit, the spacing between IEEE exponent fields
  integer(kind=int_kind), parameter :: log_table_step_bits = two_to_expbit / log_ndiv
    !< Spacing in integer representation between adjacent log table entries
  integer(kind=int_kind), parameter :: log_offset_bits = 4604367669032910848_int_kind
    !< Bit pattern for the lower end of the log table range, 0x1.6p-1

  integer(kind=int_kind) :: xb, mb
    ! Bit representations of x and its normalized significand
  integer(kind=int_kind) :: tmp, table_bin
    ! Offset from log table base and table-scale bin number
  integer(kind=int_kind) :: raw_exp
    ! Biased IEEE exponent field
  integer(kind=int_kind) :: K
    ! Binary exponent in x = 2**K m [nondim]
  integer :: idiv
    ! Lookup table subdivision index
  real :: xs
    ! Input value, possibly scaled to normalize subnormal numbers [nondim]
  real :: z
    ! Significand of x, adjusted into the log table interval [nondim]
  real :: r
    ! Reduced argument, r = z / c - 1 [nondim]
  real :: w, hi, lo
    ! Double-real partial sums for log(x) [nondim]
  logical :: scaled_subnormal
    ! True if x was scaled up to normalize a subnormal input

  xb = transfer(x, int_mold)

  ! Handle exceptional values before arithmetic range reduction.  This gives
  ! log(0) = -Inf with divide-by-zero, log(negative) = NaN with invalid,
  ! log(+Inf) = +Inf, and NaNs pass through with the usual signaling behavior.
  if (x == 0.) then
    a = -1. / abs(x)
    return
  endif

  if (x == 1.) then
    a = 0.
    return
  endif

  if (iand(xb, sign_mask) /= 0_int_kind) then
    a = (x - x) / (x - x)
    return
  endif

  raw_exp = iand(ishft(xb, -expbit), expmask)

  if (raw_exp == expmask) then
    a = x + x
    return
  endif

  ! Range reduction: decompose x = 2**K z with z in the lookup table interval.
  ! Subnormals are first multiplied by an exact power of two, then compensated
  ! in K after the table decomposition.
  scaled_subnormal = raw_exp == 0_int_kind
  if (scaled_subnormal) then
    xs = x * scale_up
    xb = transfer(xs, int_mold)
    raw_exp = iand(ishft(xb, -expbit), expmask)
  endif

  tmp = xb - log_offset_bits
  table_bin = floor_div_int(tmp, log_table_step_bits)
  idiv = int(modulo(table_bin, int(log_ndiv, int_kind)))
  K = floor_div_int(tmp, two_to_expbit)

  mb = xb - K * two_to_expbit
  z = transfer(mb, real_mold)
  if (scaled_subnormal) K = K - Kbias

  ! log(x) = K*ln2 + log(c) + log1p(r), where r = z / c - 1.
  ! The table makes |r| small; splitting c into high and low parts keeps the
  ! reduction accurate when z is close to c.
  r = ((z - log_chi_lookup(idiv)) - log_clo_lookup(idiv)) * log_invc_lookup(idiv)
  w = real(K) * ln2_hi + logc_lookup(idiv)
  hi = w + r
  lo = (w - hi + r) + real(K) * ln2_lo

  a = hi + (lo + log1p_taylor_tail_6(r))
end procedure log_repro


!> Floor division of signed integers for a positive denominator.
pure function floor_div_int(n, d) result(q)
  integer(kind=int_kind), intent(in) :: n
    !< Numerator
  integer(kind=int_kind), intent(in) :: d
    !< Positive denominator
  integer(kind=int_kind) :: q
    !< floor(n / d)

  q = (n - modulo(n, d)) / d
end function floor_div_int


!> Taylor estimate of log1p(x) - x over the log_repro() table-reduced range.
pure function log1p_taylor_tail_6(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-0.004, 0.004] [nondim]
  real :: a
    !< Approximation of log1p(x) - x [nondim]

  a = x * x * (-0.5 + x * (1. / 3. + x * (-0.25 + x * (0.2 + x * (-1. / 6.)))))
end function log1p_taylor_tail_6

end submodule MOM_log
