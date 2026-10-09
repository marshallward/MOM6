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
  integer(kind=int_kind), parameter :: log_table_step_bits = exp_stride / log_ndiv
    !< Spacing in integer representation between adjacent log table entries
  integer(kind=int_kind), parameter :: log_offset_bits = 4604367669032910848_int_kind
    !< Bit pattern for the lower end of the log table range, 0.6875 (=11/16)

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
  real :: u
    ! Difference from 1 for near-one log1p path [nondim]
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

  ! Avoid table-reduction cancellation for values very close to 1.  The
  ! existing log1p Remez tail is accurate over approximately [-1/128, 1/128],
  ! and x - 1 is exact in this range by Sterbenz's lemma.
  u = x - 1.
  if (abs(u) <= 1. / 128.) then
    a = u + log1p_remez_tail_8(u)
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

  ! Table reduction uses the ordered bit pattern of positive normal floats,
  ! xb = biased_exp * 2**52 + frac, so one exponent step is an integer stride
  ! of 2**52 and the table bins subdivide that stride.

  ! Use the lower bound of the table as an offset to determine the bin
  tmp = xb - log_offset_bits

  ! TODO: Replace these modulo() calls
  table_bin = floor_div_int(tmp, log_table_step_bits)
  idiv = int(modulo(table_bin, int(log_ndiv, int_kind)))
  K = floor_div_int(tmp, exp_stride)

  mb = xb - K * exp_stride
  z = transfer(mb, real_mold)
  if (scaled_subnormal) K = K - Kbias

  ! log(x) = K*ln2 + log(c) + log(1+r), where r = z/c - 1.

  ! Compute r as ((z - c_hi) - c_lo)) * (1/c) to avoid precision loss near c.
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


!> Remez estimate of log1p(x) - x over the near-one log_repro() range.
pure function log1p_remez_tail_6(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-1/512, 1/512] [nondim]
  real :: a
    !< Approximation of log1p(x) - x [nondim]

  real, parameter :: c(2:6) = [ &
      -0.500000000000000001531992077481654475668128157188753876631778885247311583, &
       0.3333333333322630456631544068975149947601565365240798082150360405372396857, &
      -0.2499999999976636390789139216550810004312450069947905784129872818608102103, &
       0.2000008112826621083191778349879309697511001440446713500798235478459693848, &
      -0.166667650690583434892505189737659407295169538338185395188807100923723888]
    !< Remez coefficients for log1p(x) - x on [-1/512, 1/512] [nondim]

  a = x * x * (c(2) + x * (c(3) + x * (c(4) + x * (c(5) + x * c(6)))))
end function log1p_remez_tail_6


!> Remez estimate of log1p(x) - x over the near-one log_repro() range.
pure function log1p_remez_tail_8(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-1/128, 1/128] [nondim]
  real :: a
    !< Approximation of log1p(x) - x [nondim]

  real, parameter :: c(2:8) = [ &
      -0.50000000000000000146880644447745834276117239589584147686810426077429923, &
       0.333333333333383439081489510558618980088213205599752827834025980978889853, &
      -0.250000000000300387046095139450341182750948342481656948559960050179908581, &
       0.199999996718934020585979325801837669966856969675133175387608887567261307, &
      -0.166666664868771236962982455327063859706607758596032524924398075146814502, &
       0.142871285441583989558603865248163259413977031840605019462415377267409527, &
      -0.125015076930830173046952639984941031775223381177899836187586194388830784]
    !< Remez coefficients for log1p(x) - x on [-1/128, 1/128] [nondim]

  a = x * x * (c(2) + x * (c(3) + x * (c(4) + x * (c(5) + &
      x * (c(6) + x * (c(7) + x * c(8)))))))
end function log1p_remez_tail_8


end submodule MOM_log
