! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This submodule provides a bitwise reproducible implementation of log().

submodule (MOM_intrinsic_functions) MOM_log

use MOM_log_data_n128, only : log_ndiv, log_table_start
use MOM_log_data_n128, only : log_invc_lookup, logc_lookup
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
  integer(kind=int_kind), parameter :: log_table_step_bits = exp_stride / log_ndiv
    !< Spacing in integer representation between adjacent log table entries
  integer(kind=int_kind), parameter :: log_offset_bits = transfer(log_table_start, int_mold)
    !< Bit pattern for the lower end of the log table range, 0.6875 (=11/16)

  integer(kind=int_kind) :: xb, mb
    ! Bit representations of x and its normalized significand
  integer(kind=int_kind) :: tmp, table_bin
    ! Offset from log table base and table-scale bin number
  integer(kind=int_kind) :: raw_exp
    ! Biased IEEE exponent field
  integer(kind=int_kind) :: frac_bits
    ! Fraction bits of a subnormal input
  integer(kind=int_kind) :: K
    ! Binary exponent in x = 2**K m [nondim]
  integer :: top_frac_bit
    ! Highest set fraction bit of a subnormal input
  integer :: idiv
    ! Lookup table subdivision index
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
  if (iand(xb, not(sign_mask)) == 0_int_kind) then
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
  ! dedicated log1p path is accurate over approximately [-1/16, 0.065], and
  ! x - 1 is exact in this range by Sterbenz's lemma.
  u = x - 1.
  if ((u >= -1. / 16.) .and. (u <= 0.064697265625)) then
    a = log1p_compensated_near(u)
    return
  endif

  ! Range reduction: decompose x = 2**K z with z in the lookup table interval.
  ! Subnormals are first multiplied by an exact power of two, then compensated
  ! in K after the table decomposition.
  scaled_subnormal = raw_exp == 0_int_kind
  if (scaled_subnormal) then
    ! Normalize subnormal inputs using integer bit operations rather than
    ! arithmetic scaling, which can be defeated by flush-to-zero modes.
    frac_bits = iand(xb, exp_stride - 1_int_kind)
    do top_frac_bit = expbit - 1, 0, -1
      if (btest(frac_bits, top_frac_bit)) exit
    enddo
    xb = ishft(int(expbias - expbit + top_frac_bit, int_kind), expbit) + &
        ishft(frac_bits - ishft(1_int_kind, top_frac_bit), expbit - top_frac_bit)
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


!> Compensated near-one estimate of log1p(x).
pure function log1p_compensated_near(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-1/16, 0.065] [nondim]
  real :: a
    !< Approximation of log1p(x) [nondim]

  real, parameter :: split_scale = 134217728.
    !< Exact power-of-two scale used to split x into high and low parts [nondim]
  real, parameter :: c(3:13) = [ &
       0.3333333333333333495663409158244431039492812795260412782820753513055676037, &
      -0.2500000000000048410898553018121243145959173159979113851695191401132202101, &
       0.1999999999999401954600182759053457391449548132751979049871662291688667322, &
      -0.1666666666544826424331902956193816320602454785603791935540999385211449934, &
       0.142857142926292331379290143400260129042455919821774060063839433541603189, &
      -0.1250000111081452467282831140237236744260745509124580241506815300101707851, &
       0.1111110780004715236340979900688580293441238517810620127033976097471253938, &
      -9.999529690174345780528722225743902151690897998496979148198959705548315274e-2, &
       9.09152833243368478159650834402877235721130565569984927237483052376511638e-2, &
      -8.42729367861804930413103634851785932895197979844461832654488007987301762e-2, &
       7.684686137510498609385153891204008164960977636361018003750735920585730455e-2]
    !< Remez coefficients for log1p(x) - x + x*x/2 on [-1/16, 265/4096] [nondim]

  real :: x2, x3
    !< Powers of the reduced argument [nondim]
  real :: x_hi, x_lo
    !< High and low parts of x [nondim]
  real :: w, hi, lo, tail
    !< Compensated partial sums and polynomial tail [nondim]
  real :: p1, p2, p3, p4
    !< Grouped polynomial partial sums [nondim]

  x2 = x * x
  x3 = x * x2

  p1 = (c(3) + (x * c(4))) + (x2 * c(5))
  p2 = (c(6) + (x * c(7))) + (x2 * c(8))
  p3 = (c(9) + (x * c(10))) + (x2 * c(11))
  p4 = c(12) + (x * c(13))
  tail = x3 * (p1 + (x3 * (p2 + (x3 * (p3 + (x3 * p4))))))

  ! Split x and compensate the dominant x - x*x/2 terms.  This follows the
  ! structure used by high-quality libm log1p kernels and reduces final-rounding
  ! error near x = 0, where the leading terms dominate the result.
  w = x * split_scale
  x_hi = (x + w) - w
  x_lo = x - x_hi
  w = -0.5 * (x_hi * x_hi)
  hi = x + w
  lo = (x - hi) + w
  lo = lo + ((-0.5 * x_lo) * (x_hi + x))

  a = (tail + lo) + hi
end function log1p_compensated_near


end submodule MOM_log
