! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This submodule provides a bitwise reproducible implementation of log().

submodule (MOM_intrinsic_functions) MOM_log

use MOM_log_data_n128, only : log_ndiv, log_table_start
use MOM_log_data_n128, only : log_invc_lookup, logc_lookup
use MOM_log_data_n128, only : log_chi_lookup, log_clo_lookup

implicit none

contains

!> Reproducible natural logarithm function
!!
!! Compute log(x) with bitwise reproducibility across platforms.
module procedure log_repro
  integer(kind=int_rsize), parameter :: near_lo = transfer(0.9375, int_mold)
  integer(kind=int_rsize), parameter :: near_hi = transfer(1.064697265625, int_mold)
  integer(kind=int_rsize) :: xb
    !< Bit representation of x

  xb = transfer(x, int_mold)
  if ((xb < expstep) .or. (xb >= ishft(expmask, expbit)) .or. &
      ((xb >= near_lo) .and. (xb <= near_hi))) then
    a = log_fallback(x)
  else
    a = log_table(xb, 0_int_rsize)
  endif
end procedure log_repro


!> Scalar fallback for exceptional, near-one, and subnormal inputs.
pure function log_fallback(x) result(a)
  !DIR$ ATTRIBUTES NOINLINE :: log_fallback
  real, intent(in) :: x
  real :: a
  integer(kind=int_rsize), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used to normalize subnormal inputs
  integer(kind=int_rsize) :: xb
    ! Bit representation of x
  integer(kind=int_rsize) :: exponent_adjustment
    ! Correction to the normalized exponent for subnormal inputs
  integer(kind=int_rsize) :: raw_exp
    ! Biased IEEE exponent field
  integer(kind=int_rsize) :: frac_bits
    ! Fraction bits of a subnormal input
  integer :: top_frac_bit
    ! Highest set fraction bit of a subnormal input
  real :: u
    ! Difference from 1 for near-one log1p path [nondim]

  xb = transfer(x, int_mold)

  ! Handle exceptional values before arithmetic range reduction.  This gives
  ! log(0) = -Inf with divide-by-zero, log(negative) = NaN with invalid,
  ! log(+Inf) = +Inf, and NaNs pass through with the usual signaling behavior.
  if (iand(xb, not(signmask)) == 0_int_rsize) then
    a = -1. / abs(x)
    return
  endif

  if (iand(xb, signmask) /= 0_int_rsize) then
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

  ! Normalize subnormals before entering the shared table solver.
  exponent_adjustment = 0_int_rsize
  if (raw_exp == 0_int_rsize) then
    ! Normalize subnormal inputs using integer bit operations rather than
    ! arithmetic scaling, which can be defeated by flush-to-zero modes.
    frac_bits = iand(xb, expstep - 1_int_rsize)
    do top_frac_bit = expbit - 1, 0, -1
      if (btest(frac_bits, top_frac_bit)) exit
    enddo
    xb = ishft(int(expbias - expbit + top_frac_bit, int_rsize), expbit) + &
        ishft(frac_bits - ishft(1_int_rsize, top_frac_bit), expbit - top_frac_bit)
    exponent_adjustment = -Kbias
  endif

  a = log_table(xb, exponent_adjustment)
end function log_fallback


!> Shared table solver for positive normal or normalized subnormal inputs.
pure function log_table(xb, exponent_adjustment) result(a)
  integer(kind=int_rsize), intent(in) :: xb
    !< Bit representation of the positive normal input
  integer(kind=int_rsize), intent(in) :: exponent_adjustment
    !< Correction to the exponent after subnormal normalization; zero for normal inputs
  real :: a
    !< Natural logarithm of the original input [nondim]
  real, parameter :: ln2_hi = 0.69314718036912381649017333984375
    !< Upper 32 bits of ln2 [nondim]
  real, parameter :: ln2_lo = 1.90821492927058770002e-10
    !< Remaining bits of ln2 [nondim]
  integer, parameter :: log_table_shift = trailz(expstep / log_ndiv)
    !< Arithmetic shift for floor division by the table spacing
  integer(kind=int_rsize), parameter :: log_index_mask = int(log_ndiv - 1, int_rsize)
    !< Mask for the table subdivision index
  integer(kind=int_rsize), parameter :: log_offset_bits = transfer(log_table_start, int_mold)
    !< Bit pattern for the lower end of the table range
  integer(kind=int_rsize) :: tmp, table_bin, K, mb
    !< Offset bits, table bin, binary exponent, and reduced significand bits
  integer :: idiv
    !< Lookup table subdivision index
  real :: z, r
    !< Reduced significand and reduced log1p argument [nondim]
  real :: w, hi, lo
    !< Compensated partial sums for log(x) [nondim]

  ! Table reduction uses the ordered bit pattern of positive normal floats,
  ! xb = biased_exp * 2**52 + frac, so one exponent step is an integer stride
  ! of 2**52 and the table bins subdivide that stride.

  ! Use the lower bound of the table as an offset to determine the bin
  tmp = xb - log_offset_bits

  ! Arithmetic shifts round down even for negative offsets.  With a power-of-two
  ! number of table subdivisions, masking gives the nonnegative modulo index.
  table_bin = shifta(tmp, log_table_shift)
  idiv = int(iand(table_bin, log_index_mask))
  K = shifta(tmp, expbit)

  mb = xb - K * expstep
  z = transfer(mb, real_mold)
  K = K + exponent_adjustment

  ! log(x) = K*ln2 + log(c) + log(1+r), where r = z/c - 1.

  ! Compute r as ((z - c_hi) - c_lo)) * (1/c) to avoid precision loss near c.
  r = ((z - log_chi_lookup(idiv)) - log_clo_lookup(idiv)) * log_invc_lookup(idiv)

  w = real(K) * ln2_hi + logc_lookup(idiv)

  hi = w + r
  ! XXX: Prevent (w - (w + r)) => r optimization by adding a subroundoff term.
  hi = hi + epsilon(real_mold) / 256.
  lo = ((w - hi) + r) + real(K) * ln2_lo

  a = hi + (lo + log1p_remez_tail_4(r))
end function log_table


!> Taylor estimate of log1p(x) - x over the log_repro() table-reduced range.
pure function log1p_taylor_tail_6(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-0.004, 0.004] [nondim]
  real :: a
    !< Approximation of log1p(x) - x [nondim]

  real, parameter :: c(2:6) = [ &
      -0.5, &
       0.33333333333333333333333333333333333333333333333333, &
      -0.25, &
       0.20000000000000000000000000000000000000000000000000, &
      -0.16666666666666666666666666666666666666666666666667]
    !< Taylor coefficients for log1p(x) - x through degree 6 [nondim]

  a = x * x * (c(2) + x * (c(3) + x * (c(4) + x * (c(5) + x * c(6)))))
end function log1p_taylor_tail_6


!> Remez estimate of log1p(x) - x over the log_repro() table-reduced range.
!! Uses a degree-4 polynomial for (log1p(x) - x) / x**2, giving a total
!! degree-6 approximation to log1p(x) - x.
pure function log1p_remez_tail_4(x) result(a)
  real, intent(in) :: x
    !< Input value; expected range is approximately [-0.004, 0.004] [nondim]
  real :: a
    !< Approximation of log1p(x) - x [nondim]

  real, parameter :: c(0:4) = [ &
      -0.50000000000000002711451064650790160934116372268686547205906, &
       0.3333333333230997752147523284178620486976229118544577885604334, &
      -0.2499999999767185361215234424897757570704485535532397300358652, &
       0.2000027036539271991906992254523276458937349319710495606164597, &
      -0.1666699786550779057543177943870606385873216854361152274042099]
    !< Sollya Remez coefficients for (log1p(x) - x) / x**2 on [-0.0038910505836871057, +0.0038910505836871057] [nondim]

  real :: x2, x4
    !< Powers of x [nondim]
  real :: p01, p23
    !< Grouped polynomial partial sums [nondim]

  x2 = x * x
  x4 = x2 * x2

  p01 = c(0) + x * c(1)
  p23 = c(2) + x * c(3)

  a = x2 * ((p01 + x2 * p23) + x4 * c(4))
end function log1p_remez_tail_4


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
  integer(kind=int_rsize), parameter :: split_mask = &
      -(2_int_rsize**((digits(real_mold) + 1) / 2))
    !< Retain at most half the significant bits, making x_hi*x_hi exact.

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

  x_hi = transfer(iand(transfer(x, int_mold), split_mask), real_mold)
  x_lo = x - x_hi

  w = -0.5 * (x_hi * x_hi)
  hi = x + w
  ! Clear the lowest bit to prevent (x - (x + w)) => -w optimization
  hi = transfer(iand(transfer(hi, int_mold), not(1_int_rsize), real_mold)
  lo = ((x - hi) + w) + ((-0.5 * x_lo) * (x_hi + x))

  a = (tail + lo) + hi
end function log1p_compensated_near


end submodule MOM_log
