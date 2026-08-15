! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!! This submodule provides a bitwise reproducible implementation of exp().

#include "MOM_exp.h"

submodule (MOM_intrinsic_functions) MOM_exp

use, intrinsic :: iso_fortran_env, only : int32, int64
use, intrinsic :: ieee_arithmetic, only : ieee_rint

! Temporarily used for table generation
use, intrinsic :: iso_fortran_env, only : real128

implicit none

! Molds for transfers and numerical format queries
real, parameter :: real_mold = 0.
  !< Real mold value

integer, parameter :: int_kind &
    = merge(int64, int32, storage_size(real_mold) > storage_size(0_int32))
  !< Integer kind with the same storage size as default real
integer(kind=int_kind), parameter :: int_mold = 0
  !< Integer mold value

! Floating point layout
integer, parameter :: expbit = digits(real_mold) - 1
  !< Position of lowest exponent bit
integer, parameter :: signbit = storage_size(real_mold) - 1
  !< Position of sign bit
integer, parameter :: expwidth = signbit - expbit
  !< Number of exponent bits
integer, parameter :: expbias = maxexponent(real_mold) - 1
  !< Exponent bias

! IEEE 754 special values
integer(kind=int_kind), parameter :: pos_inf_bits &
    = ishft(2_int_kind**expwidth - 1_int_kind, expbit)
  !< IEEE +Inf bit pattern
integer(kind=int_kind), parameter :: neg_inf_bits &
    = ior(pos_inf_bits, ishft(-1_int_kind, signbit))
  !< IEEE -Inf bit pattern

! Fast integer rounding offset
real, parameter :: round_bias = 1.5 * 2_int_kind**(digits(real_mold) - 1)
  !< Binary offset used to trigger rounding of fractional values

contains

!> Reproducible exponential function
!!
!! Compute exp(x) with bitwise reproducibility across platforms.
module procedure exp_repro
  ! ln2 estimates
  real, parameter :: ln2 = 0.693147180559945309417232121458176568
    !< log(2): 0.693147180559945309417232... [nondim]
  real, parameter :: I_ln2 = 1.44269504088896340735992468100189214
    !< 1 / ln2: 1.4426950408889634073599... [nondim]

  ! The max and min x values between which x remains finite.
  real, parameter :: xmax &
      = real(maxexponent(real_mold) + 1, kind(real_mold)) * ln2
    !< Largest x value before exp(x) overflow [nondim]
  real, parameter :: xmin &
      = real(minexponent(real_mold) - digits(real_mold) - 1, kind(real_mold)) * ln2
    !< Smallest x value before exp(x) underflow [nondim]

  ! Double-double precision of ln2 used in range reduction (Cody-Waite)
  ! NOTE: This split assumes double-precision.
  real, parameter :: ln2_hi = 0.69314718036912381649017333984375
    !< Upper 32 bits of ln2: 6.93147180369123816490e-01 [nondim]
  real, parameter :: ln2_lo = 1.90821492927058770002e-10
    !< Lower precision bits of ln2: 1.90821492927058770002e-10 [nondim]

  ! Table constants (may not be used)
  integer, parameter :: NTABLE = 128
  real, parameter :: I_NTABLE = 1. / real(NTABLE)
  real, parameter :: TABLE_INV_LN2 = NTABLE * I_ln2
  real, parameter :: TABLE_LN2_HI = I_NTABLE * ln2_hi
  real, parameter :: TABLE_LN2_LO = I_NTABLE * ln2_lo
  integer(int_kind), parameter :: index_mask = int(NTABLE - 1, int_kind)
  integer(int_kind), parameter :: abs_mask = int(z'7FFFFFFFFFFFFFFF', int_kind)

  ! More table stuff
  real :: Z
  integer(kind=int_kind) :: iz, kz
  integer(int32) :: table_index

  integer :: i

  ! TODO: These are convenient for table generation but may not be reproducible
  !   and should be replaced in the final version.

  ! Split 2**(i/N) into real precision and a scaled residual:
  !   2**(i/N) = [2**(i/N)] (1 + tail)

  ! 2**(i/N)
  real, parameter :: exp2_table(0:NTABLE-1) &
      = [(2.**(real(i) / real(NTABLE)), i=0,NTABLE-1)]

  ! tail
  real, parameter :: exp2_table_tail(0:NTABLE-1) &
      = [(2._real128**(real(i, real128) / real(NTABLE, real128)) &
          / 2.**(real(i) / real(NTABLE)), i=0,NTABLE-1)] - 1.

  ! Range of K = nint(x / ln2) for which direct exponent scaling is safe.
  ! Beyond this range, a bias is applied to handle subnormals and overflow.
  ! NOTE: Fortran exponent is defined as one less than IEEE exponent.
  real, parameter :: Kmin = real(minexponent(real_mold) + 1)
    !< Minimum K before subnormal scaling is needed
    !! Kmin = (minexponent() - 1) + 1 (for min exp(r)) + 1 (safety buffer)
  real, parameter :: Kmax = real(maxexponent(real_mold) - 2)
    !< Maximum K before overflow scaling is needed
    !! Kmax = (maxexponent() - 1) - 0 (max exp(r)) - 1 (safety buffer)
  integer(kind=int_kind), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used for overflow and subnormal scaling
    !! Any bias which rescales 2**K exp(r) to O(1) works here.

  ! Range-reduction variables
  real :: xc
    ! x after being bounded by xmin and xmax [nondim]
  real :: x_ln2
    ! Cached value of x / ln2 [nondim]
  real :: K
    ! Nearest IEEE-rounded integer to (x / ln2) [nondim]
    ! NOTE: K is stored as real to avoid int/real type conversions
  real :: r
    ! Range-reduced input, r = x - K ln2 [nondim]
  real :: e
    ! Exponent of range-reduced input, e = 2**(table_index/NTABLE) exp(r) [nondim]
  real :: expm1_r
    ! Approximation to exp(r) - 1 [nondim]
  real :: scale
    ! Table value for 2**(table_index/NTABLE) [nondim]
  real :: tail
    ! Relative table correction for scale [nondim]

  integer(kind=int_kind) :: xb, Kb, eb
    ! Bit representations of x, K, e

  integer(kind=int_kind) :: j
    ! Bias added to K to compensate for exponent K beyond {-1022,..,+1023}.
  integer(kind=int_kind) :: fb
    ! Bit representation of the j-bias, 2**j

  logical :: nonfinite
    ! True if input is a nonfinite float (+/-Inf, NaN)

  ! 1. Nonfinite handling
  ! ---------------------
  ! Nonfinites must be handled first to prevent their appearance in
  ! calculations, which may raise unwanted floating point signals.

  xb = transfer(x, int_mold)
  nonfinite = iand(xb, pos_inf_bits) == pos_inf_bits

  if (nonfinite) then
    ! exp(-Inf) = 0, otherwise pass-through +Inf and +/-NaN values
    ! Compute x + x to trigger `Invalid` for signaled NaNs.
    a = merge(0., x + x, xb == neg_inf_bits)
    return
  endif

  ! 2. Range Reduction
  ! ------------------
  ! Apply a range reduction of r = x - K ln2, so that exp(x) = 2**K exp(r).
  ! If K = nint(x / ln2) then r is in [-1/2 ln2, 1/2 ln2] and exp(r) can be
  ! estimated by a highly accurate polynomial.

  ! First limit x to the precision-based numerical range of exp().
  ! This allows for safe detection of over/underflow in the subnormal handler.
  xc = min(max(x, xmin), xmax)
  !xc = x

  ! Compute K = nint(x / ln2)
  !*!x_ln2 = xc * I_ln2
  !*!K = NEAREST_INT(x_ln2)
  !*!!K = anint(x_ln2)

  ! Table method
  x_ln2 = xc * TABLE_INV_LN2

  ! Z = N*K + j
  !   N = table
  !   K = nint(x/ln2), the ln2 integral step
  !   j = the j'th (of N) ln2 integral substep

  !! This is slower?
  !Z = NEAREST_INT(x_ln2)
  !iz = transfer(Z + round_bias, int_mold)

  ! Faster, and does not even require the NEAREST_INT() macro
  Z = x_ln2 + round_bias
  iz = transfer(Z, int_mold)
  Z = Z - round_bias

  table_index = iand(iz, index_mask)
  !table_index = modulo(iz, NTABLE)

  kz = iand(iz, not(index_mask))
  K = (transfer(kz, real_mold) - round_bias) * I_NTABLE

  ! NOTE: Performance may degrade if NEAREST_INT is not inlined.

  ! Since K ~ x/ln2, the terms in r = x - K ln2 will nearly cancel and there is
  ! some expected loss of precision.

  ! To compensate, we use a Cody-Waite correction which separates ln2 into its
  ! upper 32 bits and a lower residual.
  !r = (xc - K * ln2_hi) - K * ln2_lo

  ! Non-cody-waite
  !r = xc - K * ln2

  ! Table Cody-waite
  r = (xc - Z * TABLE_LN2_HI) - Z * TABLE_LN2_LO

  ! NOTE: Aggressive optimizers may reduce this to x - K * (ln2_hi + ln2_lo)
  ! which is no better than x - K ln2.

  ! 3. Polynomial approximation
  ! ---------------------------
  ! Estimate exp(r) = 1 + r * P(r) where P(r) is a 10th order Remez minimax
  ! polynomial of (exp(r) - 1)/r.  This form ensures that exp(0) = 1 exactly.

  !e = exp_remez_horner_10(r)
  !e = exp_remez_estrin_10(r)
  !e = exp_remez_estrin_full_10(r)

  ! Table evaluation
  scale = exp2_table(table_index)
  tail = exp2_table_tail(table_index)
  expm1_r = exp_remez_expm1_estrin_4(r)
  ! Evaluate the small correction before the final addition to scale.
  e = scale + (scale * (tail + expm1_r))

  ! 4. Unscaling
  ! ------------
  ! Compute exp(x) = 2**K e, an exact power-of-2 calculation.
  ! Adjust scaling to compensate for subnormal output.

  ! exp(r) has range [0.707, 1.414], so it shifts K by either 0 or -1.
  ! A resolved exponent is in the range {-1022,1023}, so K must be in
  ! {-1021,1023}.  (For safety, we further reduce the range by 1).

  ! Determine if K is outside the supported exponent range.
  ! If so, then apply a bias j to normalize the exponent.
  ! Kbias is chosen so that the exponent is "something near 1".
  j = merge(Kbias, 0_int_kind, K < Kmin) + merge(-Kbias, 0_int_kind, K > Kmax)

  ! Get the bit representation of exp(r)
  eb = transfer(e, int_mold)

  ! This is a fast alternative to int(K), similar to fast_rint().
  ! Kb includes the biased 2**52 in its exponent field, but these are dropped
  ! during the later ishft() operation.
  Kb = transfer(K + round_bias, int_mold)

  ! Rescale to exp(r) to exp(x), possibly including the j bias.
  eb = eb + ishft(Kb + j, expbit)
  a = transfer(eb, real_mold)

  ! Undo the 2**j bias as floating point multiplication.
  ! - For "normals", this has no effect.
  ! - For subnormals, this will force subnormal estimation (if enabled).
  ! - For resolvable K beyond this range, it triggers an over/underflow.
  ! - Extreme values of K have already been filtered out by the min/max step.
  fb = ishft(int(expbias, int_kind) - j, expbit)
  a = a * transfer(fb, real_mold)
end procedure exp_repro


!> Remez minimax polynomial for exp(x) over [-1/2 ln2, 1/2 ln2] in Horner form
!!
!! Approximate (exp(x) - 1) / x, then compute 1 + x * P(x).
!! Coefficients generated by Sollya fpminimax.
pure function exp_remez_horner_10(x) result(e)
  real, intent(in) :: x
    !< Input value; expected range is [-1/2 ln2, 1/2 ln2]
  real :: e
    !< Approximation of exp(x)

  real :: p
    !< Polynomial approximation to (exp(x) - 1) / x [nondim]

  !> fpminimax coefficients for (exp(x) - 1) / x on [-1/2 ln2, 1/2 ln2]
  real, parameter :: c(0:10) = [ &
      1.0, &
      0.50000000000000055511151231257827021181583404541015625, &
      0.1666666666666660745477201999165117740631103515625, &
      4.166666666657388440331288848028634674847126007080078125e-2, &
      8.333333333377164475752607586400699801743030548095703125e-3, &
      1.38888889322647565531532176663631616975180804729461669921875e-3, &
      1.984126974698501824807828075591942251776345074176788330078125e-4, &
      2.480150459644261430602885099006016389466822147369384765625e-5, &
      2.755738179851631320657146320685093598967796424403786659240723e-6, &
      2.76262647076892519593864332161370356288898619823157787322998e-7, &
      2.506210200218863960750001593138364119894845316594000905752182e-8 &
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


pure function exp_remez_estrin_10(x) result(e)
  real, intent(in) :: x
    !< Input value; expected range is [-1/2 ln2, 1/2 ln2]
  real :: e
    !< Approximation of exp(x)

  e = 1 + exp_remez_expm1_estrin_10(x)
end function exp_remez_estrin_10


pure function exp_remez_expm1_estrin_10(x) result(e)
  real, intent(in) :: x
    !< Input value; expected range is [-1/2 ln2, 1/2 ln2]
  real :: e
    !< Approximation of exp(x) - 1 [nondim]

  real :: po
    !< Odd polynomial: c1 + c3*x² + c5*x⁴ + c7*x⁶ + c9*x⁸ [nondim]
  real :: pe
    !< Even polynomial: c0 + c2*x² + c4*x⁴ + c6*x⁶ + c8*x⁸ + c10*x¹⁰ [nondim]
  real :: x2
    !< x squared [nondim]
  real :: p
    !< Combined polynomial (exp(x) - 1) / x [nondim]

  !> fpminimax coefficients for (exp(x) - 1) / x on [-1/2 ln2, 1/2 ln2]
  real, parameter :: c(0:10) = [ &
      1.0, &
      0.50000000000000055511151231257827021181583404541015625, &
      0.1666666666666660745477201999165117740631103515625, &
      4.166666666657388440331288848028634674847126007080078125e-2, &
      8.333333333377164475752607586400699801743030548095703125e-3, &
      1.38888889322647565531532176663631616975180804729461669921875e-3, &
      1.984126974698501824807828075591942251776345074176788330078125e-4, &
      2.480150459644261430602885099006016389466822147369384765625e-5, &
      2.755738179851631320657146320685093598967796424403786659240723e-6, &
      2.76262647076892519593864332161370356288898619823157787322998e-7, &
      2.506210200218863960750001593138364119894845316594000905752182e-8 &
  ]

  x2 = x * x

  ! Odd coefficients: c1 + c3*x² + c5*x⁴ + c7*x⁶ + c9*x⁸
  po = c(9)
  po = c(7) + x2 * po
  po = c(5) + x2 * po
  po = c(3) + x2 * po
  po = c(1) + x2 * po

  ! Even coefficients: c0 + c2*x² + c4*x⁴ + c6*x⁶ + c8*x⁸ + c10*x¹⁰
  pe = c(10)
  pe = c(8) + x2 * pe
  pe = c(6) + x2 * pe
  pe = c(4) + x2 * pe
  pe = c(2) + x2 * pe
  pe = c(0) + x2 * pe

  ! Combine: p = pe + x * po
  p = pe + x * po

  ! Final assembly: exp(x) - 1 = x * p(x)
  e = x * p
end function exp_remez_expm1_estrin_10


pure function exp_remez_expm1_estrin_4(x) result(e)
  real, intent(in) :: x
    !< Input value; expected range is [-ln2/256, ln2/256]
  real :: e
    !< Approximation of exp(x) - 1 [nondim]

  real :: x2, x4
    !< Powers of x [nondim]
  real :: p01, p23, p
    !< Polynomial partial sums [nondim]

  !> fpminimax coefficients for (exp(x) - 1) / x on [-ln2/256, ln2/256]
  real, parameter :: c(0:4) = [ &
      1., &
      0.4999999999999766853164828717126511037349700927734375, &
      0.166666666666670015839457619222230277955532073974609375, &
      4.1666679392304360740606483659576042555272579193115234375e-2, &
      8.3333340579158539374038383584775147028267383575439453125e-3 &
  ]

  x2 = x * x
  x4 = x2 * x2

  p01 = c(0) + (c(1) * x)
  p23 = c(2) + (c(3) * x)
  p = (p01 + (x2 * p23)) + (x4 * c(4))

  ! Final assembly: exp(x) - 1 = x * p(x)
  e = x * p
end function exp_remez_expm1_estrin_4


pure function exp_remez_estrin_full_10(x) result(e)
  real, intent(in) :: x
    !< Input value; expected range is [-1/2 ln2, 1/2 ln2]
  real :: e
    !< Approximation of exp(x)

  real :: x2, x4, x8
    !< Powers of x [nondim]
  real :: p01, p23, p45, p67, p89
    !< Paired terms [nondim]
  real :: p03, p47, p8_10
    !< Quad terms [nondim]
  real :: p07, p
    !< Combined polynomial [nondim]

  !> fpminimax coefficients for (exp(x) - 1) / x on [-1/2 ln2, 1/2 ln2]
  real, parameter :: c(0:10) = [ &
      1.0, &
      0.50000000000000055511151231257827021181583404541015625, &
      0.1666666666666660745477201999165117740631103515625, &
      4.166666666657388440331288848028634674847126007080078125e-2, &
      8.333333333377164475752607586400699801743030548095703125e-3, &
      1.38888889322647565531532176663631616975180804729461669921875e-3, &
      1.984126974698501824807828075591942251776345074176788330078125e-4, &
      2.480150459644261430602885099006016389466822147369384765625e-5, &
      2.755738179851631320657146320685093598967796424403786659240723e-6, &
      2.76262647076892519593864332161370356288898619823157787322998e-7, &
      2.506210200218863960750001593138364119894845316594000905752182e-8 &
  ]

  ! Powers of x
  x2 = x * x
  x4 = x2 * x2
  x8 = x4 * x4

  ! Level 1: pairs (all independent)
  p01 = c(0) + c(1) * x
  p23 = c(2) + c(3) * x
  p45 = c(4) + c(5) * x
  p67 = c(6) + c(7) * x
  p89 = c(8) + c(9) * x

  ! Level 2: quads (all independent, depend on level 1)
  p03 = p01 + x2 * p23
  p47 = p45 + x2 * p67
  p8_10 = p89 + x2 * c(10)

  ! Level 3: octets
  p07 = p03 + x4 * p47

  ! Level 4: final combine
  p = p07 + x8 * p8_10

  ! Final assembly: exp(x) = 1 + x * p(x)
  e = 1 + x * p
end function exp_remez_estrin_full_10


!> An optimized nearest-integer function for floating point reals.
!!
!! The value x is shifted from 2**K (1 + a) to 2**(digits-1+K) (1.5 + [a]),
!! causing the fractional part to be rounded according to the current IEEE
!! settings.  In almost all cases, this is nearest ties-to-even.
!!
!! The behavior of this function does not match nint() or anint().  The nint()
!! function always ties away from zero, e.g. nint(2.5) = 3.
!!
!! The +0.5 ensures that the biased exponent of negative numbers does not drop
!! by one, which can cause half-value rounding.
!!
!! It is essential that compilers not reduce (x+b)-b to x.  This can typically
!! be ensured as long as parentheses are respected.  This is managed by the
!! ENABLE_FAST_RINT macro in MOM_exp.h and assigned to NEAREST_INT().  If
!! unset, then ieee_rint() is used.
pure function fast_rint(x) result(n)
  real, intent(in) :: x
    !< Real value to be rounded to the nearest integer
  real :: n
    !< Nearest integer to x, stored as a real

  n = (x + round_bias) - round_bias
end function fast_rint

end submodule MOM_exp
