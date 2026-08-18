! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!! This submodule provides a bitwise reproducible implementation of exp().

#include "MOM_exp.h"

submodule (MOM_intrinsic_functions) MOM_exp

use, intrinsic :: iso_fortran_env, only : int32, int64
use, intrinsic :: ieee_arithmetic, only : ieee_rint

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
integer(int_kind), parameter :: round_bias_bits &
    = transfer(round_bias, int_mold)
  !< Bit representation of rounding bias

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

  ! Threshold for special-case handling (nonfinite, overflow, underflow).
  ! exp(x) overflows for x > ~709.78 and underflows for x < ~-745.13.
  ! We use 1024.0 as the threshold (exponent = 10); any |x| >= 1024 will
  ! definitely overflow or underflow, so we handle those specially.
  ! Values in [512, 1024) are handled by the normal path with j-bias.
  integer(int_kind), parameter :: exp_threshold = int(expbias + 10, int_kind)
    !< Exponent of 1024.0 (2^10); values with larger exponent need special handling

  ! Double-real precision of ln2 used in Cody-Waite range reduction
  ! NOTE: This split assumes real64 precision.
  real, parameter :: ln2_hi = 0.69314718036912381649017333984375
    !< Upper 32 bits of ln2: 6.93147180369123816490e-01 [nondim]
  real, parameter :: ln2_lo = 1.90821492927058770002e-10
    !< Lower precision bits of ln2: 1.90821492927058770002e-10 [nondim]

  ! Table constants (may not be used)
  ! Subdivide [-1/2 ln2, 1/2 ln2] to further reduce the approximation range.
  ! This allows for a smaller polynomial at the cost of a lookup table.
  integer, parameter :: ndiv = 128
    !< Number of domain subdivisions
  real, parameter :: I_ndiv = 1. / real(ndiv)
    !< 1 / ndiv [nondim]
  real, parameter :: n_ln2 = ndiv * I_ln2
    !< ndiv / ln2 [nondim]
  real, parameter :: ln2_ndiv_hi = I_ndiv * ln2_hi
    !< Upper 32 bits of ln2 / ndiv [nondim]
  real, parameter :: ln2_ndiv_lo = I_ndiv * ln2_lo
    !< Lower precision bits of ln2 / ndiv [nondim]
  integer(int_kind), parameter :: index_mask = int(ndiv - 1, int_kind)
    ! TODO

  ! The following lookup tables are used to reproducibily compute 2**(j/n),
  ! where 2**(j/n) is split into real precision and a scaled residual.
  !   2**(j/ndiv) = real(2**(j/N)) * (1 + residual)

  ! real, parameter :: idiv_scale_lookup(0:ndiv-1) &
  !     = [(2.**(real(i) / real(ndiv)), i=0,ndiv-1)]
  real, parameter :: idiv_scale_lookup(0:ndiv-1) = [ &
    1.00000000000000000e+00, 1.00542990111280273e+00, 1.01088928605170048e+00, 1.01637831491095310e+00, &
    1.02189714865411663e+00, 1.02744594911876375e+00, 1.03302487902122841e+00, 1.03863410196137873e+00, &
    1.04427378242741375e+00, 1.04994408580068721e+00, 1.05564517836055716e+00, 1.06137722728926209e+00, &
    1.06714040067682370e+00, 1.07293486752597556e+00, 1.07876079775711986e+00, 1.08461836221330921e+00, &
    1.09050773266525769e+00, 1.09642908181637688e+00, 1.10238258330784089e+00, 1.10836841172367873e+00, &
    1.11438674259589243e+00, 1.12043775240960675e+00, 1.12652161860824185e+00, 1.13263851959871920e+00, &
    1.13878863475669156e+00, 1.14497214443180417e+00, 1.15118922995298267e+00, 1.15744007363375112e+00, &
    1.16372485877757748e+00, 1.17004376968325019e+00, 1.17639699165028122e+00, 1.18278471098434101e+00, &
    1.18920711500272103e+00, 1.19566439203982733e+00, 1.20215673145270308e+00, 1.20868432362658162e+00, &
    1.21524735998046896e+00, 1.22184603297275762e+00, 1.22848053610687002e+00, 1.23515106393693341e+00, &
    1.24185781207348400e+00, 1.24860097718920482e+00, 1.25538075702469110e+00, 1.26219735039425074e+00, &
    1.26905095719173322e+00, 1.27594177839639200e+00, 1.28287001607877826e+00, 1.28983587340666572e+00, &
    1.29683955465100964e+00, 1.30388126519193581e+00, 1.31096121152476441e+00, 1.31807960126606405e+00, &
    1.32523664315974132e+00, 1.33243254708316150e+00, 1.33966752405330292e+00, 1.34694178623294580e+00, &
    1.35425554693689265e+00, 1.36160902063822475e+00, 1.36900242297459052e+00, 1.37643597075453017e+00, &
    1.38390988196383202e+00, 1.39142437577192624e+00, 1.39897967253831124e+00, 1.40657599381901544e+00, &
    1.41421356237309515e+00, 1.42189260216916558e+00, 1.42961333839197002e+00, 1.43737599744898237e+00, &
    1.44518080697704665e+00, 1.45302799584905262e+00, 1.46091779418064704e+00, 1.46885043333698184e+00, &
    1.47682614593949935e+00, 1.48484516587275239e+00, 1.49290772829126484e+00, 1.50101406962642558e+00, &
    1.50916442759342284e+00, 1.51735904119821474e+00, 1.52559815074453842e+00, 1.53388199784095591e+00, &
    1.54221082540794074e+00, 1.55058487768499997e+00, 1.55900440023783693e+00, 1.56746963996555300e+00, &
    1.57598084510788650e+00, 1.58453826525249375e+00, 1.59314215134226700e+00, 1.60179275568269341e+00, &
    1.61049033194925428e+00, 1.61923513519486373e+00, 1.62802742185734783e+00, 1.63686744976696441e+00, &
    1.64575547815396495e+00, 1.65469176765619430e+00, 1.66367658032673638e+00, 1.67271017964159663e+00, &
    1.68179283050742900e+00, 1.69092479926930528e+00, 1.70010635371852348e+00, 1.70933776310046293e+00, &
    1.71861929812247793e+00, 1.72795123096183767e+00, 1.73733383527370622e+00, 1.74676738619916905e+00, &
    1.75625216037329945e+00, 1.76578843593327273e+00, 1.77537649252652119e+00, 1.78501661131893496e+00, &
    1.79470907500310717e+00, 1.80445416780662393e+00, 1.81425217550039886e+00, 1.82410338540705341e+00, &
    1.83400808640934243e+00, 1.84396656895862598e+00, 1.85397912508338547e+00, 1.86404604839778898e+00, &
    1.87416763411029996e+00, 1.88434417903233453e+00, 1.89457598158696561e+00, 1.90486334181767414e+00, &
    1.91520656139714740e+00, 1.92560594363612503e+00, 1.93606179349229435e+00, 1.94657441757923322e+00, &
    1.95714412417540018e+00, 1.96777122323317588e+00, 1.97845602638795093e+00, 1.98919884696726634e+00  &
  ]   !< Lookup table of 2**(j/n)

  ! real, parameter :: idiv_tail_lookup(0:ndiv-1) &
  !     = [(2._real128**(real(i, real128) / real(ndiv, real128)) &
  !         / 2.**(real(i) / real(ndiv)), i=0,ndiv-1)] - 1.
  real, parameter :: idiv_tail_lookup(0:ndiv-1) = [ &
    0.00000000000000000e+00,  9.44788545172706630e-17, -1.50706697692603887e-17, -5.67915508282501219e-17, &
    4.99974487227263259e-17, -4.82368359999489520e-17,  7.35784687124741823e-18,  5.77323022374195081e-17, &
    8.18931763819551480e-17,  5.32689113998087777e-17,  1.66658814423267469e-18, -1.12811324546182800e-17, &
   -7.40282530942617744e-17, -3.57865976730956277e-18, -6.17065474560869496e-17,  2.91913999994927872e-17, &
   -2.79391148595157333e-17, -5.39928535518428500e-17,  4.77695942525622331e-17, -7.92770143233847309e-17, &
    9.34171060990504558e-17, -5.53452067570747201e-17,  4.58567032666235109e-17,  2.85824304111161432e-17, &
    7.82657325863607615e-17,  4.05362690676921607e-17,  2.82378442595106130e-17, -7.88280226248799127e-17, &
    3.29047266460084165e-17, -1.57920947033478823e-18,  4.72136812117012811e-17,  1.30452770969196587e-17, &
    3.34846233362515241e-17,  3.86111997749256638e-17,  5.52755004850524930e-17, -3.92718417244523391e-17, &
   -6.34655210672948264e-17, -8.68441761486594369e-17, -1.54563428193977334e-17, -8.70763476495455021e-17, &
    3.75085420130312718e-17, -6.61685450352648820e-17, -5.34609900919875074e-18, -2.44372632101501814e-17, &
    2.10230496752157140e-18,  7.77106793750106477e-17,  1.33575100888345415e-17,  6.93829169695920377e-17, &
    1.95725852931120358e-17,  6.63225696167580396e-17, -5.47806912392677801e-17, -4.14083931039262377e-17, &
   -2.15714772512087524e-17, -3.82877665521205353e-17,  6.66380458923219517e-17,  2.39361874002852819e-17, &
    5.68648095791173997e-17,  1.12645233545216843e-18,  7.00787504690699447e-17, -5.01192142783812542e-17, &
   -4.89230675135227563e-17, -3.52600899532694341e-17, -6.87230372090201806e-17,  5.00144666413353212e-18, &
   -6.83580865766192197e-17, -1.13073440929102116e-17, -8.41601163471715593e-18, -2.92479770354365663e-17, &
   -2.09230438184335296e-17, -3.97786458754271236e-17, -3.83346496865429512e-17,  5.76361116448089415e-17, &
   -2.35910947708500527e-17,  7.26007466109857455e-17,  9.50689710108795955e-18, -4.27295613383990616e-17, &
   -6.73521923237468289e-17, -2.83960444104309359e-17, -7.22663547210125660e-17,  5.78612100339591767e-17, &
    5.15483011707867833e-17, -9.41625756887815240e-18,  2.42539857666898063e-17, -6.60431405170770719e-17, &
   -6.43213177542418901e-18, -1.22040076018639209e-17, -6.33616186340129307e-17, -3.78008792463738146e-17, &
    1.53414100536037229e-17,  1.29328555804274517e-17, -4.12336733066114881e-17,  4.70308397446345575e-17, &
   -6.15260289155026467e-17,  5.82784932619527933e-17,  3.54094826264618290e-17, -3.27415713209387642e-17, &
    4.87516052622706170e-17, -5.71856979007783759e-17, -4.71953966459097162e-18, -5.77345195880570602e-17, &
   -1.07724870789340559e-17, -6.22180861853391081e-17,  1.82140544036225897e-17, -6.15553654622763901e-17, &
    1.68548729062897313e-17,  5.35812491776948158e-17,  3.62161593533689420e-17,  8.58837952757414421e-18, &
    1.01562190116415007e-17, -2.86913488918724381e-17, -5.49511896612200476e-17, -5.56965572431627048e-17, &
    1.79012690760451314e-17, -3.22117663462001635e-17,  5.26537076855627407e-17,  3.50898664024030334e-17, &
   -3.26692410090131783e-17, -4.36575930080793745e-17,  1.79639326598330223e-17,  3.43009252752141664e-17, &
   -5.54506561863942674e-17, -5.14900974545773279e-17,  5.33680587851415070e-17,  3.49897866119297325e-17, &
    4.57849152770600949e-17, -5.24193457539389921e-17,  2.04142788975783032e-17,  4.12484284860648776e-18  &
  ]   !< Lookup table of r, where 2^(j/N) = real64(2^(j/N)) * (1 + r)

  ! Range of K = nint(x / ln2) for which direct exponent scaling is safe.
  ! Beyond this range, a bias is applied to handle subnormals and overflow.
  ! NOTE: Fortran exponent is defined as one less than IEEE exponent.
  integer, parameter :: Kmin = minexponent(real_mold) + 1
    !< Minimum K before subnormal scaling is needed
    !! Kmin = (minexponent() - 1) + 1 (min exp(r)) (+1 safety buffer)
  integer, parameter :: Kmax = maxexponent(real_mold) - 2
    !< Maximum K before overflow scaling is needed
    !! Kmax = (maxexponent() - 1) - 0 (max exp(r)) (-1 safety buffer)
  integer(kind=int_kind), parameter :: Kbias = maxexponent(real_mold) - 2
    !< Exponent adjustment used for overflow and subnormal scaling
    !! Any bias which rescales 2**K exp(r) to O(1) works here.

  ! More table stuff

  ! Range-reduction variables
  integer :: idiv
    ! Subdivision index
  integer(kind=int_kind) :: K
    ! Nearest IEEE-rounded integer to x/ln2 [nondim]
  real :: Z
    ! Nearest integer to ndiv*K + idiv for ndiv subdivisions
    ! NOTE: Z is stored as real to avoid int-real type conversions
  integer(kind=int_kind) :: Zi
    ! Integer representation of Z
  real :: r
    ! Range-reduced input, r = x - Z ln2/N [nondim]
  real :: e
    ! Exponent of range-reduced input, e = 2**(idiv/ndiv) exp(r) [nondim]
  real :: expm1_r
    ! Approximation to exp(r) - 1 [nondim]
  real :: idiv_scale
    ! Table value for 2**(idiv/ndiv) [nondim]
  real :: idiv_tail
    ! Relative table correction for idiv_scale [nondim]
  integer(kind=int_kind) :: xb, eb
    ! Bit representations of x, e
  integer(kind=int_kind) :: x_exp
    ! Exponent bits of x

  integer(kind=int_kind) :: j
    ! Bias added to K to compensate for exponent K beyond {-1022,..,+1023}.
  integer(kind=int_kind) :: fb
    ! Bit representation 2**j, the K exponent rescale

  ! 1. Special case handling
  ! ------------------------
  ! Handle nonfinite inputs and extreme values that would overflow/underflow.
  ! Extract exponent bits and check if |x| is too large or nonfinite.

  xb = transfer(x, int_mold)
  x_exp = iand(ishft(xb, -expbit), int(2**expwidth - 1, int_kind))

  if (x_exp >= exp_threshold) then
    ! Either nonfinite (exp = 2047) or |x| >= 512
    if (x_exp >= int(2**expwidth - 1, int_kind)) then
      ! Nonfinite: +/-Inf or NaN
      ! exp(-Inf) = 0, otherwise pass-through +Inf and +/-NaN values
      ! Compute x + x to trigger `Invalid` for signaled NaNs.
      a = merge(0., x + x, xb == neg_inf_bits)
    elseif (xb >= 0) then
      ! Large positive x -> overflow to +Inf
      a = x * huge(x)
    else
      ! Large negative x -> underflow to 0
      a = tiny(x) * (tiny(x) + 0. * x)
    endif
    return
  endif

  ! 2. Range Reduction
  ! ------------------
  ! Apply a range reduction of r = x - K ln2 - (j / N) ln2, so that
  !     exp(x) = 2**K 2**(j/N) exp(r).
  ! If K = nint(x / ln2) then r is in [-1/2N ln2, 1/2N ln2] and exp(r) can be
  ! estimated by a sufficiently accurate polynomial.

  ! Compute Z = N K + j, where r = x - Z (ln 2 / N)
  Z = NEAREST_INT(x * n_ln2)
  Zi = transfer(Z + round_bias, int_mold) - round_bias_bits

  idiv = iand(Zi, index_mask)
  K = (Zi - int(idiv, int_kind)) / ndiv

  ! Since K ~ x/ln2, the terms in r = x - K ln2 will nearly cancel and there is
  ! some expected loss of precision.  To compensate, we use a Cody-Waite
  ! correction which separates ln2 into its upper 32 bits and a lower residual.

  r = (x - Z * ln2_ndiv_hi) - Z * ln2_ndiv_lo

  ! NOTE: Aggressive optimizers may reduce this to x - K * (ln2_hi + ln2_lo)
  ! which is no better than x - K ln2.

  ! 3. Polynomial approximation
  ! ---------------------------
  ! Estimate exp(r) = 1 + r * P(r) where P(r) is a 10th order Remez minimax
  ! polynomial of (exp(r) - 1)/r.  This form ensures that exp(0) = 1 exactly.

  ! Table evaluation
  idiv_scale = idiv_scale_lookup(idiv)
  idiv_tail = idiv_tail_lookup(idiv)
  expm1_r = exp_remez_expm1_estrin_4(r)
  ! Evaluate the small correction before the final addition to idiv_scale
  e = idiv_scale + idiv_scale * (idiv_tail + expm1_r)

  ! 4. Unscaling
  ! ------------
  ! Compute exp(x) = 2**K e, an exact power-of-2 calculation.
  ! Adjust scaling to compensate for subnormal output.

  ! exp(r) has range [0.707/ndiv, 1.414/ndiv], so K shifts by either 0 or -1.
  ! A resolved exponent is in the range {-1022..1023}, so K must be in
  ! {-1021,1023}.  (For safety, we reduce the range to {-1020,1022}.)

  ! Determine if K is outside the supported exponent range.
  ! If so, then apply a bias j to normalize the exponent.
  ! Kbias is chosen so that the exponent is "something near 1".
  j = merge(Kbias, 0_int_kind, K < Kmin) + merge(-Kbias, 0_int_kind, K > Kmax)

  ! Get the bit representation of exp(r)
  eb = transfer(e, int_mold)

  ! Rescale to exp(r) to exp(x), possibly including the j bias.
  eb = eb + ishft(k + j, expbit)
  a = transfer(eb, real_mold)

  ! Undo the 2**j bias as floating point multiplication.
  ! - For "normals", this has no effect.
  ! - For subnormals, this will force subnormal estimation (if enabled).
  ! - For resolvable K beyond this range, it triggers an over/underflow.
  ! - Extreme values of K have already been filtered out by the min/max step.
  fb = ishft(int(expbias, int_kind) - j, expbit)
  a = a * transfer(fb, real_mold)
end procedure exp_repro


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
  real, parameter :: c(0:3) = [ &
      0.4999999999999766853164828717126511037349700927734375, &
      0.166666666666670015839457619222230277955532073974609375, &
      4.1666679392304360740606483659576042555272579193115234375e-2, &
      8.3333340579158539374038383584775147028267383575439453125e-3 &
  ]

  x2 = x * x
  x4 = x2 * x2

  p01 = c(0) + c(1) * x
  p23 = c(2) + c(3) * x

  ! Final assembly: exp(x) - 1 = x + x2*(C2+x*C3) + x4*(C4+x*C5)
  e = (x + x2 * p01) + x4 * p23

end function exp_remez_expm1_estrin_4


!> An optimized nearest-integer function for floating point reals.
!!
!! The value x is shifted from 2**K (1 + a) to 2**(digits-1+K) (1.5 + [a]),
!! causing the fractional part to be rounded according to the current IEEE
!! settings.  In almost all cases, this is nearest ties-to-even.
!!
!! The +0.5 ensures that the biased exponent of negative numbers does not drop
!! by one, which can cause half-value rounding.
!!
!! The behavior of this function does not match nint() or anint().  The nint()
!! function always ties away from zero, e.g. nint(2.5) = 3.
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
