! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This submodule provides a bitwise reproducible implementation of nth_root().

submodule (MOM_intrinsic_functions) MOM_nth_root

use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_fma

implicit none

contains

!> Reproducible nth root function for positive arguments.
!!
!! The input is rescaled by an integer power of 2 so that the iterative solve is
!! done on an argument in the range [2**(-n), 1).  This preserves exact
!! dimensional rescaling by powers of 2**n, following the cuberoot() pattern.
module procedure nth_root
  integer, parameter :: max_nth_root = 32
    !< Largest root degree supported by nth_root [nondim]
  integer, parameter :: halley_iterations_by_n(max_nth_root) = [ &
      0, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, &
      5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8]
    !< Number of Halley iterations used by each supported root degree [nondim]

  real :: xr
    !< The rescaled value of x in the range [2**(-n), 1) [B^n]
  real :: root_xr
    !< The nth root of xr [B]
  real :: root_n
    !< root_xr raised to the nth power [B^n]
  real :: root_nm1
    !< root_xr raised to the n-1 power [B^(n-1)]
  real :: num
    !< Numerator in the Halley iteration [B^n]
  real :: den
    !< Denominator in the Halley iteration [B^n]
  real :: rn_minus
    !< n-1 represented as a real number [nondim]
  real :: rn_plus
    !< n+1 represented as a real number [nondim]
  integer(kind=int64) :: e_x
    !< Integral component of the nth-root exponent of x [nondim]
  integer :: itt
    !< Iteration counter [nondim]

  if ((n < 1) .or. (n > max_nth_root)) then
    root = ieee_value(root, ieee_quiet_nan)
  elseif (x < 0.0) then
    root = ieee_value(root, ieee_quiet_nan)
  elseif (is_nonfinite(x)) then
    ! Pass through +Inf and NaN values without attempting exponent reduction.
    root = x + x
  elseif (x == 0.0) then
    root = x
  elseif (n == 1) then
    root = x
  else
    call rescale_nth_root(x, n, xr, e_x)

    ! This first estimate is centered in the possible root interval [0.5, 1).
    root_xr = 0.7071067811865475
    rn_minus = real(n - 1)
    rn_plus = real(n + 1)

    do itt=1,halley_iterations_by_n(n)
      root_n = integer_power(root_xr, n)
      num = (rn_minus * root_n) + (rn_plus * xr)
      den = (rn_plus * root_n) + (rn_minus * xr)
      root_xr = root_xr * (num / den)
    enddo

    ! One Newton iteration with a compensated residual polishes the result
    ! after Halley convergence.
    root_nm1 = integer_power(root_xr, n - 1)
    root_xr = root_xr + (power_residual(xr, root_xr, n) / (real(n) * root_nm1))

    root = descale_nth_root(root_xr, e_x)
  endif
end procedure nth_root


!> Rescale x to the range [2**(-n), 1) and compute its nth-root exponent.
pure subroutine rescale_nth_root(x, n, r, e_a)
  real, intent(in) :: x
    !< The number to be rescaled for nth-root computation [A^n]
  integer, intent(in) :: n
    !< The degree of the root [nondim]
  real, intent(out) :: r
    !< The rescaled value of x in the range [2**(-n), 1) [B^n]
  integer(kind=int64), intent(out) :: e_a
    !< The integral component of the nth-root exponent of x [nondim]

  integer(kind=int64) :: xb
    !< Floating point integer representation of x [nondim]
  integer(kind=int64) :: e_x
    !< Exponent of `a` [nondim]
  integer(kind=int64) :: e_r
    !< Exponent of `x` [nondim]
  integer(kind=int64) :: e_shift
    !< Normalizing exponent shift applied to subnormal inputs [nondim]
  integer(kind=int64) :: n64
    !< The root degree promoted to the exponent integer kind [nondim]

  xb = transfer(x, 1_int64)
  e_shift = 0_int64

  ! Rescale subnormal inputs to normal form
  if (ibits(xb, expbit, expwidth) == 0_int64) then
    e_shift = int(digits(real_mold), int64)
    xb = transfer(scale(x, int(e_shift)), 1_int64)
  endif

  e_x = ibits(xb, expbit, expwidth) - expbias
  e_x = e_x - e_shift
  n64 = int(n, int64)

  ! Use floor(e_x/n) + 1 so that the residual exponent is in {-n,...,-1}.
  if (e_x >= 0_int64) then
    e_a = (e_x + n64) / n64
  else
    e_a = (e_x + 1_int64) / n64
  endif
  e_r = e_x - e_a * n64

  ! Insert the new exponent and clear the sign bit so x is positive.
  call mvbits(e_r + expbias, 0, expwidth + 1, xb, expbit)
  r = transfer(xb, 1.)
end subroutine rescale_nth_root


!> Undo the rescaling of a real number back to its original base.
pure function descale_nth_root(x, e_a) result(a)
  real, intent(in) :: x
    !< The rescaled value which is to be restored in ambiguous units [B]
  integer(kind=int64), intent(in) :: e_a
    !< Exponent of the unscaled value [nondim]
  real :: a
    !< Restored value with the corrected exponent in arbitrary units [A]

  integer(kind=int64) :: xb
    !< Bit-packed real number into integer form [nondim]
  integer(kind=int64) :: e_x
    !< Biased exponent of x [nondim]

  xb = transfer(x, 1_int64)
  e_x = ibits(xb, expbit, expwidth)
  call mvbits(e_a + e_x, 0, expwidth, xb, expbit)
  a = transfer(xb, 1.)
end function descale_nth_root


!> Return true if a real value is an IEEE Inf or NaN value.
pure function is_nonfinite(x) result(nonfinite)
  real, intent(in) :: x
    !< A value to test [nondim]
  logical :: nonfinite
    !< True if x is an IEEE Inf or NaN value

  integer(kind=int64) :: xb
    !< Bit-packed real number into integer form [nondim]

  xb = transfer(x, 1_int64)
  nonfinite = (ibits(xb, expbit, expwidth) == (ishft(1_int64, expwidth) - 1_int64))
end function is_nonfinite


!> Compute x - root**n with a compensated product for root**n.
pure function power_residual(x, root, n) result(resid)
  real, intent(in) :: x
    !< Value whose nth root is being estimated [A^n]
  real, intent(in) :: root
    !< Root estimate [A]
  integer, intent(in) :: n
    !< The degree of the root [nondim]
  real :: resid
    !< Compensated residual, x - root**n [A^n]

  real :: p_hi
    !< High part of the accumulated product [A^m]
  real :: p_lo
    !< Low part of the accumulated product [A^m]
  real :: p_new
    !< Rounded high part of the next product [A^m]
  real :: p_err
    !< Low-order residual of the next product [A^m]
  integer :: m
    !< Power counter [nondim]

  p_hi = 1.0
  p_lo = 0.0
  do m=1,n
    p_new = p_hi * root
    p_err = ieee_fma(p_hi, root, -p_new) + (p_lo * root)
    p_hi = p_new + p_err
    p_lo = p_err - (p_hi - p_new)
  enddo

  resid = (x - p_hi) - p_lo
end function power_residual


!> Raise a real number to a small positive integer power with ordered products.
pure function integer_power(x, n) result(xn)
  real, intent(in) :: x
    !< The value to raise to an integer power [A]
  integer, intent(in) :: n
    !< The exponent of x [nondim]
  real :: xn
    !< x raised to the nth power [A^n]

  integer :: m
    !< Power counter [nondim]

  select case (n)
  case (0)
    xn = 1.0
  case (1)
    xn = x
  case (2)
    xn = x * x
  case (3)
    xn = x * x * x
  case (4)
    xn = (x * x) * (x * x)
  case (5)
    xn = ((x * x) * (x * x)) * x
  case default
    xn = 1.0
    do m=1,n
      xn = xn * x
    enddo
  end select
end function integer_power

end submodule MOM_nth_root
