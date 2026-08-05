! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Unit tests for MOM_intrinsic_functions module, including exp_repro()
module MOM_intrinsic_functions_tests

use, intrinsic :: iso_fortran_env, only : real128
use, intrinsic :: ieee_arithmetic, only : ieee_value
use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_signaling_nan
use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
use, intrinsic :: ieee_exceptions, only : ieee_set_flag, ieee_get_flag
use, intrinsic :: ieee_exceptions, only : ieee_support_flag
use, intrinsic :: ieee_exceptions, only : ieee_all
use, intrinsic :: ieee_exceptions, only : ieee_invalid, ieee_overflow
use, intrinsic :: ieee_exceptions, only : ieee_underflow, ieee_inexact
use, intrinsic :: ieee_exceptions, only : ieee_divide_by_zero

use MOM_error_handler, only : assert
use MOM_unit_testing, only : TestSuite
use MOM_intrinsic_functions, only : exp_repro

implicit none ; private

public :: run_intrinsic_functions_tests

! Mold value for default real
real, parameter :: rmold = 0.

! Mathematical constants.
real, parameter :: ln2 = 0.693147180559945309417232121458176568
  !< ln(2) = 0.6931471805599453094172321...
real, parameter :: half_ln2 = 0.346573590279972654708616060729088284
  !< ln(2)/2 = 0.3465735902799726547086160...
real, parameter :: sqrt2 = 1.41421356237309504880168872420969808
  !< sqrt(2) = 1.4142135623730950488016887...
real, parameter :: e_val = 2.71828182845904523536028747135266250
  !< e = 2.7182818284590452353602874...
real, parameter :: inv_e = 0.367879441171442321595523770161460867
  !< 1/e = 0.3678794411714423215955237...

! Module-level flag to enable/disable IEEE exception tests
logical :: ieee_flags_supported = .false.
  !< True if the platform supports IEEE exception flag queries

contains

!> Check if IEEE exception flags are supported on this platform
subroutine check_ieee_support
  ieee_flags_supported = ieee_support_flag(ieee_overflow) &
      .and. ieee_support_flag(ieee_underflow) &
      .and. ieee_support_flag(ieee_inexact) &
      .and. ieee_support_flag(ieee_invalid) &
      .and. ieee_support_flag(ieee_divide_by_zero)
end subroutine check_ieee_support


!> Test that exp_repro(0) = 1 exactly
subroutine test_exp_zero
  real :: val

  val = exp_repro(0.)

  print '(a,ES22.15)', "  exp_repro(0) = ", val

  call assert(val == 1., "exp_repro(0) should equal 1 exactly")
end subroutine test_exp_zero


!> Test that exp_repro(-0) = 1 exactly
subroutine test_exp_neg_zero
  real :: val

  val = exp_repro(-0.)

  print '(a,ES22.15)', "  exp_repro(-0) = ", val

  call assert(val == 1., "exp_repro(-0) should equal 1 exactly")
end subroutine test_exp_neg_zero


!> Test exp_repro(1) against reference value
subroutine test_exp_one
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(1.)
  err = abs(val - e_val) / e_val

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(1) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(1) relative error exceeds tolerance")
end subroutine test_exp_one


!> Test exp_repro(-1) against reference value
subroutine test_exp_neg_one
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(-1.)
  err = abs(val - inv_e) / inv_e

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(-1) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(-1) relative error exceeds tolerance")
end subroutine test_exp_neg_one


!> Test exp_repro at a general value (0.33)
subroutine test_exp_general
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-15

  x = 0.33
  ref = exp(x)  ! Use intrinsic as reference for general case
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(0.33) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(0.33) relative error exceeds tolerance")
end subroutine test_exp_general


!> Test exp_repro at range reduction boundary: ln(2)
subroutine test_exp_ln2
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(ln2)
  err = abs(val - 2.0) / 2.0

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(ln2) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(ln(2)) should be close to 2")
end subroutine test_exp_ln2


!> Test exp_repro at range reduction boundary: -ln(2)
subroutine test_exp_neg_ln2
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(-ln2)
  err = abs(val - 0.5) / 0.5

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(-ln2) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(-ln(2)) should be close to 0.5")
end subroutine test_exp_neg_ln2


!> Test exp_repro at half the range reduction boundary: 0.5*ln(2)
subroutine test_exp_half_ln2
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = exp_repro(half_ln2)
  err = abs(val - sqrt2) / sqrt2

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(ln2/2) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(0.5*ln(2)) should be close to sqrt(2)")
end subroutine test_exp_half_ln2


!> Test exp_repro with tiny argument (tests 1 + x accuracy)
subroutine test_exp_tiny_pos
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-15

  x = 1.0e-15
  ref = exp(x)
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(1e-15) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(1e-15) relative error exceeds tolerance")
end subroutine test_exp_tiny_pos


!> Test exp_repro with tiny negative argument
subroutine test_exp_tiny_neg
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-15

  x = -1.0e-15
  ref = exp(x)
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(-1e-15) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(-1e-15) relative error exceeds tolerance")
end subroutine test_exp_tiny_neg


!> Test exp_repro overflow behavior: exp(1000)
subroutine test_exp_overflow
  real :: x, val, ref

  x = 1000.
  ref = exp(x)  ! Should be +Inf
  val = exp_repro(x)

  print '(a,ES22.15)', "  exp_repro(1000) = ", val

  call assert(val == ref, "exp_repro(1000) should overflow to +Inf")
end subroutine test_exp_overflow


!> Test exp_repro underflow behavior: exp(-1000)
subroutine test_exp_underflow
  real :: x, val, ref

  x = -1000.
  ref = exp(x)  ! Should be 0
  val = exp_repro(x)

  print '(a,ES22.15)', "  exp_repro(-1000) = ", val

  call assert(val == ref, "exp_repro(-1000) should underflow to 0")
end subroutine test_exp_underflow


!> Test exp_repro near overflow boundary (709.78)
subroutine test_exp_near_overflow
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-14

  x = 709.78
  ref = exp(x)
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(709.78) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(709.78) relative error exceeds tolerance")
end subroutine test_exp_near_overflow


!> Test exp_repro near underflow boundary (-708.39)
subroutine test_exp_near_underflow
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-14

  x = -708.39
  ref = exp(x)
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(-708.39) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(-708.39) relative error exceeds tolerance")
end subroutine test_exp_near_underflow


!> Test exp_repro in subnormal region (-745.13)
subroutine test_exp_subnormal
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-10  ! Relaxed for subnormals

  x = -745.13
  ref = exp(x)
  val = exp_repro(x)

  print '(a,ES22.15)', "  exp_repro(-745.13) = ", val

  ! For subnormals, check absolute error or that both are very small
  if (ref > 0.) then
    err = abs(val - ref) / ref
    call assert(err < tol, "exp_repro(-745.13) relative error exceeds tolerance")
  else
    call assert(val == 0., "exp_repro(-745.13) should be zero or subnormal")
  endif
end subroutine test_exp_subnormal


!> Test exp_repro at largest representable float boundary
subroutine test_exp_largest_float
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-13  ! Relaxed due to log/exp round-trip error
  real, parameter :: LOG_HUGE = log(huge(rmold))

  x = LOG_HUGE
  ref = exp(x)
  val = exp_repro(x)
  err = abs(val - ref) / ref

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(log(huge)) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(log(huge)) relative error exceeds tolerance")
end subroutine test_exp_largest_float


!> Test exp_repro at smallest normal float boundary
subroutine test_exp_smallest_float
  real :: val, ref, err
  real, parameter :: tol = 5.e-14
  real, parameter :: LOG_TINY = log(tiny(rmold))

  ref = tiny(1.)
  val = exp_repro(LOG_TINY)
  err = abs(val - ref) / ref

  print '(a,ES22.15,a,ES9.2)', "  exp_repro(log(tiny)) = ", val, ", rel err = ", err

  call assert(err < tol, "exp_repro(log(tiny)) relative error exceeds tolerance")
end subroutine test_exp_smallest_float


!> Test exp_repro(+Inf) = +Inf
subroutine test_exp_pos_inf
  real :: x, val, ref

  x = ieee_value(0., ieee_positive_inf)
  ref = ieee_value(0., ieee_positive_inf)
  val = exp_repro(x)

  print '(a,ES22.15)', "  exp_repro(+Inf) = ", val

  call assert(val == ref, "exp_repro(+Inf) should equal +Inf")
end subroutine test_exp_pos_inf


!> Test exp_repro(-Inf) = 0
subroutine test_exp_neg_inf
  real :: x, val

  x = ieee_value(0., ieee_negative_inf)
  val = exp_repro(x)

  print '(a,ES22.15)', "  exp_repro(-Inf) = ", val

  call assert(val == 0., "exp_repro(-Inf) should equal 0")
end subroutine test_exp_neg_inf


!> Test exp_repro(NaN) = NaN
subroutine test_exp_nan
  real :: x, val

  x = ieee_value(0., ieee_quiet_nan)
  val = exp_repro(x)

  print '(a,ES22.15)', "  exp_repro(NaN) = ", val

  call assert(ieee_is_nan(val), "exp_repro(NaN) should be NaN")
end subroutine test_exp_nan


!> Test exp_repro(-NaN) = NaN
subroutine test_exp_neg_nan
  real :: x, val

  x = -ieee_value(0., ieee_quiet_nan)
  val = exp_repro(x)

  print '(a,ES22.15)', "  exp_repro(-NaN) = ", val

  call assert(ieee_is_nan(val), "exp_repro(-NaN) should be NaN")
end subroutine test_exp_neg_nan


!> Test ULP accuracy over a range of values [-10, 10]
!!
!! exp_repro should be within 1.5 ULP of the true value (1 ULP with FMA).
!! We test against 2 ULP to allow for combined error from exp_repro and the
!! intrinsic exp() reference.
subroutine test_exp_ulp_accuracy
  integer, parameter :: npts = 10000
  real, parameter :: xmin = -10.
  real, parameter :: xmax = 10.
  real, parameter :: max_ulp_tol = 2.

  ! Input axis
  real :: x(npts)
  real :: I_npts

  real :: val(npts)
  real(real128) :: ref(npts)

  real(real128) :: err, rel_err, ulp_err
  real :: ulp_val

  real :: max_abs_err, max_rel_err, max_ulp
  real :: sum_abs_err, sum_rel_err, sum_sq_err

  real :: x_max_abs, x_max_rel, x_max_ulp
  integer :: i, count_exact, count_half_ulp, count_one_ulp

  ! Quad precision reference
  real :: max_ulp_quad, ulp_err_quad, x_max_ulp_quad

  ! Generate test points
  I_npts = 1. / (npts - 1)
  do i = 1, npts
    x(i) = xmin + (i - 1) * ((xmax - xmin) * I_npts)
  enddo

  ! Compute test and reference values
  do i = 1, npts
    val(i) = exp_repro(x(i))
    ref(i) = exp(real(x(i), real128))
  enddo

  ! Initialize statistics
  max_abs_err = 0.
  max_rel_err = 0.
  max_ulp = 0.
  sum_abs_err = 0.
  sum_rel_err = 0.
  sum_sq_err = 0.
  count_exact = 0
  count_half_ulp = 0
  count_one_ulp = 0
  max_ulp_quad = 0.

  ! Calculate statistics
  do i = 1, npts
    ! Absolute error

    err = abs(real(val(i), real128) - ref(i))

    sum_abs_err = sum_abs_err + err
    sum_sq_err = sum_sq_err + err * err

    ! Update max error
    if (err > max_abs_err) then
      max_abs_err = err
      x_max_abs = x(i)
    endif

    ! Relative error

    if (ref(i) /= 0.) then
      rel_err = err / abs(ref(i))

      sum_rel_err = sum_rel_err + rel_err

      if (err > max_rel_err) then
        max_rel_err = rel_err
        x_max_rel = x(i)
      endif
    endif

    ! Report error relative to ULP

    ulp_val = spacing(real(ref(i), kind(rmold)))
    ulp_err = abs(val(i) - ref(i)) / ulp_val

    if (ulp_err > max_ulp) then
      max_ulp = ulp_err
      x_max_ulp = x(i)
    endif

    if (ulp_err < 0.5) count_exact = count_exact + 1
    if (ulp_err >= 0.5) count_half_ulp = count_half_ulp + 1
    if (ulp_err >= 1.0) count_one_ulp = count_one_ulp + 1
  enddo

  print '(a,1x,i0,1x,a)', "Tested", npts, "points in [-10, 10]"

  print '(a21,1x,ES12.5,1x,a,1x,f10.4)', "max abs err:", max_abs_err, &
      "at x =", x_max_abs
  print '(a21,1x,ES12.5,1x,a,1x,f10.4)', "max rel err:", max_rel_err, &
      "at x =", x_max_rel
  print '(a21,1x,f12.10,1x,a,1x,f10.4)', "max ULP err (vs exp):", max_ulp, &
      "at x =", x_max_ulp

  print '(a21,1x,ES12.5)', "mean abs err:", sum_abs_err / npts
  print '(a21,1x,ES12.5)', "mean rel err:", sum_rel_err / npts
  print '(a21,1x,ES12.5)', "RMS err:", sqrt(sum_sq_err / npts)

  print '(a21,1x,i0,1x,a,1x,f6.2,a)', "correct (<0.5 ULP):", count_exact, &
      "(", 100. * count_exact / npts, "%)"
  print '(a21,1x,i0,1x,a,1x,f6.2,a)', "above 0.5 ULP:", count_half_ulp, &
      "(", 100. * count_half_ulp / npts, "%)"
  print '(a21,1x,i0,1x,a,1x,f6.2,a)', "above 1 ULP:", count_one_ulp, &
      "(", 100. * count_one_ulp / npts, "%)"

  ! exp_repro should be within 2 ULP of intrinsic exp
  call assert(max_ulp < max_ulp_tol, "exp_repro max ULP error exceeds 2 over [-10, 10]")
end subroutine test_exp_ulp_accuracy


!> Test that exp_repro is elemental (works on arrays)
subroutine test_exp_elemental
  real :: x(5), val(5), ref(5), err
  real, parameter :: tol = 1.e-14
  integer :: i

  x(:) = [-2., -1., 0., 1., 2.]

  do i = 1, 5
    ref(i) = exp(x(i))
    val(i) = exp_repro(x(i))
  enddo

  do i = 1, 5
    if (ref(i) /= 0.) then
      err = abs(val(i) - ref(i)) / abs(ref(i))
      call assert(err < tol, "exp_repro elemental test failed")
    endif
  enddo
end subroutine test_exp_elemental


!> Test IEEE exception flags for normal input (should raise inexact only)
subroutine test_exp_flags_normal
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = 0.33

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(0.33) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(0.33) should not raise overflow")
  call assert(.not. flag_underflow, "exp_repro(0.33) should not raise underflow")
  call assert(flag_inexact, "exp_repro(0.33) should raise inexact")
end subroutine test_exp_flags_normal


!> Test IEEE exception flags for overflow (should raise overflow + inexact)
subroutine test_exp_flags_overflow
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = 1000.

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(1000) should not raise invalid")
  call assert(flag_overflow, "exp_repro(1000) should raise overflow")
  call assert(.not. flag_underflow, "exp_repro(1000) should not raise underflow")
  call assert(flag_inexact, "exp_repro(1000) should raise inexact")
end subroutine test_exp_flags_overflow


!> Test IEEE exception flags for underflow (should raise underflow + inexact)
subroutine test_exp_flags_underflow
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = -1000.

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(-1000) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(-1000) should not raise overflow")
  call assert(flag_underflow, "exp_repro(-1000) should raise underflow")
  call assert(flag_inexact, "exp_repro(-1000) should raise inexact")
end subroutine test_exp_flags_underflow


!> Test IEEE exception flags for +Inf input (should raise no flags)
subroutine test_exp_flags_pos_inf
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = ieee_value(0., ieee_positive_inf)

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(+Inf) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(+Inf) should not raise overflow")
  call assert(.not. flag_underflow, "exp_repro(+Inf) should not raise underflow")
  call assert(.not. flag_inexact, "exp_repro(+Inf) should not raise inexact")
end subroutine test_exp_flags_pos_inf


!> Test IEEE exception flags for -Inf input (should raise no flags)
subroutine test_exp_flags_neg_inf
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = ieee_value(0., ieee_negative_inf)

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(-Inf) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(-Inf) should not raise overflow")
  call assert(.not. flag_underflow, "exp_repro(-Inf) should not raise underflow")
  call assert(.not. flag_inexact, "exp_repro(-Inf) should not raise inexact")
end subroutine test_exp_flags_neg_inf


!> Test IEEE exception flags for NaN input (should raise no flags for quiet NaN)
subroutine test_exp_flags_nan
  real, volatile :: x, val
  logical :: flag_invalid, flag_overflow, flag_underflow, flag_inexact

  if (.not. ieee_flags_supported) return

  x = ieee_value(0., ieee_quiet_nan)

  call ieee_set_flag(ieee_all, .false.)
  val = exp_repro(x)
  call ieee_get_flag(ieee_invalid, flag_invalid)
  call ieee_get_flag(ieee_overflow, flag_overflow)
  call ieee_get_flag(ieee_underflow, flag_underflow)
  call ieee_get_flag(ieee_inexact, flag_inexact)

  call assert(.not. flag_invalid, "exp_repro(NaN) should not raise invalid")
  call assert(.not. flag_overflow, "exp_repro(NaN) should not raise overflow")
  call assert(.not. flag_underflow, "exp_repro(NaN) should not raise underflow")
  call assert(.not. flag_inexact, "exp_repro(NaN) should not raise inexact")
end subroutine test_exp_flags_nan


!> Print a summary table of IEEE exception flags for various inputs
!!
!! This is informational only - it shows which flags are raised by exp() and
!! exp_repro() for different input categories. Not enforced by assertions.
subroutine print_ieee_flags_summary
  real, volatile :: x, val_intrinsic, val_repro
  logical :: flags_intrinsic(5), flags_repro(5)
  character(len=5) :: flag_str_intrinsic, flag_str_repro
  integer :: i

  ! Test cases: name, input value
  character(len=20) :: test_names(12)
  real :: test_values(12)

  if (.not. ieee_flags_supported) then
    print '(a)', "  IEEE flags not supported on this platform"
    return
  endif

  ! Define test cases
  test_names(1) = "exact (0)"; test_values(1) = 0.
  test_names(2) = "normal"; test_values(2) = 0.33
  test_names(3) = "overflow"; test_values(3) = 1000.
  test_names(4) = "underflow"; test_values(4) = -1000.
  test_names(5) = "near overflow"; test_values(5) = 709.78
  test_names(6) = "near underflow"; test_values(6) = -708.39
  test_names(7) = "largest float"; test_values(7) = log(huge(rmold))
  test_names(8) = "smallest normal"; test_values(8) = log(tiny(rmold))
  test_names(9) = "+Inf"; test_values(9) = ieee_value(0., ieee_positive_inf)
  test_names(10) = "-Inf"; test_values(10) = ieee_value(0., ieee_negative_inf)
  test_names(11) = "NaN"; test_values(11) = ieee_value(0., ieee_quiet_nan)
  test_names(12) = "sNaN"; test_values(12) = ieee_value(0., ieee_signaling_nan)

  print '(/a,2x,a,3x,a)', "IEEE flags summary:", "exp()", "exp_repro"
  print '(a18,2x,a5,3x,a5)', "", "IOUXZ", "IOUXZ"

  do i = 1, 12
    x = test_values(i)

    ! Test intrinsic exp()
    call ieee_set_flag(ieee_all, .false.)
    val_intrinsic = exp(x)
    call get_flag_string(flags_intrinsic, flag_str_intrinsic)

    ! Test exp_repro()
    call ieee_set_flag(ieee_all, .false.)
    val_repro = exp_repro(x)
    call get_flag_string(flags_repro, flag_str_repro)

    print '(a18,a,a5,a,a5)', trim(test_names(i)), ": ", flag_str_intrinsic, "   ", flag_str_repro
  enddo

contains

  subroutine get_flag_string(flags, str)
    logical, intent(out) :: flags(5)
    character(len=5), intent(out) :: str

    call ieee_get_flag(ieee_invalid, flags(1))
    call ieee_get_flag(ieee_overflow, flags(2))
    call ieee_get_flag(ieee_underflow, flags(3))
    call ieee_get_flag(ieee_inexact, flags(4))
    call ieee_get_flag(ieee_divide_by_zero, flags(5))

    str(1:1) = merge('I', '.', flags(1))
    str(2:2) = merge('O', '.', flags(2))
    str(3:3) = merge('U', '.', flags(3))
    str(4:4) = merge('X', '.', flags(4))
    str(5:5) = merge('Z', '.', flags(5))
  end subroutine get_flag_string

end subroutine print_ieee_flags_summary


!> Run all intrinsic function tests
subroutine run_intrinsic_functions_tests
  type(TestSuite) :: suite

  ! Check IEEE exception support before running tests
  call check_ieee_support

  suite = TestSuite()

  ! Basic value tests
  call suite%add(test_exp_zero, "test_exp_zero")
  call suite%add(test_exp_neg_zero, "test_exp_neg_zero")
  call suite%add(test_exp_one, "test_exp_one")
  call suite%add(test_exp_neg_one, "test_exp_neg_one")
  call suite%add(test_exp_general, "test_exp_general")

  ! Range reduction boundary tests
  call suite%add(test_exp_ln2, "test_exp_ln2")
  call suite%add(test_exp_neg_ln2, "test_exp_neg_ln2")
  call suite%add(test_exp_half_ln2, "test_exp_half_ln2")

  ! Tiny argument tests
  call suite%add(test_exp_tiny_pos, "test_exp_tiny_pos")
  call suite%add(test_exp_tiny_neg, "test_exp_tiny_neg")

  ! Overflow/underflow tests
  call suite%add(test_exp_overflow, "test_exp_overflow")
  call suite%add(test_exp_underflow, "test_exp_underflow")
  call suite%add(test_exp_near_overflow, "test_exp_near_overflow")
  call suite%add(test_exp_near_underflow, "test_exp_near_underflow")
  call suite%add(test_exp_subnormal, "test_exp_subnormal")

  ! Boundary tests
  call suite%add(test_exp_largest_float, "test_exp_largest_float")
  call suite%add(test_exp_smallest_float, "test_exp_smallest_float")

  ! Special value tests
  call suite%add(test_exp_pos_inf, "test_exp_pos_inf")
  call suite%add(test_exp_neg_inf, "test_exp_neg_inf")
  call suite%add(test_exp_nan, "test_exp_nan")
  call suite%add(test_exp_neg_nan, "test_exp_neg_nan")

  ! Evaluate error if quad precision is available
  if (real128 >= 0) then
    call suite%add(test_exp_ulp_accuracy, "test_exp_ulp_accuracy")
  endif

  ! Elemental test
  call suite%add(test_exp_elemental, "test_exp_elemental")

  ! IEEE exception flag tests (skipped if not supported)
  call suite%add(test_exp_flags_normal, "test_exp_flags_normal")
  call suite%add(test_exp_flags_overflow, "test_exp_flags_overflow")
  call suite%add(test_exp_flags_underflow, "test_exp_flags_underflow")
  call suite%add(test_exp_flags_pos_inf, "test_exp_flags_pos_inf")
  call suite%add(test_exp_flags_neg_inf, "test_exp_flags_neg_inf")
  call suite%add(test_exp_flags_nan, "test_exp_flags_nan")

  call suite%run()

  ! Print IEEE flags summary table (informational, not enforced)
  call print_ieee_flags_summary()
end subroutine run_intrinsic_functions_tests

end module MOM_intrinsic_functions_tests
