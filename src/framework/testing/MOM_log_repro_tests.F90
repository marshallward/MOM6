! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Unit tests for log_repro
module MOM_log_repro_tests

use, intrinsic :: ieee_arithmetic, only : ieee_value
use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_signaling_nan
use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

use MOM_error_handler, only : assert
use MOM_unit_testing, only : TestSuite
use MOM_intrinsic_functions, only : exp_repro, log_repro

implicit none ; private

public :: add_log_repro_tests

integer, parameter :: realquad = selected_real_kind(p=30, r=300)
  !< Potential real128 precision kind.  If unavailable, this will be negative.
integer, parameter :: realq = merge(realquad, kind(1.), realquad >= 0.)
  !< Placeholder real128 for declarations.  Unused if real128 is unavailable.

real, parameter :: rmold = 0.
  !< Mold value for default real [nondim]

! Mathematical constants
real, parameter :: ln2 = 0.69314718055994531
  !< ln(2) [nondim]
real, parameter :: half_ln2 = 0.34657359027997265
  !< ln(2)/2 [nondim]
real, parameter :: sqrt2 = 1.4142135623730950
  !< sqrt(2) [nondim]
real, parameter :: log_huge = 709.78271289338397
  !< log(huge()) [nondim]
real, parameter :: log_tiny = -708.39641853226408
  !< log(tiny()) [nondim]
real, parameter :: log_033 = -1.1086626245216111
  !< log(0.33) [nondim]
real, parameter :: log_2 = 0.69314718055994531
  !< log(2) [nondim]
real, parameter :: log_10 = 2.3025850929940457
  !< log(10) [nondim]
real, parameter :: log_smallest_subnormal = -744.44007192138122
  !< log(nearest(0, 1)) [nondim]
real, parameter :: log_one_plus_eps = 2.2204460492503128e-16
  !< log(1 + epsilon(1)) [nondim]
real, parameter :: log_one_minus_eps = -2.2204460492503136e-16
  !< log(1 - epsilon(1)) [nondim]

contains

!> Test that log_repro is elemental (works on arrays)
subroutine test_log_elemental
  real :: x(4), val(4), ref(4), err
  real, parameter :: tol = 1.e-14
  integer :: i

  x(:) = [0.5, 1., 2., 10.]
  ref(:) = [-ln2, 0., ln2, log_10]

  val = log_repro(x)

  do i = 1, 4
    if (ref(i) /= 0.) then
      err = abs(val(i) - ref(i)) / abs(ref(i))
    else
      err = abs(val(i) - ref(i))
    endif
    call assert(err < tol, "log_repro elemental test failed")
  enddo
end subroutine test_log_elemental

!> Test that log_repro(1) = 0 exactly
subroutine test_log_one
  real :: val

  val = log_repro(1.)

  print '(2x, "log_repro(1) = ", ES22.15)', val

  call assert(val == 0., "log_repro(1) should equal 0 exactly")
end subroutine test_log_one

!> Test log_repro at a general value (0.33)
subroutine test_log_general
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-12

  x = 0.33
  ref = log_033
  val = log_repro(x)
  err = abs(val - ref) / abs(ref)

  print '(2x, "log_repro(0.33) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(0.33) relative error exceeds tolerance")
end subroutine test_log_general

!> Test log_repro(2) against ln(2)
subroutine test_log_two
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = log_repro(2.)
  err = abs(val - log_2) / log_2

  print '(2x, "log_repro(2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(2) relative error exceeds tolerance")
end subroutine test_log_two

!> Test log_repro at range reduction boundary: 0.5
subroutine test_log_half
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = log_repro(0.5)
  err = abs(val + ln2) / ln2

  print '(2x, "log_repro(0.5) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(0.5) should be close to -ln(2)")
end subroutine test_log_half

!> Test log_repro at half the range reduction boundary: sqrt(2)
subroutine test_log_sqrt2
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = log_repro(sqrt2)
  err = abs(val - half_ln2) / half_ln2

  print '(2x, "log_repro(sqrt2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(sqrt(2)) should be close to ln(2)/2")
end subroutine test_log_sqrt2

!> Test log_repro at half the range reduction boundary: 1 / sqrt(2)
subroutine test_log_inv_sqrt2
  real :: val, err
  real, parameter :: tol = 1.e-15

  val = log_repro(1. / sqrt2)
  err = abs(val + half_ln2) / half_ln2

  print '(2x, "log_repro(1/sqrt2) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(1/sqrt(2)) should be close to -ln(2)/2")
end subroutine test_log_inv_sqrt2

!> Test log_repro(10) against reference value
subroutine test_log_ten
  real :: val, err
  real, parameter :: tol = 1.e-12

  val = log_repro(10.)
  err = abs(val - log_10) / log_10

  print '(2x, "log_repro(10) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(10) relative error exceeds tolerance")
end subroutine test_log_ten

!> Test log_repro(+Inf) = +Inf
subroutine test_log_pos_inf
  real :: x, val, ref

  x = ieee_value(0., ieee_positive_inf)
  ref = ieee_value(0., ieee_positive_inf)
  val = log_repro(x)

  print '(2x, "log_repro(+Inf) = ", ES22.15)', val

  call assert(val == ref, "log_repro(+Inf) should equal +Inf")
end subroutine test_log_pos_inf

!> Test log_repro(0) = -Inf
subroutine test_log_zero
  real :: val, ref

  ref = ieee_value(0., ieee_negative_inf)
  val = log_repro(0.)

  print '(2x, "log_repro(0) = ", ES22.15)', val

  call assert(val == ref, "log_repro(0) should equal -Inf")
end subroutine test_log_zero

!> Test log_repro(-0) = -Inf
subroutine test_log_neg_zero
  real :: val, ref

  ref = ieee_value(0., ieee_negative_inf)
  val = log_repro(-0.)

  print '(2x, "log_repro(-0) = ", ES22.15)', val

  call assert(val == ref, "log_repro(-0) should equal -Inf")
end subroutine test_log_neg_zero

!> Test log_repro(-Inf) = NaN
subroutine test_log_neg_inf
  real :: x, val

  x = ieee_value(0., ieee_negative_inf)
  val = log_repro(x)

  print '(2x, "log_repro(-Inf) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(-Inf) should be NaN")
end subroutine test_log_neg_inf

!> Test log_repro(negative) = NaN
subroutine test_log_negative
  real :: val

  val = log_repro(-1.)

  print '(2x, "log_repro(-1) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(-1) should be NaN")
end subroutine test_log_negative

!> Test log_repro(NaN) = NaN
subroutine test_log_nan
  real :: x, val

  x = ieee_value(0., ieee_quiet_nan)
  val = log_repro(x)

  print '(2x, "log_repro(NaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(NaN) should be NaN")
end subroutine test_log_nan

!> Test log_repro(-NaN) = NaN
subroutine test_log_neg_nan
  real :: x, val

  x = -ieee_value(0., ieee_quiet_nan)
  val = log_repro(x)

  print '(2x, "log_repro(-NaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(-NaN) should be NaN")
end subroutine test_log_neg_nan

!> Test log_repro(sNaN) = NaN
subroutine test_log_signaling_nan
  real :: x, val

  x = ieee_value(0., ieee_signaling_nan)
  val = log_repro(x)

  print '(2x, "log_repro(sNaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(sNaN) should be NaN")
end subroutine test_log_signaling_nan

!> Test log_repro(-sNaN) = NaN
subroutine test_log_neg_signaling_nan
  real :: x, val

  x = -ieee_value(0., ieee_signaling_nan)
  val = log_repro(x)

  print '(2x, "log_repro(-sNaN) = ", ES22.15)', val

  call assert(ieee_is_nan(val), "log_repro(-sNaN) should be NaN")
end subroutine test_log_neg_signaling_nan

!> Test log_repro at largest representable float boundary
subroutine test_log_largest_float
  real :: val, ref, err
  real, parameter :: tol = 1.e-15

  ref = log_huge
  val = log_repro(huge(1.))
  err = abs(val - ref) / ref

  print '(2x, "log_repro(huge) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(huge) relative error exceeds tolerance")
end subroutine test_log_largest_float

!> Test log_repro at smallest normal float boundary
subroutine test_log_smallest_float
  real :: val, ref, err
  real, parameter :: tol = 1.e-15

  ref = log_tiny
  val = log_repro(tiny(1.))
  err = abs(val - ref) / abs(ref)

  print '(2x, "log_repro(tiny) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(tiny) relative error exceeds tolerance")
end subroutine test_log_smallest_float

!> Test log_repro at smallest positive subnormal float boundary
subroutine test_log_smallest_subnormal
  real :: val, ref, err
  real, parameter :: tol = 1.e-15

  ref = log_smallest_subnormal
  val = log_repro(nearest(0., 1.))
  err = abs(val - ref) / abs(ref)

  print '(2x, "log_repro(nearest(0,1)) = ", ES22.15, ", rel err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(nearest(0,1)) relative error exceeds tolerance")
end subroutine test_log_smallest_subnormal

!> Test log_repro near 1 from above
subroutine test_log_near_one_above
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-30

  x = 1. + epsilon(1.)
  ref = log_one_plus_eps
  val = log_repro(x)
  err = abs(val - ref)

  print '(2x, "log_repro(1+eps) = ", ES22.15, ", abs err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(1+eps) absolute error exceeds tolerance")
end subroutine test_log_near_one_above

!> Test log_repro near 1 from below
subroutine test_log_near_one_below
  real :: x, val, ref, err
  real, parameter :: tol = 1.e-30

  x = 1. - epsilon(1.)
  ref = log_one_minus_eps
  val = log_repro(x)
  err = abs(val - ref)

  print '(2x, "log_repro(1-eps) = ", ES22.15, ", abs err = ", ES9.2)', val, err

  call assert(err < tol, "log_repro(1-eps) absolute error exceeds tolerance")
end subroutine test_log_near_one_below

!> Test ULP accuracy over a wide range of positive values.
!!
!! log_repro should be within a few ULP of the true value.
subroutine test_log_ulp_accuracy
  integer, parameter :: npts = 100000
  real, parameter :: ymin = -700.
  real, parameter :: ymax = 700.

  ! Input axis
  real :: y(npts), x(npts)
  real :: I_npts

  real :: val(npts), val_vec(npts)
  real :: val_log(npts), val_log_vec(npts)
  real(kind=realq) :: val_quad(npts), val_quad_vec(npts)

  integer :: i

  ! Generate positive test points with logarithms spanning [ymin,ymax].
  I_npts = 1. / (npts - 1)
  do i = 1, npts
    y(i) = ymin + (i - 1) * ((ymax - ymin) * I_npts)
    x(i) = exp(y(i))
  enddo

  ! Several libraries have scalar and vector implementations, chosen at the
  ! discretion of the compiler.  The following attempts to test each case.

  ! Scalar evaluations
  do i = 1, npts
    val_log(i) = log(x(i))
    val(i) = log_repro(x(i))
    val_quad(i) = log(real(x(i), realq))

    ! Impossible branch to prevent vectorization
    if (val(i) > huge(val(i))) exit
  enddo

  ! Vector-favorable evaluation
  val_log_vec = log(x)
  val_vec = log_repro(x)
  val_quad_vec = log(real(x, realq))

  ! Assert that log_repro() is within 8 ULP.
  print '(1x,a)', '=== scalar log_repro() accuracy'
  call check_ulp_accuracy(x, val, val_quad, max_ulp_tol=8.)

  ! We expect scalar and vector implementations to agree.
  call assert(all(val == val_vec), 'Scalar and vector log_repro() do not agree')
  print '(1x,a)', '=== vector log_repro() matches scalar'

  ! log() accuracy is provided for comparison.
  print '(1x,a)', '=== scalar log() accuracy'
  call check_ulp_accuracy(x, val_log, val_quad)

  if (all(val_log == val_log_vec)) then
    print '(1x,a)', '=== vector log() matches scalar'
  else
    print '(1x,a)', '=== vector log() accuracy'
    call check_ulp_accuracy(x, val_log_vec, val_quad_vec)
  endif
end subroutine test_log_ulp_accuracy

!> Compute the function accuracy relative to a real128-precision reference.
!! Absolute, relative, and ULP error is computed, as well as the number of
!! points above 0.5 and 1 ULP.  An error is raised if any point exceeds the
!! optionally prescribed maximum ULP tolerance.
subroutine check_ulp_accuracy(x, val, ref, max_ulp_tol)
  real, intent(in) :: x(:)
    !< Input grid
  real, intent(in) :: val(:)
    !< Output estimates
  real(kind=realq), intent(in) :: ref(:)
    !< Reference estimates in real128 precision
  real, optional, intent(in) :: max_ulp_tol
    !< Maximum ULP tolerance

  real(kind=realq) :: err, rel_err, ulp_err
    ! Absolute, relative, and ULP error of val relative to ref.

  real :: max_abs_err, max_rel_err, max_ulp
  real :: sum_abs_err, sum_rel_err, sum_sq_err
  real :: x_max_abs, x_max_rel, x_max_ulp
  real :: ulp_val

  integer :: count_exact, count_half_ulp, count_one_ulp
  integer :: i, npts

  npts = size(x)

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

  do i = 1, npts
    ! Absolute error
    err = abs(real(val(i), realq) - ref(i))

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

      if (rel_err > max_rel_err) then
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

  print '(2x,"Tested ", i0, " points")', npts

  ! NOTE: Floats use t25 assuming a positive sign as blank
  print '(2x,"max abs err:", t25, ES12.5, " at input = ", ES14.5E3)', &
      max_abs_err, x_max_abs
  print '(2x,"max rel err:", t25, ES12.5, " at input = ", ES14.5E3)', &
      max_rel_err, x_max_rel
  print '(2x,"max ULP err (vs quad):", t26, f12.10, " at input = ", ES14.5E3)', &
      max_ulp, x_max_ulp
  print '(2x,"mean abs err:", t25, ES12.5)', sum_abs_err / npts
  print '(2x,"mean rel err:", t25, ES12.5)', sum_rel_err / npts
  print '(2x,"RMS err:", t25, ES12.5)', sqrt(sum_sq_err / npts)

  print '(2x,"correct (<0.5 ULP):", t25, i0, 1x, "(", f6.2, "%)")', &
      count_exact, 100. * count_exact / npts
  print '(2x,"above 0.5 ULP:", t26, i0, 1x, "(", f6.2, "%)")', &
      count_half_ulp, 100. * count_half_ulp / npts
  print '(2x,"above 1 ULP:", t26, i0, 1x, "(", f6.2, "%)")', &
      count_one_ulp, 100. * count_one_ulp / npts

  ! exp_repro should be within 2 ULP of the quad-precision reference.
  if (present(max_ulp_tol)) then
    call assert(max_ulp < max_ulp_tol, "Max ULP error exceeds tolerance")
  endif
end subroutine check_ulp_accuracy

!> Test the logarithmic property: log(x*y) = log(x) + log(y)
!!
!! This property test verifies that log_repro produces consistent results
!! by checking the logarithm of a product against the sum of logarithms,
!! within floating-point tolerance.
subroutine test_log_product_property
  integer, parameter :: npts = 1000
  real, parameter :: tol = 1.e-13

  real :: r1, r2
  real :: x1, x2
  real :: log_product, log_sum
  real :: err, max_err
  real :: r1_max, r2_max
  integer :: i
  real :: seed

  max_err = 0.

  ! Use a simple deterministic sequence for reproducibility.
  seed = 0.246813579

  do i = 1, npts
    ! Generate pseudo-random logarithms in a range where the product is finite.
    ! Keep r1, r2 in [-5, 5] so x1*x2 is in exp([-10, 10]).
    seed = mod(seed * 1103515245. + 12345., 2.**31)
    r1 = (seed / 2.**31) * 10. - 5.

    seed = mod(seed * 1103515245. + 12345., 2.**31)
    r2 = (seed / 2.**31) * 10. - 5.

    x1 = exp_repro(r1)
    x2 = exp_repro(r2)

    log_product = log_repro(x1 * x2)
    log_sum = log_repro(x1) + log_repro(x2)

    if (log_product /= 0.) then
      err = abs(log_sum - log_product) / abs(log_product)
    else
      err = abs(log_sum - log_product)
    endif

    if (err > max_err) then
      max_err = err
      r1_max = r1
      r2_max = r2
    endif
  enddo

  print '("Tested ", i0, " random positive (x1, x2) pairs")', npts
  print '("max rel err in log(x1*x2) vs log(x1)+log(x2):", t52, ES12.5)', max_err
  print '("  at log(x1) = ", f10.6, ", log(x2) = ", f10.6)', r1_max, r2_max

  call assert(max_err < tol, "log_repro product property test failed")
end subroutine test_log_product_property

!> Test the logarithmic property: log(1/x) = -log(x)
!!
!! This property test verifies that log_repro(1/x) equals the negation of
!! log_repro(x), within floating-point tolerance.
subroutine test_log_reciprocal_property
  integer, parameter :: npts = 1000
  real, parameter :: tol = 1.e-13

  real :: r
  real :: x
  real :: log_inv_x, neg_log_x
  real :: err, max_err
  real :: r_max
  integer :: i
  real :: seed

  max_err = 0.

  ! Use a simple deterministic sequence for reproducibility.
  seed = 0.135792468

  do i = 1, npts
    ! Generate pseudo-random logarithms that avoid overflow and underflow.
    ! Keep r in [-300, 300] so x and 1/x are both representable.
    seed = mod(seed * 1103515245. + 12345., 2.**31)
    r = (seed / 2.**31) * 600. - 300.

    x = exp_repro(r)

    log_inv_x = log_repro(1. / x)
    neg_log_x = -log_repro(x)

    if (neg_log_x /= 0.) then
      err = abs(log_inv_x - neg_log_x) / abs(neg_log_x)
    else
      err = abs(log_inv_x - neg_log_x)
    endif

    if (err > max_err) then
      max_err = err
      r_max = r
    endif
  enddo

  print '("Tested ", i0, " random positive x values")', npts
  print '("max rel err in log(1/x) vs -log(x):", t52, ES12.5)', max_err
  print '("  at log(x) = ", f10.6)', r_max

  call assert(max_err < tol, "log_repro reciprocal property test failed")
end subroutine test_log_reciprocal_property


!> Add log_repro tests to the intrinsic function test suite.
subroutine add_log_repro_tests(suite)
  type(TestSuite), intent(inout) :: suite

  ! Elemental test
  call suite%add(test_log_elemental, "test_log_elemental")

  ! Basic value tests
  call suite%add(test_log_one, "test_log_one")
  call suite%add(test_log_general, "test_log_general")
  call suite%add(test_log_two, "test_log_two")
  call suite%add(test_log_half, "test_log_half")
  call suite%add(test_log_sqrt2, "test_log_sqrt2")
  call suite%add(test_log_inv_sqrt2, "test_log_inv_sqrt2")
  call suite%add(test_log_ten, "test_log_ten")
  call suite%add(test_log_pos_inf, "test_log_pos_inf")
  call suite%add(test_log_zero, "test_log_zero")
  call suite%add(test_log_neg_zero, "test_log_neg_zero")
  call suite%add(test_log_neg_inf, "test_log_neg_inf")
  call suite%add(test_log_negative, "test_log_negative")
  call suite%add(test_log_nan, "test_log_nan")
  call suite%add(test_log_neg_nan, "test_log_neg_nan")
  call suite%add(test_log_signaling_nan, "test_log_signaling_nan")
  call suite%add(test_log_neg_signaling_nan, "test_log_neg_signaling_nan")
  call suite%add(test_log_largest_float, "test_log_largest_float")
  call suite%add(test_log_smallest_float, "test_log_smallest_float")
  call suite%add(test_log_smallest_subnormal, "test_log_smallest_subnormal")
  call suite%add(test_log_near_one_above, "test_log_near_one_above")
  call suite%add(test_log_near_one_below, "test_log_near_one_below")

  ! Log property tests
  call suite%add(test_log_product_property, "test_log_product_property")
  call suite%add(test_log_reciprocal_property, "test_log_reciprocal_property")

  ! Evaluate log_repro error if quad precision is available
  if (realquad >= 0) then
    call suite%add(test_log_ulp_accuracy, "test_log_ulp_accuracy")
  else
    print '(1x,a)', 'Skipping log_repro ULP accuracy test: quad precision is unavailable.'
  endif
end subroutine add_log_repro_tests

end module MOM_log_repro_tests
