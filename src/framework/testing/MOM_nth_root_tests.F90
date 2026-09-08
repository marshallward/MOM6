! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Unit tests for nth_root()
module MOM_nth_root_tests

use, intrinsic :: ieee_arithmetic, only : ieee_value
use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

use MOM_error_handler, only : assert
use MOM_intrinsic_functions, only : nth_root
use MOM_unit_testing, only : TestSuite

implicit none ; private

public :: add_nth_root_tests

integer, parameter :: realquad = selected_real_kind(p=30, r=300)
  !< Potential real128 precision kind.  If unavailable, this will be negative.
integer, parameter :: realq = merge(realquad, kind(1.), realquad >= 0.)
  !< Placeholder real128 for declarations.  Unused if real128 is unavailable.

real, parameter :: rmold = 0.
  !< Mold value for default real [nondim]
real, parameter :: fifth_root_3 = 1.2457309396155174
  !< The fifth root of 3 [nondim]

contains

!> Add all nth_root tests to a test suite.
subroutine add_nth_root_tests(suite)
  type(TestSuite), intent(inout) :: suite

  call suite%add(test_nth_root_values, "test_nth_root_values")
  call suite%add(test_nth_root_self_consistency, "test_nth_root_self_consistency")
  if (realquad >= 0) then
    call suite%add(test_nth_root_ulp_accuracy, "test_nth_root_ulp_accuracy")
  else
    print '(1x,a)', 'Skipping nth_root ULP accuracy test: quad precision is unavailable.'
  endif
  call suite%add(test_nth_root_rescaling, "test_nth_root_rescaling")
  call suite%add(test_nth_root_special, "test_nth_root_special")
end subroutine add_nth_root_tests


!> Test nth_root against known values.
subroutine test_nth_root_values
  real :: val
    !< A computed nth root [nondim]
  real :: err
    !< A relative nth root error [nondim]
  real, parameter :: tol = 2.e-15
    !< Error tolerance for nth_root value tests [nondim]

  val = nth_root(16.0, 4)
  err = abs(val - 2.0) / 2.0
  print '(2x, "nth_root(16,4) = ", ES22.15, ", rel err = ", ES9.2)', val, err
  call assert(err < tol, "nth_root(16,4) relative error exceeds tolerance")

  val = nth_root(3.0, 5)
  err = abs(val - fifth_root_3) / fifth_root_3
  print '(2x, "nth_root(3,5) = ", ES22.15, ", rel err = ", ES9.2)', val, err
  call assert(err < tol, "nth_root(3,5) relative error exceeds tolerance")
end subroutine test_nth_root_values


!> Test nth_root self consistency across all supported root degrees.
subroutine test_nth_root_self_consistency
  real :: x
    !< Value whose nth root is being tested [nondim]
  real :: root
    !< The computed nth root of x [nondim]
  real :: root_n
    !< The nth power of the computed root [nondim]
  real :: err
    !< A relative nth root reconstruction error [nondim]
  real, parameter :: tol = 1.e-13
    !< Error tolerance for nth_root self consistency [nondim]
  integer :: n
    !< Root degree [nondim]
  integer :: m
    !< Power counter [nondim]

  do n=2,32
    x = 0.03125 + (0.0625 * real(n))
    root = nth_root(x, n)
    root_n = 1.0
    do m=1,n
      root_n = root_n * root
    enddo
    err = abs(root_n - x) / x
    call assert(err < tol, "nth_root self-consistency error exceeds tolerance")
  enddo
end subroutine test_nth_root_self_consistency


!> Test nth_root ULP accuracy against a quad-precision reference.
subroutine test_nth_root_ulp_accuracy
  integer, parameter :: npts = 100000
    !< Number of test points [nondim]
  real, parameter :: xmin = 1.e-6
    !< Minimum test value [nondim]
  real, parameter :: xmax = 10.
    !< Maximum test value [nondim]

  real :: x(npts)
    !< Test input values [nondim]
  real :: val(npts)
    !< Scalar nth_root estimates [nondim]
  real :: val_vec(npts)
    !< Vector nth_root estimates [nondim]
  real :: val_pow(npts)
    !< Scalar exponent-operator estimates [nondim]
  real :: val_pow_vec(npts)
    !< Vector exponent-operator estimates [nondim]
  real(kind=realq) :: val_quad(npts)
    !< Quad-precision reference values [nondim]
  real :: I_npts
    !< Inverse number of intervals in the test grid [nondim]
  integer :: i
    !< Test point counter [nondim]
  integer :: n
    !< Root degree [nondim]

  I_npts = 1. / (npts - 1)
  do i = 1,npts
    x(i) = xmin + (i - 1) * ((xmax - xmin) * I_npts)
  enddo

  do n=2,32
    do i = 1,npts
      val(i) = nth_root(x(i), n)
      val_pow(i) = x(i)**(1.0 / real(n))
      val_quad(i) = real(x(i), realq)**(1._realq / real(n, realq))

      ! Impossible branch to prevent vectorization.
      if (val(i) < 0.) exit
    enddo

    val_vec = nth_root(x, n)
    val_pow_vec = x**(1.0 / real(n))

    print '(1x,a,i0,a)', '=== scalar nth_root(x,', n, ') accuracy'
    call check_nth_root_ulp_accuracy(x, val, val_quad, max_ulp_tol=2.)
    call assert(all(val == val_vec), "Scalar and vector nth_root() do not agree")

    print '(1x,a,i0,a)', '=== scalar x**(1/', n, ') accuracy'
    call check_nth_root_ulp_accuracy(x, val_pow, val_quad)
    if (.not. all(val_pow == val_pow_vec)) then
      print '(1x,a,i0,a)', '=== vector x**(1/', n, ') accuracy'
      call check_nth_root_ulp_accuracy(x, val_pow_vec, val_quad)
    endif
  enddo
end subroutine test_nth_root_ulp_accuracy


!> Compute nth_root accuracy relative to a real128-precision reference.
subroutine check_nth_root_ulp_accuracy(x, val, ref, max_ulp_tol)
  real, intent(in) :: x(:)
    !< Input grid [nondim]
  real, intent(in) :: val(:)
    !< Output estimates [nondim]
  real(kind=realq), intent(in) :: ref(:)
    !< Reference estimates in real128 precision [nondim]
  real, optional, intent(in) :: max_ulp_tol
    !< Maximum ULP tolerance [nondim]

  real(kind=realq) :: err, rel_err, ulp_err
    !< Absolute, relative, and ULP error of val relative to ref [nondim]
  real :: max_abs_err, max_rel_err, max_ulp
    !< Maximum absolute, relative, and ULP errors [nondim]
  real :: sum_abs_err, sum_rel_err, sum_sq_err
    !< Sums used to report average errors [nondim]
  real :: x_max_abs, x_max_rel, x_max_ulp
    !< Inputs where the maximum errors occur [nondim]
  real :: ulp_val
    !< ULP spacing of the rounded reference value [nondim]
  integer :: count_exact, count_half_ulp, count_one_ulp
    !< Counters of ULP accuracy categories [nondim]
  integer :: i, npts
    !< Loop counter and number of test points [nondim]

  npts = size(x)
  max_abs_err = 0.
  max_rel_err = 0.
  max_ulp = 0.
  sum_abs_err = 0.
  sum_rel_err = 0.
  sum_sq_err = 0.
  count_exact = 0
  count_half_ulp = 0
  count_one_ulp = 0

  do i = 1,npts
    err = abs(real(val(i), realq) - ref(i))
    sum_abs_err = sum_abs_err + err
    sum_sq_err = sum_sq_err + err * err

    if (err > max_abs_err) then
      max_abs_err = err
      x_max_abs = x(i)
    endif

    rel_err = err / abs(ref(i))
    sum_rel_err = sum_rel_err + rel_err
    if (rel_err > max_rel_err) then
      max_rel_err = rel_err
      x_max_rel = x(i)
    endif

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

  print '(2x,"Tested ", i0, " points in [", ES10.3, ", ", ES10.3, "]")', &
      npts, minval(x), maxval(x)
  print '(2x,"max abs err:", t25, ES12.5, " at x = ", f10.4)', max_abs_err, x_max_abs
  print '(2x,"max rel err:", t25, ES12.5, " at x = ", f10.4)', max_rel_err, x_max_rel
  print '(2x,"max ULP err (vs quad):", t26, f12.10, " at x = ", f10.4)', max_ulp, x_max_ulp
  print '(2x,"mean abs err:", t25, ES12.5)', sum_abs_err / npts
  print '(2x,"mean rel err:", t25, ES12.5)', sum_rel_err / npts
  print '(2x,"RMS err:", t25, ES12.5)', sqrt(sum_sq_err / npts)
  print '(2x,"correct (<0.5 ULP):", t25, i0, 1x, "(", f6.2, "%)")', &
      count_exact, 100. * count_exact / npts
  print '(2x,"above 0.5 ULP:", t26, i0, 1x, "(", f6.2, "%)")', &
      count_half_ulp, 100. * count_half_ulp / npts
  print '(2x,"above 1 ULP:", t26, i0, 1x, "(", f6.2, "%)")', &
      count_one_ulp, 100. * count_one_ulp / npts

  if (present(max_ulp_tol)) then
    call assert(max_ulp < max_ulp_tol, "nth_root max ULP error exceeds tolerance")
  endif
end subroutine check_nth_root_ulp_accuracy


!> Test nth_root dimensional rescaling by powers of 2**n.
subroutine test_nth_root_rescaling
  real :: base
    !< Base value for the rescaling test [nondim]
  real :: ref
    !< Reference value before rescaling [nondim]
  real :: val
    !< Value after rescaling [nondim]
  real :: scaled
    !< Base value rescaled by an integer power of 2**n [nondim]
  integer :: n
    !< Root degree [nondim]
  integer :: k
    !< Power-of-two root rescaling exponent [nondim]
  integer :: r
    !< Residual exponent class [nondim]

  do n=2,8
    do r=0,n - 1
      base = scale(1.0 + (0.125 * real(n)), r)
      ref = nth_root(base, n)
      do k=-6,6
        scaled = scale(base, n * k)
        val = nth_root(scaled, n)
        call assert(val == scale(ref, k), "nth_root failed exact power-of-two rescaling")
      enddo
    enddo
  enddo
end subroutine test_nth_root_rescaling


!> Test nth_root special cases.
subroutine test_nth_root_special
  real :: val
    !< A computed nth root [nondim]
  real :: subnormal
    !< A subnormal positive input value [nondim]

  call assert(nth_root(0.0, 5) == 0.0, "nth_root(0,n) should be 0")
  call assert(nth_root(2.0, 1) == 2.0, "nth_root(x,1) should be x")

  subnormal = nearest(0.0, 1.0)
  val = nth_root(subnormal, 2)
  call assert(abs(val - sqrt(subnormal)) <= spacing(sqrt(subnormal)), &
      "nth_root of a positive subnormal value is inaccurate")

  val = nth_root(-1.0, 3)
  call assert(ieee_is_nan(val), "nth_root of a negative value should be NaN")

  val = nth_root(ieee_value(0., ieee_positive_inf), 3)
  call assert(val > huge(val), "nth_root of +Inf should be +Inf")

  val = nth_root(ieee_value(0., ieee_negative_inf), 3)
  call assert(ieee_is_nan(val), "nth_root of -Inf should be NaN")

  val = nth_root(1.0, 0)
  call assert(ieee_is_nan(val), "nth_root with n < 1 should be NaN")

  val = nth_root(1.0, 33)
  call assert(ieee_is_nan(val), "nth_root with n > 32 should be NaN")
end subroutine test_nth_root_special

end module MOM_nth_root_tests
