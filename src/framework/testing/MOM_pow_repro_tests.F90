! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Tests for reproducible real powers, including accuracy against real128.
module MOM_pow_repro_tests

use, intrinsic :: iso_fortran_env, only : int64
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_is_nan, ieee_is_negative
use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_signaling_nan
use, intrinsic :: ieee_exceptions, only : ieee_set_flag, ieee_get_flag, ieee_support_flag
use, intrinsic :: ieee_exceptions, only : ieee_all, ieee_invalid, ieee_divide_by_zero, ieee_overflow, ieee_underflow
use MOM_error_handler, only : assert
use MOM_unit_testing, only : TestSuite
use MOM_intrinsic_functions, only : pow_repro

implicit none ; private
public :: add_pow_repro_tests

integer, parameter :: realquad = selected_real_kind(p=30, r=300)
  !< Quad precision if available
integer, parameter :: realq = merge(realquad, kind(1.), realquad > 0)
  !< Valid fallback kind for declarations

contains

!> Verify identities, small integer powers, and the elemental interface.
subroutine test_pow_exact
  real :: a(9), b(9), expected(9), p(9)
    ! Bases, exponents, exact references and results [nondim]

  a(:) = [0., -0., 2., -2., -2., -2., 4., 16., 256.]
  b(:) = [0., 0., 10., 3., 4., -3., 0.5, -0.25, 0.125]
  expected(:) = [1., 1., 1024., -8., 16., -0.125, 2., 0.5, 2.]
  p(:) = pow_repro(a(:), b(:))
  call assert(all(p == expected), 'pow_repro exact powers or elemental evaluation failed')
  call assert(pow_repro(huge(1.), 1.) == huge(1.), 'pow_repro identity at huge')
  call assert(pow_repro(tiny(1.), 1.) == tiny(1.), 'pow_repro identity at tiny')
  call assert(pow_repro(nearest(0., 1.), 1.) == nearest(0., 1.), 'pow_repro subnormal identity')
  call assert(pow_repro(-1., 2.**53) == 1., 'pow_repro large even exponent')
  call assert(pow_repro(-1., (2.**53) - 1.) == -1., 'pow_repro large odd exponent')
  call assert(pow_repro(-1., huge(1.)) == 1., 'pow_repro huge integral exponent')
  call assert(abs(pow_repro(-2., 17.) + 131072.) < 1.e-8, 'pow_repro general odd negative-base power')
  call assert(abs(pow_repro(-2., 18.) - 262144.) < 1.e-8, 'pow_repro general even negative-base power')
  call assert(pow_repro(2.**512, -2.) == 2.**(-1024), 'pow_repro negative power intermediate overflow')
end subroutine test_pow_exact

!> Check special values against IEEE-style power semantics, including signs.
subroutine test_pow_special
  real :: inf, nan, negzero, p
    ! Special inputs and result [nondim]
  integer :: n

  inf = ieee_value(0., ieee_positive_inf)
  nan = ieee_value(0., ieee_quiet_nan)
  negzero = -0.
  call assert(pow_repro(nan, 0.) == 1., 'NaN**0 should be one')
  call assert(pow_repro(1., nan) == 1., '1**NaN should be one')
  call assert(ieee_is_nan(pow_repro(nan, 2.)), 'NaN base should propagate')
  call assert(ieee_is_nan(pow_repro(2., nan)), 'NaN exponent should propagate')
  call assert(ieee_is_nan(pow_repro(-1., nan)), '(-1)**NaN should be NaN')
  call assert(ieee_is_nan(pow_repro(-2., 0.5)), 'Negative base and fractional exponent should be NaN')
  call assert(pow_repro(-1., inf) == 1., '(-1)**Inf should be one')
  call assert(pow_repro(-1., -inf) == 1., '(-1)**(-Inf) should be one')
  call assert(pow_repro(2., inf) == inf, '2**Inf should be Inf')
  call assert(pow_repro(0.5, inf) == 0., '0.5**Inf should be zero')
  call assert(pow_repro(2., -inf) == 0., '2**(-Inf) should be zero')
  call assert(pow_repro(0.5, -inf) == inf, '0.5**(-Inf) should be Inf')
  call assert(pow_repro(0., -inf) == inf, '0**(-Inf) should be Inf')
  call assert(pow_repro(inf, -inf) == 0., 'Inf**(-Inf) should be zero')
  call assert(pow_repro(inf, 0.) == 1., 'Inf**0 should be one')
  call assert(pow_repro(-inf, 0.5) == inf, '(-Inf)**0.5 should be Inf')
  do n=2,3
    p = pow_repro(negzero, real(n))
    call assert(p == 0., 'Negative zero to positive integer should be zero')
    call assert(ieee_is_negative(p) .eqv. (n == 3), 'Incorrect sign of zero power')
    p = pow_repro(negzero, -real(n))
    call assert(abs(p) == inf, 'Negative zero to negative integer should be infinite')
    call assert(ieee_is_negative(p) .eqv. (n == 3), 'Incorrect sign of infinite power')
    p = pow_repro(-inf, real(n))
    call assert(abs(p) == inf, 'Negative infinity to positive integer should be infinite')
    call assert(ieee_is_negative(p) .eqv. (n == 3), 'Incorrect sign of infinite-base power')
    p = pow_repro(-inf, -real(n))
    call assert(p == 0., 'Negative infinity to negative integer should be zero')
    call assert(ieee_is_negative(p) .eqv. (n == 3), 'Incorrect sign of infinite-base zero')
  enddo
  p = pow_repro(negzero, 0.5)
  call assert(p == 0. .and. .not. ieee_is_negative(p), '(-0)**0.5 should be positive zero')
end subroutine test_pow_special

!> Check exception flags, especially very large finite exponents.
subroutine test_pow_flags
  real, volatile :: a, b, p
    ! Inputs and result kept observable for exception tests [nondim]
  logical :: flag, other

  if (.not. ieee_support_flag(ieee_invalid)) return
  a = -2.
  b = 0.25
  call ieee_set_flag(ieee_all, .false.)
  p = pow_repro(a, b)
  call ieee_get_flag(ieee_invalid, flag)
  call assert(flag .and. ieee_is_nan(p), 'Fractional power of negative base must raise invalid')
  a = ieee_value(0., ieee_signaling_nan)
  b = 0.
  call ieee_set_flag(ieee_all, .false.)
  p = pow_repro(a, b)
  call ieee_get_flag(ieee_invalid, flag)
  call assert(flag .and. ieee_is_nan(p), 'Signaling NaN must raise invalid even with zero exponent')

  if (ieee_support_flag(ieee_divide_by_zero)) then
    a = -0.
    b = -3.
    call ieee_set_flag(ieee_all, .false.)
    p = pow_repro(a, b)
    call ieee_get_flag(ieee_divide_by_zero, flag)
    call assert(flag, 'Zero to negative power must raise divide by zero')
    b = 0.
    call ieee_set_flag(ieee_all, .false.)
    p = pow_repro(a, b)
    call ieee_get_flag(ieee_divide_by_zero, flag)
    call ieee_get_flag(ieee_invalid, other)
    call assert(p == 1. .and. .not. (flag .or. other), '0**0 must return one without invalid or divide by zero')
  endif
  if (ieee_support_flag(ieee_overflow) .and. ieee_support_flag(ieee_underflow)) then
    a = 0.125
    b = huge(1.)
    call ieee_set_flag(ieee_all, .false.)
    p = pow_repro(a, b)
    call ieee_get_flag(ieee_underflow, flag)
    call ieee_get_flag(ieee_overflow, other)
    call assert(p == 0. .and. flag .and. .not. other, 'Huge exponent should underflow without spurious overflow')
    a = 8.
    call ieee_set_flag(ieee_all, .false.)
    p = pow_repro(a, b)
    call ieee_get_flag(ieee_overflow, flag)
    call assert(flag .and. p > huge(1.), 'Huge exponent should overflow')
  endif
  call ieee_set_flag(ieee_all, .false.)
end subroutine test_pow_flags

!> Test subnormal inputs and results, and both ends of the finite range.
subroutine test_pow_range
  real :: bases(5), powers(8), p, expected, tolerance
    ! Inputs, computed and rounded reference results, absolute tolerance [nondim]
  real(kind=realq) :: ref
    ! Quad reference, including results outside the binary64 range [nondim]
  integer :: i, j

  bases(:) = [nearest(0., 1.), tiny(1.), nearest(tiny(1.), -1.), huge(1.), 2.]
  powers(:) = [-1.01, -0.5, -0.125, 0.37, 0.99, 1.01, -1074., -1075.]
  do j=1,size(powers)
    do i=1,size(bases)
      ! Restrict the very large exponents to base two (quad range is finite).
      if (j >= 7 .and. i /= 5) cycle
      ref = real(bases(i), realq)**real(powers(j), realq)
      p = pow_repro(bases(i), powers(j))
      if (ref > real(huge(1.), realq)) then
        call assert(p > huge(1.), 'pow_repro range overflow')
      else
        expected = real(ref)
        tolerance = max(2. * nearest(0., 1.), (2.e-13 * expected))
        call assert(abs(p - expected) <= tolerance, 'pow_repro range accuracy')
      endif
    enddo
  enddo
  call assert(pow_repro(2., -1074.) == nearest(0., 1.), 'pow_repro smallest subnormal result')
  call assert(pow_repro(2., -1076.) == 0., 'pow_repro below subnormal range')
end subroutine test_pow_range

!> Compare a wide dynamic range, near-one bases, and fast paths with quad power.
!! Report ULP error as well as relative error and enforce a near-roundoff bound.
subroutine test_pow_accuracy
  integer, parameter :: npts = 12000
  real :: a(npts), b(npts), scalar(npts), vector(npts)
    ! Input pairs and scalar/elemental results [nondim]
  real :: relerr, max_relerr, max_ulp, tolerance
    ! Relative error, maxima, and ULP tolerance [nondim]
  real(kind=realq) :: ref, err
    ! Quad reference and absolute error [nondim]
  integer :: i, j

  do i=1,npts
    j = mod(i, 2001) - 1000
    a(i) = scale(1. + real(mod(i, 97)) / 128., j)
    b(i) = real(mod(i, 199) - 99) / 111.
    if (mod(i, 3) == 0) then
      a(i) = 1. + real(j) * epsilon(1.)
      b(i) = real(mod(i, 199) - 99) * 1.e10
    endif
    if (mod(i, 5) == 0) then
      a(i) = 0.5 + real(mod(i, 97)) / 64.
      b(i) = real(mod(i, 33) - 16)
    endif
    if (mod(i, 7) == 0) b(i) = sign(0.125 * (2.**mod(i, 3)), b(i))
  enddo
  do i=1,npts
    scalar(i) = pow_repro(a(i), b(i))
    if (scalar(i) < 0.) exit ! Prevent vectorization of this loop.
  enddo
  vector(:) = pow_repro(a(:), b(:))
  call assert(all(transfer(scalar, [0_int64], npts) == transfer(vector, [0_int64], npts)), &
      'pow_repro scalar and elemental results differ')
  max_relerr = 0.
  max_ulp = 0.
  do i=1,npts
    ref = real(a(i), realq)**real(b(i), realq)
    err = abs(real(scalar(i), realq) - ref)
    relerr = real(err / abs(ref))
    max_relerr = max(max_relerr, relerr)
    max_ulp = max(max_ulp, real(err / real(spacing(real(ref)), realq)))
    tolerance = 3. * spacing(real(ref))
    call assert(err <= real(tolerance, realq), 'pow_repro error exceeds 3 ULP')
  enddo
  print '(a,es12.4,a,f10.2)', ' pow_repro max relative error: ', max_relerr, ', max ULP error: ', max_ulp
end subroutine test_pow_accuracy

!> Compare scalar and array power accuracy with a quad-precision reference.
subroutine test_pow_accuracy_comparison
  integer, parameter :: npts = 100000
  real, parameter :: amin = 1.e-10, amax = 1.e10, bmin = -10., bmax = 10.
    ! Base and exponent interval endpoints [nondim]
  real :: a(npts), b(npts), val(npts), val_intrinsic(npts), val_vec(npts), intrinsic_vec(npts)
    ! Inputs and scalar/array results from both implementations [nondim]
  real(kind=realq) :: ref(npts)
    ! Quad-precision reference powers [nondim]
  integer(kind=int64) :: first_bits, span, step, remainder, offset, denom, bits
    ! Binary64 input grid coordinates and integer interpolation terms
  integer :: i

  ! Sample base bit coordinates uniformly for broad scale-aware coverage,
  ! without using either implementation to generate its own test inputs.
  first_bits = transfer(amin, 0_int64)
  span = transfer(amax, 0_int64) - first_bits
  denom = int(npts - 1, int64)
  step = span / denom
  remainder = mod(span, denom)
  do i=1,npts
    offset = int(i - 1, int64)
    bits = first_bits + (offset * step + (offset * remainder) / denom)
    a(i) = transfer(bits, 0.)
    ! A coprime permutation spreads exponents across the base interval rather
    ! than correlating large bases with large exponents. Results stay normal.
    b(i) = bmin + real(mod((i - 1) * 7919, npts)) * ((bmax - bmin) / real(npts - 1))
  enddo

  do i=1,npts
    val_intrinsic(i) = a(i)**b(i)
    val(i) = pow_repro(a(i), b(i))
    ref(i) = real(a(i), realq)**real(b(i), realq)
    if (val(i) < 0.) exit ! Impossible branch to prevent vectorization.
  enddo
  val_vec(:) = pow_repro(a(:), b(:))
  intrinsic_vec(:) = a(:)**b(:)
  call assert(all(transfer(val, [0_int64], npts) == transfer(val_vec, [0_int64], npts)), &
      'pow_repro scalar and array accuracy samples differ bitwise')

  print '(1x,a)', '=== scalar pow_repro() and intrinsic ** accuracy'
  print '(2x,"Tested ",i0," pairs across a interval [",ES12.5,",",ES12.5,"] and b interval [",F6.1,",",F6.1,"]")', &
      npts, amin, amax, bmin, bmax
  call report_pow_accuracy('pow_repro', a, b, val, ref, check_bound=.true.)
  call report_pow_accuracy('intrinsic **', a, b, val_intrinsic, ref)
  print '(1x,a)', '=== vector pow_repro() and intrinsic ** accuracy'
  print '(2x,"Tested ",i0," pairs across a interval [",ES12.5,",",ES12.5,"] and b interval [",F6.1,",",F6.1,"]")', &
      npts, amin, amax, bmin, bmax
  call report_pow_accuracy('pow_repro', a, b, val_vec, ref, check_bound=.true.)
  call report_pow_accuracy('intrinsic **', a, b, intrinsic_vec, ref)
end subroutine test_pow_accuracy_comparison

!> Stress large exponents near one and both sides of logarithm table boundaries.
subroutine test_pow_compensation
  integer, parameter :: npts = 6000
  real :: a, b, p, max_ulp
    ! Base, exponent, result and maximum ULP error [nondim]
  real(kind=realq) :: ref, ulp_error
    ! Quad reference and error in binary64 ULPs [nondim]
  integer :: i, j

  max_ulp = 0.
  do i=1,npts
    j = mod(i, 1000) + 1
    select case (mod(i, 3))
    case (0)
      a = 1. + real(j) * epsilon(1.)
    case (1)
      a = 1. - real(j) * epsilon(1.)
    case (2)
      a = 1. + real(j) / 262144.
    end select
    b = real(mod(i * 17, 1361) - 680) / (a - 1.)
    p = pow_repro(a, b)
    ref = real(a, realq)**real(b, realq)
    ulp_error = abs(real(p, realq) - ref) / real(spacing(real(ref)), realq)
    max_ulp = max(max_ulp, real(ulp_error))
    call assert(ulp_error <= 2._realq, 'pow_repro near-one compensation exceeds 2 ULP')
  enddo
  do i=0,128
    do j=-1,1,2
      a = nearest(1. + real(i) / 128., real(j))
      b = 999.25
      p = pow_repro(a, b)
      ref = real(a, realq)**real(b, realq)
      ulp_error = abs(real(p, realq) - ref) / real(spacing(real(ref)), realq)
      max_ulp = max(max_ulp, real(ulp_error))
      call assert(ulp_error <= 2._realq, 'pow_repro table boundary exceeds 2 ULP')
    enddo
  enddo
  print '(2x,a,f12.8)', 'pow_repro near-one/table-boundary max ULP error: ', max_ulp
end subroutine test_pow_compensation

!> Report power errors and both inputs at each maximum using a quad reference.
!! Intrinsic accuracy is informational; the reproducible implementation must
!! remain within 2 ULP on this grid of normal results.
subroutine report_pow_accuracy(label, a, b, val, ref, check_bound)
  character(len=*), intent(in) :: label !< Implementation name
  real, intent(in) :: a(:), b(:), val(:)
    !< Bases, exponents, and computed powers [nondim]
  real(kind=realq), intent(in) :: ref(:)
    !< Quad reference powers, all positive and normal in default real [nondim]
  logical, optional, intent(in) :: check_bound !< Enforce the reproducible ULP error bound
  real(kind=realq) :: err(3), max_err(3)
    ! Absolute, relative and ULP errors and maxima [nondim]
  integer :: max_index(3), count_half_ulp, i, j
    ! Indices of maxima, count above half an ULP, and loop indices

  max_err(:) = 0.
  max_index(:) = 1
  count_half_ulp = 0
  do i=1,size(a)
    err(1) = abs(real(val(i), realq) - ref(i))
    err(2) = err(1) / abs(ref(i))
    err(3) = err(1) / real(spacing(real(ref(i))), realq)
    do j=1,3
      if (err(j) > max_err(j)) then
        max_err(j) = err(j)
        max_index(j) = i
      endif
    enddo
    if (err(3) > 0.5_realq) count_half_ulp = count_half_ulp + 1
    if (present(check_bound)) then
      if (check_bound) then
        call assert(err(3) <= 2._realq, 'pow_repro comparison exceeds 2 ULP')
      endif
    endif
  enddo ! sample loop

  print '(2x,a," max abs err:",t36,ES12.5," at (a, b) = (",ES14.5E3,",",ES14.5E3,")")', &
      label, max_err(1), a(max_index(1)), b(max_index(1))
  print '(2x,a," max rel err:",t36,ES12.5," at (a, b) = (",ES14.5E3,",",ES14.5E3,")")', &
      label, max_err(2), a(max_index(2)), b(max_index(2))
  print '(2x,a," max ULP err:",t36,F16.10," at (a, b) = (",ES14.5E3,",",ES14.5E3,")")', &
      label, max_err(3), a(max_index(3)), b(max_index(3))
  print '(2x,a," points > 0.5 ULP:",t36,i12)', label, count_half_ulp
end subroutine report_pow_accuracy

!> Register power tests in the intrinsic suite.
subroutine add_pow_repro_tests(suite)
  type(TestSuite), intent(inout) :: suite !< Test suite to extend

  call suite%add(test_pow_exact, 'test_pow_exact')
  call suite%add(test_pow_special, 'test_pow_special')
  call suite%add(test_pow_flags, 'test_pow_flags')
  if (realquad > 0) then
    call suite%add(test_pow_accuracy, 'test_pow_accuracy')
    call suite%add(test_pow_accuracy_comparison, 'test_pow_accuracy_comparison')
    call suite%add(test_pow_compensation, 'test_pow_compensation')
    call suite%add(test_pow_range, 'test_pow_range')
  endif
end subroutine add_pow_repro_tests

end module MOM_pow_repro_tests
