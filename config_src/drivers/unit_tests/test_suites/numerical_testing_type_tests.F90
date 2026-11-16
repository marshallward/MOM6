! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Tests for the numerical testing type.
module numerical_testing_type_tests

use, intrinsic :: iso_fortran_env, only : iostat_end
use numerical_testing_type, only : testing
use testdrive, only : check, error_type, new_unittest, unittest_type

implicit none ; private

public :: collect_tests

contains

!> Collect tests of passing, ignored, and failed checks.
subroutine collect_tests(tests)
  type(unittest_type), allocatable, intent(out) :: tests(:) !< Named tests

  tests = [ &
    new_unittest("test state", test_state), &
    new_unittest("real scalar", test_real_scalar), &
    new_unittest("real array", test_real_array), &
    new_unittest("integer array", test_integer_array), &
    new_unittest("summary", test_summary) ]
end subroutine collect_tests

!> Open a scratch stream for expected failure messages.
subroutine open_output(unit, error)
  integer, intent(out) :: unit !< Scratch stream
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  integer :: status ! I/O status

  open(newunit=unit, status="scratch", action="readwrite", iostat=status)
  call check(error, status, 0, "Opening scratch output")
end subroutine open_output

!> Configure a testing instance to use the scratch stream.
subroutine setup_test(instance, unit)
  type(testing), intent(inout) :: instance !< Instance to configure
  integer, intent(in) :: unit !< Scratch stream

  call instance%set(verbose=.false., stdout=unit, stderr=unit, stop_instantly=.false.)
end subroutine setup_test

!> Check and report the outcomes of passing, ignored, and failed checks.
subroutine check_outcomes(error, unit, passed, ignored, failed)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  integer, intent(in) :: unit !< Scratch stream
  type(testing), intent(inout) :: passed, ignored, failed !< Testing instances
  logical :: passed_state, ignored_state, failed_state ! Results of summaries
  integer :: status ! I/O status

  passed_state = passed%summarize("passing")
  ignored_state = ignored%summarize("ignored")
  failed_state = failed%summarize("failing")
  close(unit, iostat=status)
  call check(error, status, 0, "Closing scratch output")
  if (allocated(error)) return
  call check(error, .not. passed_state, "A passing check must not fail")
  if (allocated(error)) return
  call check(error, .not. ignored_state, "An ignored check must not fail")
  if (allocated(error)) return
  call check(error, failed_state, "A failing check must report failure")
end subroutine check_outcomes

!> Check that explicit test results update the state.
subroutine test_state(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(testing) :: passed, ignored, failed ! Independent test states
  integer :: unit ! Scratch stream

  call open_output(unit, error)
  if (allocated(error)) return
  call setup_test(passed, unit)
  call setup_test(ignored, unit)
  call setup_test(failed, unit)
  call passed%test(.false., "passing")
  call ignored%test(.true., "ignored", ignore=.true.)
  call failed%test(.true., "failing")
  call check_outcomes(error, unit, passed, ignored, failed)
end subroutine test_state

!> Check scalar real comparisons.
subroutine test_real_scalar(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(testing) :: passed, ignored, failed ! Independent test states
  integer :: unit ! Scratch stream

  call open_output(unit, error)
  if (allocated(error)) return
  call setup_test(passed, unit)
  call setup_test(ignored, unit)
  call setup_test(failed, unit)
  call passed%real_scalar(1., 1., "passing", robits=0, tol=0.)
  call ignored%real_scalar(1., 2., "ignored", ignore=.true.)
  call failed%real_scalar(1., 2., "failing")
  call check_outcomes(error, unit, passed, ignored, failed)
end subroutine test_real_scalar

!> Check real array comparisons.
subroutine test_real_array(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(testing) :: passed, ignored, failed ! Independent test states
  integer :: unit ! Scratch stream

  call open_output(unit, error)
  if (allocated(error)) return
  call setup_test(passed, unit)
  call setup_test(ignored, unit)
  call setup_test(failed, unit)
  call passed%real_arr(2, [1., 2.], [1., 2.], "passing", robits=0, tol=0.)
  call ignored%real_arr(2, [1., 2.], [3., 4.], "ignored", ignore=.true.)
  call failed%real_arr(2, [1., 2.], [3., 4.], "failing")
  call check_outcomes(error, unit, passed, ignored, failed)
end subroutine test_real_array

!> Check integer array comparisons.
subroutine test_integer_array(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(testing) :: passed, ignored, failed ! Independent test states
  integer :: unit ! Scratch stream

  call open_output(unit, error)
  if (allocated(error)) return
  call setup_test(passed, unit)
  call setup_test(ignored, unit)
  call setup_test(failed, unit)
  call passed%int_arr(2, [1, 2], [1, 2], "passing")
  call ignored%int_arr(2, [1, 2], [3, 4], "ignored", ignore=.true.)
  call failed%int_arr(2, [1, 2], [3, 4], "failing")
  call check_outcomes(error, unit, passed, ignored, failed)
end subroutine test_integer_array

!> Check summaries for an empty instance and a sequence of checks.
subroutine test_summary(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(testing) :: instance ! Test state
  logical :: empty_state, failed_state ! Results of summaries
  logical :: found_counts, found_first ! Whether the summary includes the expected details
  integer :: unit, status, close_status ! Scratch stream and I/O statuses
  character(len=256) :: line ! Output line

  call open_output(unit, error)
  if (allocated(error)) return
  call setup_test(instance, unit)
  empty_state = instance%summarize("empty")
  call instance%test(.false., "passing")
  call instance%test(.true., "ignored", ignore=.true.)
  call instance%test(.true., "first failure")
  call instance%real_scalar(1., 1., "passing", robits=0, tol=0.)
  call instance%real_scalar(1., 2., "ignored", ignore=.true.)
  call instance%real_scalar(1., 2., "failing")
  call instance%real_arr(2, [1., 2.], [1., 2.], "passing", robits=0, tol=0.)
  call instance%real_arr(2, [1., 2.], [3., 4.], "ignored", ignore=.true.)
  call instance%real_arr(2, [1., 2.], [3., 4.], "failing")
  call instance%int_arr(2, [1, 2], [1, 2], "passing")
  call instance%int_arr(2, [1, 2], [3, 4], "ignored", ignore=.true.)
  call instance%int_arr(2, [1, 2], [3, 4], "failing")
  failed_state = instance%summarize("failed")

  found_counts = .false.
  found_first = .false.
  rewind(unit, iostat=status)
  if (status == 0) then
    do
      read(unit, '(a)', iostat=status) line
      if (status /= 0) exit
      if (index(line, "4 failed of   12 tested") > 0) found_counts = .true.
      if (index(line, "First failed test: first failure") > 0) found_first = .true.
    enddo
    if (status == iostat_end) status = 0
  endif
  close(unit, iostat=close_status)
  call check(error, close_status, 0, "Closing scratch output")
  if (allocated(error)) return
  call check(error, status, 0, "Reading scratch output")
  if (allocated(error)) return
  call check(error, .not. empty_state, "An empty summary must pass")
  if (allocated(error)) return
  call check(error, failed_state, "A summary must report failure")
  if (allocated(error)) return
  call check(error, found_counts, "Summary must report four failures in twelve tests")
  if (allocated(error)) return
  call check(error, found_first, "Summary must report the first failure")
end subroutine test_summary

end module numerical_testing_type_tests
