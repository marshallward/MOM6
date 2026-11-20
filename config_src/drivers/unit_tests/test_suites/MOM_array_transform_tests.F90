! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Test symmetric sums independently of model startup.
module MOM_array_transform_tests

use testdrive, only : check, error_type, new_unittest, test_failed, unittest_type
use MOM_array_transform, only : rotate_array, symmetric_sum

implicit none ; private

public :: collect_tests

integer, parameter :: max_size = 13 ! Maximum array dimension

contains

!> Collect the symmetric-sum tests.
subroutine collect_tests(tests)
  type(unittest_type), allocatable, intent(out) :: tests(:) !< Named tests

  tests = [ &
    new_unittest("sum compared with naive sum", test_sum_accuracy), &
    new_unittest("90-degree rotation", test_rotation_90), &
    new_unittest("180-degree rotation", test_rotation_180), &
    new_unittest("270-degree rotation", test_rotation_270) ]
end subroutine collect_tests

!> Fill an array with alternating positive and negative inexact values.
subroutine fill_array(array)
  real, intent(out) :: array(max_size,max_size) !< Test values [A]
  integer :: i, j ! Array indices

  do j = 1, max_size ; do i = 1, max_size
    array(i,j) = 1.0 / (2.0*(j*max_size + i) + 1.0)
    if (modulo(i+j, 2) == 0) array(i,j) = -array(i,j)
  enddo ; enddo
end subroutine fill_array

!> Compare every subarray's symmetric sum against a naive sum.
subroutine test_sum_accuracy(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  real :: array(max_size,max_size) ! Test values [A]
  real :: actual, naive, abs_sum, tol ! Sums and their tolerance [A]
  integer :: i, j, n, m ! Array indices and subarray dimensions
  character(len=160) :: message ! Failure context

  call fill_array(array)
  do n = 1, max_size ; do m = 1, max_size
    actual = symmetric_sum(array(1:n,1:m))
    naive = 0.0
    abs_sum = 0.0
    do j = 1, m ; do i = 1, n
      naive = naive + array(i,j)
      abs_sum = abs_sum + abs(array(i,j))
    enddo ; enddo
    tol = (2.0 * abs_sum) * epsilon(abs_sum)
    write(message,'(i0," x ",i0," symmetric vs naive sum, sum=",es13.5," diff=",es13.5)') &
        n, m, actual, naive - actual
    call check(error, abs(actual - naive) <= tol, trim(message))
    if (allocated(error)) return
  enddo ; enddo
end subroutine test_sum_accuracy

!> Check invariance under a 90-degree rotation.
subroutine test_rotation_90(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any

  call check_rotation(error, 1)
end subroutine test_rotation_90

!> Check invariance under a 180-degree rotation.
subroutine test_rotation_180(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any

  call check_rotation(error, 2)
end subroutine test_rotation_180

!> Check invariance under a 270-degree rotation.
subroutine test_rotation_270(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any

  call check_rotation(error, 3)
end subroutine test_rotation_270

!> Compare every subarray's symmetric sum with its rotated counterpart.
subroutine check_rotation(error, turns)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  integer, intent(in) :: turns !< Number of quarter turns
  real :: array(max_size,max_size), rotated(max_size,max_size) ! Test values [A]
  real :: actual, rotated_sum ! Original and rotated sums [A]
  integer :: n, m ! Subarray dimensions
  character(len=160) :: message ! Failure context

  call fill_array(array)
  call rotate_array(array, turns, rotated)
  do n = 1, max_size ; do m = 1, max_size
    actual = symmetric_sum(array(1:n,1:m))
    select case (turns)
      case (1)
        rotated_sum = symmetric_sum(rotated(max_size+1-m:max_size,1:n))
      case (2)
        rotated_sum = symmetric_sum(rotated(max_size+1-n:max_size,max_size+1-m:max_size))
      case (3)
        rotated_sum = symmetric_sum(rotated(1:m,max_size+1-n:max_size))
      case default
        call test_failed(error, "Invalid number of quarter turns")
        return
    end select
    write(message,'(i0," x ",i0," with ",i0," degree rotation, sum=",es13.5," diff=",es13.5)') &
        n, m, 90*turns, actual, rotated_sum - actual
    call check(error, actual == rotated_sum, trim(message))
    if (allocated(error)) return
  enddo ; enddo
end subroutine check_rotation

end module MOM_array_transform_tests
