! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Focused profiling driver for nth_root().
program time_nth_root_perf

use MOM_intrinsic_functions, only : nth_root

implicit none

integer, parameter :: default_degree = 5
  !< Default degree of the root to profile [nondim]
integer, parameter :: default_niter = 500
  !< Default number of times to repeat the profiling loop [nondim]
integer, parameter :: default_npts = 262144
  !< Default number of values in the input array [nondim]
real, parameter :: xmin = 1.e-4
  !< Minimum input value for the array profiling loop [nondim]
real, parameter :: xmax = 10.
  !< Maximum input value for the array profiling loop [nondim]

character(len=32) :: mode
  !< Profiling mode, either "array" or "scalar"
real, allocatable :: x(:)
  !< Input values for array profiling [nondim]
real, allocatable :: val(:)
  !< Output values for array profiling [nondim]
real :: checksum
  !< Accumulated value used to prevent optimization elision [nondim]
real :: I_npts
  !< Inverse number of intervals in the input grid [nondim]
integer :: degree
  !< Degree of the root to profile [nondim]
integer :: niter
  !< Number of times to repeat the profiling loop [nondim]
integer :: npts
  !< Number of values in the input array [nondim]
integer :: i
  !< Loop counter [nondim]

degree = default_degree
niter = default_niter
npts = default_npts
mode = 'scalar'

call read_integer_arg(1, degree)
call read_integer_arg(2, niter)
call read_integer_arg(3, npts)
call read_mode_arg(4, mode)

if ((degree < 1) .or. (degree > 32)) error stop 'degree must be in [1,32]'
if (niter < 1) error stop 'niter must be positive'
if (npts < 2) error stop 'npts must be at least 2'
if ((trim(mode) /= 'array') .and. (trim(mode) /= 'scalar')) &
  error stop 'mode must be array or scalar'

allocate(x(npts), val(npts))

I_npts = 1. / real(npts - 1)
do i = 1,npts
  x(i) = xmin + real(i - 1) * ((xmax - xmin) * I_npts)
  val(i) = 0.
enddo

if (trim(mode) == 'array') then
  call run_array_profile(x, val, degree, 3)
  call run_array_profile(x, val, degree, niter)
  checksum = array_checksum(val)
else
  call run_scalar_profile(degree, npts, 3, checksum)
  call run_scalar_profile(degree, npts, niter, checksum)
endif

print '("mode: ", a)', trim(mode)
print '("degree: ", i0)', degree
print '("npts: ", i0)', npts
print '("niter: ", i0)', niter
print '("checksum: ", ES22.15)', checksum

deallocate(x, val)

contains

!> Read an integer command-line argument if it is present.
subroutine read_integer_arg(arg_index, value)
  integer, intent(in) :: arg_index
    !< Command-line argument index [nondim]
  integer, intent(inout) :: value
    !< Parsed integer value [nondim]

  character(len=32) :: arg
    !< Command-line argument text [nondim]
  integer :: status
    !< Argument-read status [nondim]

  call get_command_argument(arg_index, arg, status=status)
  if (status == 0) read(arg, *) value
end subroutine read_integer_arg


!> Read the profiling mode command-line argument if it is present.
subroutine read_mode_arg(arg_index, value)
  integer, intent(in) :: arg_index
    !< Command-line argument index [nondim]
  character(len=*), intent(inout) :: value
    !< Parsed profiling mode [nondim]

  character(len=32) :: arg
    !< Command-line argument text [nondim]
  integer :: status
    !< Argument-read status [nondim]

  call get_command_argument(arg_index, arg, status=status)
  if (status == 0) value = trim(adjustl(arg))
end subroutine read_mode_arg


!> Profile repeated nth_root() calls over an array of inputs.
subroutine run_array_profile(x, val, degree, niter)
  real, intent(in) :: x(:)
    !< Input values for array profiling [nondim]
  real, intent(out) :: val(:)
    !< Output values for array profiling [nondim]
  integer, intent(in) :: degree
    !< Degree of the root to profile [nondim]
  integer, intent(in) :: niter
    !< Number of times to repeat the profiling loop [nondim]

  integer :: i
    !< Array index [nondim]
  integer :: itt
    !< Iteration counter [nondim]

  do itt = 1,niter
    do i = 1,size(x)
      val(i) = nth_root(x(i), degree)
    enddo
  enddo
end subroutine run_array_profile


!> Profile repeated nth_root() calls in a scalar dependency chain.
subroutine run_scalar_profile(degree, npts, niter, result)
  integer, intent(in) :: degree
    !< Degree of the root to profile [nondim]
  integer, intent(in) :: npts
    !< Number of values in each profiling loop [nondim]
  integer, intent(in) :: niter
    !< Number of times to repeat the profiling loop [nondim]
  real, intent(out) :: result
    !< Final scalar result used to prevent optimization elision [nondim]

  real, volatile :: scalar_x
    !< Loop-carried scalar input value [nondim]
  real, volatile :: scalar_val
    !< Computed scalar nth root [nondim]
  integer :: i
    !< Inner loop counter [nondim]
  integer :: itt
    !< Iteration counter [nondim]

  scalar_x = 0.125
  scalar_val = 0.
  do itt = 1,niter
    do i = 1,npts
      scalar_val = nth_root(scalar_x, degree)
      scalar_x = scalar_x + (1.e-7 * scalar_val)
      if (scalar_x > 0.5) scalar_x = scalar_x - 0.375
    enddo
  enddo
  result = scalar_val
end subroutine run_scalar_profile


!> Compute a checksum with explicit summation order.
function array_checksum(val) result(checksum)
  real, intent(in) :: val(:)
    !< Values to accumulate [nondim]
  real :: checksum
    !< Accumulated values [nondim]

  integer :: i
    !< Array index [nondim]

  checksum = 0.
  do i = 1,size(val)
    checksum = checksum + val(i)
  enddo
end function array_checksum

end program time_nth_root_perf
