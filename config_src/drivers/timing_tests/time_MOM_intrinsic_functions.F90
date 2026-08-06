! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Timing tests for MOM_intrinsic_functions (exp_repro)
program time_MOM_intrinsic_functions

use, intrinsic :: iso_fortran_env, only : int64, real64, real128
use MOM_intrinsic_functions, only : exp_repro

implicit none

integer, parameter :: npts = 10000000
  !< Number of test points
integer, parameter :: niter = 20
  !< Number of timing iterations
real, parameter :: xmin = -10.
  !< Minimum x value
real, parameter :: xmax = 10.
  !< Maximum x value

real, allocatable :: x(:), val(:)
integer(int64) :: count_rate, c1, c2
real(real64) :: clock_rate, time_intrinsic, time_repro
real :: I_npts
integer :: i, j

call system_clock(count_rate=count_rate)
clock_rate = real(count_rate, real64)

allocate(x(npts), val(npts))
I_npts = 1. / (npts - 1)
do i = 1, npts
  x(i) = xmin + (i - 1) * ((xmax - xmin) * I_npts)
enddo

print '("=== MOM_intrinsic_functions timing ===")'
print '("npts = ", i0, ", niter = ", i0)', npts, niter
print '("x range: [", f6.1, ", ", f6.1, "]")', xmin, xmax
print *

! Warm-up and time intrinsic exp()
do j = 1, 3
  val = exp(x)
enddo
call system_clock(count=c1)
do j = 1, niter
  val = exp(x)
enddo
call system_clock(count=c2)
time_intrinsic = real(c2 - c1, real64) / clock_rate / niter / npts * 1e9

print '("exp() time/elem:", t30, f8.2, " ns")', time_intrinsic
print '("  sum (to prevent elision):", t30, ES12.5)', sum(val)

! Warm-up and time exp_repro()
do j = 1, 3
  val = exp_repro(x)
enddo
call system_clock(count=c1)
do j = 1, niter
  val = exp_repro(x)
enddo
call system_clock(count=c2)
time_repro = real(c2 - c1, real64) / clock_rate / niter / npts * 1e9

print '("exp_repro() time/elem:", t30, f8.2, " ns")', time_repro
print '("  sum (to prevent elision):", t30, ES12.5)', sum(val)

print *
print '("slowdown factor:", t30, f8.2, "x")', time_repro / time_intrinsic

deallocate(x, val)

end program time_MOM_intrinsic_functions
