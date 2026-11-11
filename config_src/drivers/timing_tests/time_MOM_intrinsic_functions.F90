! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Timing tests for MOM_intrinsic_functions (exp_repro, log_repro, pow_repro)
program time_MOM_intrinsic_functions

use, intrinsic :: iso_fortran_env, only : int64, real64
use MOM_intrinsic_functions, only : exp_repro, log_repro, pow_repro

implicit none

integer, parameter :: npts = 100000
  !< Number of test points
integer, parameter :: niter = 200
  !< Number of timing iterations
real, parameter :: xmin = -10.
  !< Minimum x value [nondim]
real, parameter :: xmax = 10.
  !< Maximum x value [nondim]

real, allocatable :: x(:), x_log(:)
  ! Normal and log-distributed sample grids [nondim]
real, allocatable :: val(:)
  ! Function evaluations along sample grids [nondim]
integer(kind=int64) :: c1, c2
  ! Clock cycle measurements
integer(kind=int64) :: count_rate
  ! Clock cycles per second
real(kind=real64) :: t_cycle
  ! Nanoseconds per nominal clock cycle [wall-clock nanoseconds]
real(kind=real64) :: t_exp, t_exp_repro
  ! exp() evaluation time per loop iteration [wall-clock nanoseconds]
real(kind=real64) :: t_log, t_log_repro
  ! log() evaluation time per loop iteration [wall-clock nanoseconds]
real(kind=real64) :: t_scalar_baseline
  ! Scalar baseline time per iteration [wall-clock nanoseconds]
real(kind=real64) :: t_exp_scalar, t_exp_repro_scalar
  ! exp() evaluation time per scalar loop iteration [wall-clock nanoseconds]
real(kind=real64) :: t_log_scalar, t_log_repro_scalar
  ! exp() evaluation time per scalar loop iteration [wall-clock nanoseconds]
real :: I_npts
  ! Bit-reproducible 1/dx used in sample grid generation [nondim]
real, volatile :: scalar_x, scalar_val
  ! Unused summations, to preserve loops during optimization [nondim]
integer :: i, j

call system_clock(count_rate=count_rate)
t_cycle = 1e9 / real(count_rate, real64)

allocate(x(npts), x_log(npts), val(npts))
I_npts = 1. / (npts - 1)
do i = 1, npts
  x(i) = xmin + (i - 1) * ((xmax - xmin) * I_npts)
  x_log(i) = exp(x(i))
enddo

print '("=== MOM_intrinsic_functions timing ===")'
print '("npts = ", i0, ", niter = ", i0)', npts, niter
print '("x range: [", f6.1, ", ", f6.1, "]")', xmin, xmax
print '("log input range: [exp(", f6.1, "), exp(", f6.1, ")]")', xmin, xmax
print *

! Warm-up and time intrinsic exp()
do j = 1, 3
  val = exp(x)
  if (val(1) < 0.) val(1) = 0.
enddo
call system_clock(count=c1)
do j = 1, niter
  val = exp(x)
  if (val(1) < 0.) val(1) = 0.
enddo
call system_clock(count=c2)
t_exp = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("exp() time/elem:", t30, f8.2, " ns")', t_exp
print '("  sum (to prevent elision):", t30, ES12.5)', sum(val)

! Warm-up and time exp_repro()
do j = 1, 3
  val = exp_repro(x)
  if (val(1) < 0.) val(1) = 0.
enddo
call system_clock(count=c1)
do j = 1, niter
  val = exp_repro(x)
  if (val(1) < 0.) val(1) = 0.
enddo
call system_clock(count=c2)
t_exp_repro = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("exp_repro() time/elem:", t30, f8.2, " ns")', t_exp_repro
print '("  sum (to prevent elision):", t30, ES12.5)', sum(val)

print *
print '("slowdown factor:", t30, f8.2, "x")', t_exp_repro / t_exp

print *

! Warm-up and time intrinsic log()
do j = 1, 3
  val = log(x_log)
  if (val(1) > huge(val(1))) val(1) = 0.
enddo
call system_clock(count=c1)
do j = 1, niter
  val = log(x_log)
  if (val(1) > huge(val(1))) val(1) = 0.
enddo
call system_clock(count=c2)
t_log = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("log() time/elem:", t30, f8.2, " ns")', t_log
print '("  sum (to prevent elision):", t30, ES12.5)', sum(val)

! Warm-up and time log_repro()
do j = 1, 3
  val = log_repro(x_log)
  if (val(1) > huge(val(1))) val(1) = 0.
enddo
call system_clock(count=c1)
do j = 1, niter
  val = log_repro(x_log)
  if (val(1) > huge(val(1))) val(1) = 0.
enddo
call system_clock(count=c2)
t_log_repro = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("log_repro() time/elem:", t30, f8.2, " ns")', t_log_repro
print '("  sum (to prevent elision):", t30, ES12.5)', sum(val)

print *
print '("log slowdown factor:", t30, f8.2, "x")', t_log_repro / t_log

print *
print '("=== scalar loop-carried timing ===")'

! Warm-up and time a scalar dependency chain without exp().
! This gives an approximate loop overhead for interpreting the scalar timings.
do j = 1, 3
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = 1. + scalar_x
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo

call system_clock(count=c1)
do j = 1, niter
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = 1. + scalar_x
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo
call system_clock(count=c2)
t_scalar_baseline = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("baseline scalar time/call:", t30, f8.2, " ns")', t_scalar_baseline
print '("  final scalar value:", t30, ES12.5)', scalar_val

! Warm-up and time intrinsic exp() in a scalar loop with a carried dependency.
! This prevents the loop from being converted into independent vector lanes.
do j = 1, 3
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = exp(scalar_x)
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo

call system_clock(count=c1)
do j = 1, niter
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = exp(scalar_x)
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo
call system_clock(count=c2)
t_exp_scalar = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("exp() scalar time/call:", t30, f8.2, " ns")', t_exp_scalar
print '("  final scalar value:", t30, ES12.5)', scalar_val

! Warm-up and time exp_repro() in the same scalar loop pattern.
do j = 1, 3
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = exp_repro(scalar_x)
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo

call system_clock(count=c1)
do j = 1, niter
  scalar_x = 0.125
  do i = 1, npts
    scalar_val = exp_repro(scalar_x)
    scalar_x = scalar_x + (1.e-7 * (scalar_val - 1.))
    if (scalar_x > 0.5) scalar_x = scalar_x - 1.
  enddo
enddo
call system_clock(count=c2)
t_exp_repro_scalar = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("exp_repro() scalar time/call:", t30, f8.2, " ns")', t_exp_repro_scalar
print '("  final scalar value:", t30, ES12.5)', scalar_val

print *
print '("scalar slowdown factor:", t30, f8.2, "x")', &
    t_exp_repro_scalar / t_exp_scalar
print '("exp() minus baseline:", t30, f8.2, " ns")', &
    t_exp_scalar - t_scalar_baseline
print '("exp_repro() minus baseline:", t30, f8.2, " ns")', &
    t_exp_repro_scalar - t_scalar_baseline
print '("adjusted slowdown factor:", t30, f8.2, "x")', &
    (t_exp_repro_scalar - t_scalar_baseline) &
      / (t_exp_scalar - t_scalar_baseline)

print *
print '("=== scalar log loop-carried timing ===")'

! Warm-up and time intrinsic log() in a scalar loop with a carried dependency.
do j = 1, 3
  scalar_x = 1.125
  do i = 1, npts
    scalar_val = log(scalar_x)
    scalar_x = scalar_x + (1.e-7 * scalar_val)
    if (scalar_x > 1.5) scalar_x = scalar_x - 0.5
  enddo
enddo

call system_clock(count=c1)
do j = 1, niter
  scalar_x = 1.125
  do i = 1, npts
    scalar_val = log(scalar_x)
    scalar_x = scalar_x + (1.e-7 * scalar_val)
    if (scalar_x > 1.5) scalar_x = scalar_x - 0.5
  enddo
enddo
call system_clock(count=c2)
t_log_scalar = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("log() scalar time/call:", t30, f8.2, " ns")', t_log_scalar
print '("  final scalar value:", t30, ES12.5)', scalar_val

! Warm-up and time log_repro() in the same scalar loop pattern.
do j = 1, 3
  scalar_x = 1.125
  do i = 1, npts
    scalar_val = log_repro(scalar_x)
    scalar_x = scalar_x + (1.e-7 * scalar_val)
    if (scalar_x > 1.5) scalar_x = scalar_x - 0.5
  enddo
enddo

call system_clock(count=c1)
do j = 1, niter
  scalar_x = 1.125
  do i = 1, npts
    scalar_val = log_repro(scalar_x)
    scalar_x = scalar_x + (1.e-7 * scalar_val)
    if (scalar_x > 1.5) scalar_x = scalar_x - 0.5
  enddo
enddo
call system_clock(count=c2)
t_log_repro_scalar = real(c2 - c1, real64) * t_cycle / (niter * npts)

print '("log_repro() scalar time/call:", t30, f8.2, " ns")', t_log_repro_scalar
print '("  final scalar value:", t30, ES12.5)', scalar_val

print *
print '("log scalar slowdown factor:", t30, f8.2, "x")', &
    t_log_repro_scalar / t_log_scalar
print '("log() minus baseline:", t30, f8.2, " ns")', &
    t_log_scalar - t_scalar_baseline
print '("log_repro() minus baseline:", t30, f8.2, " ns")', &
    t_log_repro_scalar - t_scalar_baseline
print '("log adjusted slowdown factor:", t30, f8.2, "x")', &
    (t_log_repro_scalar - t_scalar_baseline) &
      / (t_log_scalar - t_scalar_baseline)

call time_power('general exponents', 0)
call time_power('eighth roots', 1)
call time_power('small integer powers', 2)

deallocate(x, x_log, val)

contains

!> Compare power implementations on general, root, and integral-power inputs.
subroutine time_power(label, category)
  character(len=*), intent(in) :: label !< Input category label
  integer, intent(in) :: category !< General (0), root (1), or integer (2)
  real :: powers(npts), results(npts)
    ! Exponents and results [nondim]
  real, volatile :: checksum, sx, sy, sb
    ! Observable checksum and scalar dependency-chain operands [nondim]
  real(kind=real64) :: elapsed
    ! Wall-clock nanoseconds per evaluation [wall-clock nanoseconds]
  integer(kind=int64) :: start, finish
    ! Clock counts
  integer :: method, iteration, k
  character(len=25), parameter :: names(4) = &
      [character(len=25) :: 'intrinsic **', 'pow_repro', 'exp_repro(b*log_repro(a))', 'sqrt(sqrt(sqrt(a)))']

  powers(:) = 0.37 + (0.01 * x(:))
  if (category == 1) powers(:) = 0.125
  if (category == 2) powers(:) = 3.
  print '(/,a,a)', '=== power timing: ', label
  do method=1,4
    if (method == 4 .and. category /= 1) cycle
    do iteration=-2,niter
      if (iteration == 1) call system_clock(count=start)
      select case (method)
      case (1)
        results(:) = x_log(:)**powers(:)
      case (2)
        results(:) = pow_repro(x_log(:), powers(:))
      case (3)
        results(:) = exp_repro(powers(:) * log_repro(x_log(:)))
      case (4)
        results(:) = sqrt(sqrt(sqrt(x_log(:))))
      end select
      checksum = results(1 + mod(abs(iteration), npts))
    enddo
    call system_clock(count=finish)
    elapsed = real(finish - start, real64) * t_cycle / (niter * npts)
    print '(a25,a,f8.2,a,es12.5)', names(method), ' array: ', elapsed, ' ns, check: ', checksum

    do iteration=-2,niter
      if (iteration == 1) call system_clock(count=start)
      sx = 1.125
      sb = powers(1)
      do k=1,npts
        select case (method)
        case (1)
          sy = sx**sb
        case (2)
          sy = pow_repro(sx, sb)
        case (3)
          sy = exp_repro(sb * log_repro(sx))
        case (4)
          sy = sqrt(sqrt(sqrt(sx)))
        end select
        sx = sx + (1.e-7 * (sy - 1.))
        if (sx > 1.5) sx = sx - 0.5
      enddo
    enddo
    call system_clock(count=finish)
    elapsed = real(finish - start, real64) * t_cycle / (niter * npts)
    print '(a25,a,f8.2,a,es12.5)', names(method), ' scalar:', elapsed, ' ns, check: ', sy
  enddo ! method loop
end subroutine time_power

end program time_MOM_intrinsic_functions
