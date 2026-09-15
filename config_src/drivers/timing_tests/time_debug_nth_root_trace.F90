! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Standalone diagnostic trace for nth_root() arithmetic.
program time_debug_nth_root_trace

use, intrinsic :: iso_fortran_env, only : int64
use, intrinsic :: ieee_arithmetic, only : ieee_fma
use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value

implicit none

real, parameter :: real_mold = 0.
  !< Real mold for numerical format queries [nondim]
integer, parameter :: expbit = digits(real_mold) - 1
  !< Position of lowest exponent bit [nondim]
integer, parameter :: signbit = storage_size(real_mold) - 1
  !< Position of sign bit [nondim]
integer, parameter :: expwidth = signbit - expbit
  !< Number of exponent bits [nondim]
integer, parameter :: expbias = maxexponent(real_mold) - 1
  !< Exponent bias [nondim]

real :: x
  !< Input value for nth-root tracing [A^n]
real :: root
  !< Traced nth root of x [A]
integer :: n
  !< Degree of the root [nondim]
integer :: forced_iterations
  !< Forced number of Halley iterations, or -1 for the table [nondim]
logical :: do_newton
  !< If true, apply the compensated Newton correction [nondim]
logical :: do_select_best
  !< If true, test the adjacent floating-point root [nondim]

call read_arguments(x, n, forced_iterations, do_newton, do_select_best)

write(*,'(a)') '=== debug_nth_root_trace ==='
call trace_real('input x', x)
call trace_int('input n', int(n, int64))
call trace_int('forced_iterations', int(forced_iterations, int64))
write(*,'(a,l1)') 'do_newton      = ', do_newton
write(*,'(a,l1)') 'do_select_best = ', do_select_best

root = nth_root_trace(x, n, forced_iterations, do_newton, do_select_best)
call trace_real('final root', root)

contains

!> Read command-line arguments for the trace program.
subroutine read_arguments(x, n, forced_iterations, do_newton, do_select_best)
  real, intent(out) :: x
    !< Input value for nth-root tracing [A^n]
  integer, intent(out) :: n
    !< Degree of the root [nondim]
  integer, intent(out) :: forced_iterations
    !< Forced number of Halley iterations, or -1 for the table [nondim]
  logical, intent(out) :: do_newton
    !< If true, apply the compensated Newton correction [nondim]
  logical, intent(out) :: do_select_best
    !< If true, test the adjacent floating-point root [nondim]

  character(len=128) :: arg
    !< Command-line argument text [nondim]
  integer :: arg_count
    !< Number of command-line arguments [nondim]
  integer :: flag
    !< Integer flag used for logical arguments [nondim]

  x = 0.802209
  n = 6
  forced_iterations = -1
  do_newton = .true.
  do_select_best = .true.

  arg_count = command_argument_count()

  if (arg_count >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) x
  endif
  if (arg_count >= 2) then
    call get_command_argument(2, arg)
    read(arg, *) n
  endif
  if (arg_count >= 3) then
    call get_command_argument(3, arg)
    read(arg, *) forced_iterations
  endif
  if (arg_count >= 4) then
    call get_command_argument(4, arg)
    read(arg, *) flag
    do_newton = (flag /= 0)
  endif
  if (arg_count >= 5) then
    call get_command_argument(5, arg)
    read(arg, *) flag
    do_select_best = (flag /= 0)
  endif
end subroutine read_arguments


!> Trace a real value in decimal and as transferred integer bits.
subroutine trace_real(label, value)
  character(len=*), intent(in) :: label
    !< Label to print before the traced value [nondim]
  real, intent(in) :: value
    !< Real value to print [nondim]

  integer(kind=int64) :: bits
    !< Bit representation of value [nondim]

  bits = transfer(value, 0_int64)
  write(*,'(a32,1x,es27.18e3,1x,z16.16)') label, value, bits
end subroutine trace_real


!> Trace an integer value in decimal and hexadecimal.
subroutine trace_int(label, value)
  character(len=*), intent(in) :: label
    !< Label to print before the traced value [nondim]
  integer(kind=int64), intent(in) :: value
    !< Integer value to print [nondim]

  write(*,'(a32,1x,i0,1x,z16.16)') label, value, value
end subroutine trace_int


!> Trace the nth-root calculation.
function nth_root_trace(x, n, forced_iterations, do_newton, do_select_best) result(root)
  integer, parameter :: max_nth_root = 32
    !< Largest root degree supported by nth_root [nondim]
  integer, parameter :: halley_iterations_by_n(max_nth_root) = [ &
      0, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, &
      5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8]
    !< Number of Halley iterations used by each supported root degree [nondim]
  real, intent(in) :: x
    !< Input value for nth-root tracing [A^n]
  integer, intent(in) :: n
    !< Degree of the root [nondim]
  integer, intent(in) :: forced_iterations
    !< Forced number of Halley iterations, or -1 for the table [nondim]
  logical, intent(in) :: do_newton
    !< If true, apply the compensated Newton correction [nondim]
  logical, intent(in) :: do_select_best
    !< If true, test the adjacent floating-point root [nondim]
  real :: root
    !< Traced nth root of x [A]

  real :: xr
    !< The rescaled value of x in the range [2**(-n),1) [B^n]
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
  real :: ratio
    !< Halley numerator-to-denominator ratio [nondim]
  real :: resid
    !< Compensated residual for the Newton correction [B^n]
  real :: correction
    !< Newton correction to root_xr [B]
  real :: rn_minus
    !< n-1 represented as a real number [nondim]
  real :: rn_plus
    !< n+1 represented as a real number [nondim]
  integer(kind=int64) :: e_x
    !< Integral component of the nth-root exponent of x [nondim]
  integer :: halley_iterations
    !< Number of Halley iterations to use [nondim]
  integer :: itt
    !< Iteration counter [nondim]

  if ((n < 1) .or. (n > max_nth_root)) then
    root = ieee_value(root, ieee_quiet_nan)
    call trace_real('root invalid n', root)
  elseif (x < 0.0) then
    root = ieee_value(root, ieee_quiet_nan)
    call trace_real('root negative x', root)
  elseif (is_nonfinite_trace(x)) then
    root = x + x
    call trace_real('root nonfinite x+x', root)
  elseif (x == 0.0) then
    root = x
    call trace_real('root zero x', root)
  elseif (n == 1) then
    root = x
    call trace_real('root n==1 x', root)
  else
    call rescale_nth_root_trace(x, n, xr, e_x)
    call trace_real('after rescale xr', xr)
    call trace_int('after rescale e_x', e_x)

    root_xr = 0.7071067811865475
    call trace_real('root_xr initial', root_xr)
    rn_minus = real(n - 1)
    call trace_real('rn_minus', rn_minus)
    rn_plus = real(n + 1)
    call trace_real('rn_plus', rn_plus)

    halley_iterations = halley_iterations_by_n(n)
    if (forced_iterations >= 0) halley_iterations = forced_iterations
    call trace_int('halley_iterations', int(halley_iterations, int64))

    do itt=1,halley_iterations
      write(*,'(/,a,i0)') '--- Halley iteration ', itt
      root_n = integer_power_trace(root_xr, n, 'halley root_n')
      call trace_real('halley root_n assigned', root_n)
      num = (rn_minus * root_n) + (rn_plus * xr)
      call trace_real('halley num', num)
      den = (rn_plus * root_n) + (rn_minus * xr)
      call trace_real('halley den', den)
      ratio = num / den
      call trace_real('halley ratio', ratio)
      root_xr = root_xr * ratio
      call trace_real('halley root_xr', root_xr)
    enddo

    if (do_newton) then
      write(*,'(/,a)') '--- Newton polish'
      root_nm1 = integer_power_trace(root_xr, n - 1, 'newton root_nm1')
      call trace_real('newton root_nm1 assigned', root_nm1)
      resid = power_residual_trace(xr, root_xr, n, 'newton residual')
      call trace_real('newton residual assigned', resid)
      correction = resid / (real(n) * root_nm1)
      call trace_real('newton correction', correction)
      root_xr = root_xr + correction
      call trace_real('newton root_xr', root_xr)
    endif

    if (do_select_best) then
      write(*,'(/,a)') '--- Select best root'
      root_xr = select_best_root_trace(root_xr, xr, n)
      call trace_real('selected root_xr', root_xr)
    endif

    write(*,'(/,a)') '--- Descale'
    root = descale_nth_root_trace(root_xr, e_x)
    call trace_real('descaled root', root)
  endif
end function nth_root_trace


!> Rescale x to the range [2**(-n),1) and trace the exponent operations.
subroutine rescale_nth_root_trace(x, n, r, e_a)
  real, intent(in) :: x
    !< The number to be rescaled for nth-root computation [A^n]
  integer, intent(in) :: n
    !< The degree of the root [nondim]
  real, intent(out) :: r
    !< The rescaled value of x in the range [2**(-n),1) [B^n]
  integer(kind=int64), intent(out) :: e_a
    !< The integral component of the nth-root exponent of x [nondim]

  integer(kind=int64) :: xb
    !< Floating point integer representation of x [nondim]
  integer(kind=int64) :: exp_bits
    !< Raw exponent bits from xb [nondim]
  integer(kind=int64) :: e_x
    !< Exponent of input x [nondim]
  integer(kind=int64) :: e_r
    !< Residual exponent of rescaled r [nondim]
  integer(kind=int64) :: e_shift
    !< Normalizing exponent shift applied to subnormal inputs [nondim]
  integer(kind=int64) :: n64
    !< The root degree promoted to the exponent integer kind [nondim]

  write(*,'(/,a)') '--- Rescale'
  xb = transfer(x, 1_int64)
  call trace_int('rescale xb transfer', xb)
  e_shift = 0_int64
  call trace_int('rescale e_shift init', e_shift)
  exp_bits = ibits(xb, expbit, expwidth)
  call trace_int('rescale exp_bits', exp_bits)

  if (exp_bits == 0_int64) then
    e_shift = int(digits(real_mold), int64)
    call trace_int('rescale e_shift subnormal', e_shift)
    xb = transfer(scale(x, int(e_shift)), 1_int64)
    call trace_int('rescale xb scaled', xb)
    exp_bits = ibits(xb, expbit, expwidth)
    call trace_int('rescale exp_bits scaled', exp_bits)
  endif

  e_x = ibits(xb, expbit, expwidth) - expbias
  call trace_int('rescale e_x biased', e_x)
  e_x = e_x - e_shift
  call trace_int('rescale e_x shifted', e_x)
  n64 = int(n, int64)
  call trace_int('rescale n64', n64)

  if (e_x >= 0_int64) then
    e_a = (e_x + n64) / n64
    call trace_int('rescale e_a nonnegative', e_a)
  else
    e_a = (e_x + 1_int64) / n64
    call trace_int('rescale e_a negative', e_a)
  endif
  e_r = e_x - e_a * n64
  call trace_int('rescale e_r residual', e_r)

  call mvbits(e_r + expbias, 0, expwidth + 1, xb, expbit)
  call trace_int('rescale xb mvbits', xb)
  r = transfer(xb, 1.)
  call trace_real('rescale r transfer', r)
end subroutine rescale_nth_root_trace


!> Undo the rescaling and trace the bit operations.
function descale_nth_root_trace(x, e_a) result(a)
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
  call trace_int('descale xb transfer', xb)
  e_x = ibits(xb, expbit, expwidth)
  call trace_int('descale e_x biased', e_x)
  call trace_int('descale e_a', e_a)
  call mvbits(e_a + e_x, 0, expwidth, xb, expbit)
  call trace_int('descale xb mvbits', xb)
  a = transfer(xb, 1.)
  call trace_real('descale a transfer', a)
end function descale_nth_root_trace


!> Return true if a real value is an IEEE Inf or NaN value and trace the test.
function is_nonfinite_trace(x) result(nonfinite)
  real, intent(in) :: x
    !< A value to test [nondim]
  logical :: nonfinite
    !< True if x is an IEEE Inf or NaN value [nondim]

  integer(kind=int64) :: xb
    !< Bit-packed real number into integer form [nondim]
  integer(kind=int64) :: exp_bits
    !< Raw exponent bits from xb [nondim]
  integer(kind=int64) :: nonfinite_bits
    !< Raw exponent bits for Inf or NaN values [nondim]

  xb = transfer(x, 1_int64)
  call trace_int('is_nonfinite xb', xb)
  exp_bits = ibits(xb, expbit, expwidth)
  call trace_int('is_nonfinite exp_bits', exp_bits)
  nonfinite_bits = ishft(1_int64, expwidth) - 1_int64
  call trace_int('is_nonfinite all ones', nonfinite_bits)
  nonfinite = (exp_bits == nonfinite_bits)
  write(*,'(a32,1x,l1)') 'is_nonfinite result', nonfinite
end function is_nonfinite_trace


!> Select the root estimate with the smallest residual from three adjacent floats.
function select_best_root_trace(root, x, n) result(best_root)
  real, intent(in) :: root
    !< Initial root estimate [A]
  real, intent(in) :: x
    !< Value whose nth root is being estimated [A^n]
  integer, intent(in) :: n
    !< The degree of the root [nondim]
  real :: best_root
    !< Adjacent root estimate with the smallest residual [A]

  real :: trial_root
    !< The adjacent root estimate being tested [A]
  real :: resid
    !< Signed residual of the initial estimate [A^n]
  real :: best_resid
    !< The smallest residual found so far [A^n]
  real :: trial_resid
    !< The residual for a trial root estimate [A^n]

  best_root = root
  call trace_real('select best_root init', best_root)
  resid = power_residual_trace(x, best_root, n, 'select residual')
  call trace_real('select resid assigned', resid)
  best_resid = abs(resid)
  call trace_real('select best_resid', best_resid)

  if (resid > 0.0) then
    trial_root = nearest(root, 1.0)
    call trace_real('select trial nearest up', trial_root)
  else
    trial_root = nearest(root, -1.0)
    call trace_real('select trial nearest down', trial_root)
  endif
  trial_resid = abs(power_residual_trace(x, trial_root, n, 'select trial residual'))
  call trace_real('select trial_resid', trial_resid)
  if (trial_resid < best_resid) then
    best_root = trial_root
    call trace_real('select accepted trial', best_root)
  else
    call trace_real('select kept original', best_root)
  endif
end function select_best_root_trace


!> Compute x - root**n with a compensated product and trace each step.
function power_residual_trace(x, root, n, context) result(resid)
  real, intent(in) :: x
    !< Value whose nth root is being estimated [A^n]
  real, intent(in) :: root
    !< Root estimate [A]
  integer, intent(in) :: n
    !< The degree of the root [nondim]
  character(len=*), intent(in) :: context
    !< Trace context label [nondim]
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
  real :: p_fma
    !< Fused multiply-add residual for p_hi * root [A^m]
  real :: p_lo_root
    !< Low-order product contribution [A^m]
  integer :: m
    !< Power counter [nondim]

  write(*,'(/,a,1x,a)') '--- power_residual', trim(context)
  p_hi = 1.0
  call trace_real('power p_hi init', p_hi)
  p_lo = 0.0
  call trace_real('power p_lo init', p_lo)
  do m=1,n
    write(*,'(a,i0)') 'power m = ', m
    p_new = p_hi * root
    call trace_real('power p_new', p_new)
    p_fma = ieee_fma(p_hi, root, -p_new)
    call trace_real('power p_fma', p_fma)
    p_lo_root = p_lo * root
    call trace_real('power p_lo_root', p_lo_root)
    p_err = p_fma + p_lo_root
    call trace_real('power p_err', p_err)
    p_hi = p_new + p_err
    call trace_real('power p_hi', p_hi)
    p_lo = p_err - (p_hi - p_new)
    call trace_real('power p_lo', p_lo)
  enddo

  resid = (x - p_hi) - p_lo
  call trace_real('power resid', resid)
end function power_residual_trace


!> Raise a real number to a small positive integer power with traced products.
function integer_power_trace(x, n, context) result(xn)
  real, intent(in) :: x
    !< The value to raise to an integer power [A]
  integer, intent(in) :: n
    !< The exponent of x [nondim]
  character(len=*), intent(in) :: context
    !< Trace context label [nondim]
  real :: xn
    !< x raised to the nth power [A^n]

  real :: x2
    !< x raised to the second power [A^2]
  real :: x4
    !< x raised to the fourth power [A^4]
  integer :: m
    !< Power counter [nondim]

  write(*,'(/,a,1x,a,1x,i0)') '--- integer_power', trim(context), n
  select case (n)
  case (0)
    xn = 1.0
    call trace_real('integer_power case 0', xn)
  case (1)
    xn = x
    call trace_real('integer_power case 1', xn)
  case (2)
    xn = x * x
    call trace_real('integer_power case 2', xn)
  case (3)
    x2 = x * x
    call trace_real('integer_power x*x', x2)
    xn = x2 * x
    call trace_real('integer_power case 3', xn)
  case (4)
    x2 = x * x
    call trace_real('integer_power x*x', x2)
    xn = x2 * x2
    call trace_real('integer_power case 4', xn)
  case (5)
    x2 = x * x
    call trace_real('integer_power x*x', x2)
    x4 = x2 * x2
    call trace_real('integer_power x4', x4)
    xn = x4 * x
    call trace_real('integer_power case 5', xn)
  case default
    xn = 1.0
    call trace_real('integer_power init', xn)
    do m=1,n
      write(*,'(a,i0)') 'integer_power m = ', m
      xn = xn * x
      call trace_real('integer_power xn', xn)
    enddo
  end select
end function integer_power_trace

end program time_debug_nth_root_trace
