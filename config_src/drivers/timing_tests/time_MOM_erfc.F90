! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Regional throughput and dependent-call latency benchmarks for erfc_repro().
!! Array timings include input addition and an integer checksum; the copy
!! baseline measures that overhead. Scalar calls depend on the preceding result
!! through a load from one of two identical input copies; both kernels visit the
!! same values in the same order. No fast-math options should be used.
program time_MOM_erfc

use, intrinsic :: iso_fortran_env, only : int64, real64
use, intrinsic :: iso_fortran_env, only : compiler_version, compiler_options
use MOM_intrinsic_functions, only : erfc_repro

implicit none

integer, parameter :: npts = 16384, trials = 5, ncases = 8
  !< Array length (power of two), timing repetitions, and workloads
real(real64), parameter :: target_seconds = 0.05_real64
  !< Minimum calibrated duration of a timing trial [s]
character(len=16), parameter :: labels(ncases) = [character(len=16) :: &
    'near zero', 'middle positive', 'tail normal', 'tail subnormal', &
    'negative', 'mixed [-6,28]', 'saturated', 'hard midpoints']
real, parameter :: lower(ncases-1) = [-0.46875, 0.46875, 4., 26.55, -6., -6., 28.]
  !< Lower workload bounds [nondim]
real, parameter :: upper(ncases-1) = [0.46875, 4., 26.5, 27.22, -0.46875, 28., 100.]
  !< Upper workload bounds [nondim]
integer(int64), parameter :: hard_bits(8) = [ &
    int(z'40348DE452FB1A15',int64), int(z'BC9C5BF891B4EF6B',int64), &
    int(z'BCCFE777A3EB8D58',int64), int(z'BDAB71F3A10691D9',int64), &
    int(z'BFFF9A4A209CA0E4',int64), int(z'3C8C5BF891B4EF6B',int64), &
    int(z'3FBD4AF8ADB90116',int64), int(z'403A8F7BFBD15495',int64)]
  !< Difficult binary64 inputs, also covered by the accuracy tests
integer, parameter :: realquad = selected_real_kind(p=30, r=300)
  !< Optional quad kind for untimed accuracy checks
integer, parameter :: realq = merge(realquad, kind(1.), realquad > 0)
  !< Valid kind even when quad is unavailable
real :: x(2*npts), fraction ! Duplicated input array and deterministic random fraction [nondim]
real, volatile :: perturbation ! Opaque zero preventing invariant-call hoisting [nondim]
real(real64) :: seconds ! Elapsed time [s]
real(real64) :: samples(trials,3), medians(3), tmp
  ! Sample/median per-call times and sorting temporary [ns]
real(real64) :: spread ! Full trial range relative to the median [%]
integer(int64) :: rate, before, after, checksum, digest, state
  ! Clock rate/ticks, anti-elision checksum, and bounded random state
integer :: workload, mode, kernel, repetition, pass, loops(3), n, i, j
  ! Workload, timing mode, kernel, repetition, calibration and array indices

if (storage_size(0.) /= 64) error stop 'time_MOM_erfc requires binary64 default real'
call system_clock(count_rate=rate)
perturbation = 0.
print '(a)', trim(compiler_version())
print '(a)', trim(compiler_options())
print '(a,i0,a,f5.2,a)', 'Median of ', trials, ' trials, each calibrated to >= ', target_seconds, ' s'
print '(a)', 'All timings are raw ns/call; copy baseline includes the common harness overhead.'
print '(a)', 'Hard-midpoint timings are an adversarial workload, not a normal fallback frequency.'

do mode=1,2
  state = 12345_int64
  if (mode == 1) then
    print '(/,a)', 'Array throughput (independent inputs)'
  else
    print '(/,a)', 'Scalar latency (result-dependent input index)'
  endif
  print '(a)', 'workload            copy ns     libm ns    repro ns     ratio   repro spread %'
  do workload=1,ncases
    digest = 0_int64
    n = npts
    if (workload == ncases) n = 128
    do i=1,n
      if (workload == ncases) then
        x(i) = transfer(hard_bits(1+modulo(i-1,size(hard_bits))), 0.)
      else
        state = modulo(16807_int64 * state, 2147483647_int64)
        fraction = real(state) / 2147483647.
        x(i) = lower(workload) + (upper(workload)-lower(workload)) * fraction
        if (workload == 7 .and. modulo(i,2) == 0) x(i) = -x(i)
      endif
      x(n+i) = x(i)
    enddo

    do kernel=1,3
      loops(kernel) = 1
      do
        call system_clock(count=before)
        call batch(kernel, mode, x(1:2*n), loops(kernel), checksum)
        call system_clock(count=after)
        digest = ieor(ishftc(digest,1), checksum)
        seconds = real(after-before,real64) / real(rate,real64)
        if (seconds >= target_seconds) exit
        loops(kernel) = 2 * loops(kernel)
      enddo
    enddo
    ! Rotate kernel order to reduce systematic bias from clock and temperature drift.
    do repetition=1,trials
      do pass=1,3
        kernel = 1 + modulo(pass+repetition-2,3)
        call system_clock(count=before)
        call batch(kernel, mode, x(1:2*n), loops(kernel), checksum)
        call system_clock(count=after)
        digest = ieor(ishftc(digest,1), checksum)
        samples(repetition,kernel) = (1.e9_real64 * real(after-before,real64)) / &
            ((real(rate,real64) * real(n,real64)) * real(loops(kernel),real64))
      enddo
    enddo
    do kernel=1,3
      do i=2,trials
        tmp = samples(i,kernel)
        j = i-1
        do while (j >= 1)
          if (samples(j,kernel) <= tmp) exit
          samples(j+1,kernel) = samples(j,kernel)
          j = j-1
        enddo
        samples(j+1,kernel) = tmp
      enddo
      medians(kernel) = samples((trials+1)/2,kernel)
    enddo
    spread = 100. * (samples(trials,3)-samples(1,3)) / medians(3)
    print '(a18,3f12.2,f10.2,f16.1)', labels(workload), medians(:), medians(3)/medians(2), spread
    ! This observable value depends on every output of every timed repetition.
    print '(a,z16.16)', '  checksum: ', digest
    if (mode == 1 .and. realquad > 0 .and. workload < ncases) call array_accuracy(x(1:n))
  enddo
enddo

contains

!> Run a timed kernel; the volatile shift and output checksum prevent elision.
subroutine batch(kernel, mode, inputs, iterations, checksum)
  integer, intent(in) :: kernel !< 1: copy, 2: libm erfc, 3: erfc_repro
  integer, intent(in) :: mode !< 1: independent arrays, 2: dependent scalar calls
  real, intent(in) :: inputs(:) !< Two identical copies of the workload inputs [nondim]
  integer, intent(in) :: iterations !< Number of batches
  integer(int64), intent(out) :: checksum !< Observable digest of every result
  real :: values(size(inputs)/2), value, shift ! Results and opaque input shift [nondim]
  integer(int64) :: bits ! Result bits used by the digest and scalar dependency
  integer :: i, j, index, offset, n ! Loop and input indices and duplicate-copy offset

  checksum = 0_int64
  n = size(inputs)/2
  offset = 0
  do j=1,iterations
    shift = perturbation
    if (mode == 1) then
      select case (kernel)
      case (1)
        values(:) = inputs(1:n) + shift
      case (2)
        values(:) = erfc(inputs(1:n) + shift)
      case (3)
        values(:) = erfc_repro(inputs(1:n) + shift)
      end select
      do i=1,n
        checksum = ieor(checksum, transfer(values(i),0_int64))
      enddo
    else
      do i=1,n
        index = i + offset
        select case (kernel)
        case (1)
          value = inputs(index) + shift
        case (2)
          value = erfc(inputs(index) + shift)
        case (3)
          value = erfc_repro(inputs(index) + shift)
        end select
        bits = transfer(value,0_int64)
        offset = n * int(iand(bits,1_int64))
        checksum = ieor(checksum, bits)
      enddo
    endif
    checksum = ieor(ishftc(checksum,1), int(j,int64))
  enddo
end subroutine batch

!> Check the actual intrinsic array path outside the timed region.
!! GNU may select libmvec here, whose accuracy differs from scalar libm erfc.
subroutine array_accuracy(inputs)
  real, intent(in) :: inputs(:) !< Benchmark input values [nondim]
  real :: values(size(inputs)), rounded, ulp ! Array results, rounded reference, and spacing [nondim]
  real(kind=realq) :: reference, error, worst ! Quad reference and fractional ULP errors [nondim]
  integer :: i, different ! Input index and count differing from the rounded quad reference

  values(:) = erfc(inputs(:))
  worst = 0.
  different = 0
  do i=1,size(inputs)
    reference = erfc(real(inputs(i),realq))
    rounded = real(reference)
    ulp = nearest(rounded,1.) - rounded
    error = abs(real(values(i),realq)-reference) / real(ulp,realq)
    worst = max(worst,error)
    if (values(i) /= rounded) different = different+1
  enddo
  print '(a,es12.4,a,i0,a,i0)', '  intrinsic array max ULP: ', real(worst,real64), &
      '; differs from rounded quad: ', different, '/', size(inputs)
end subroutine array_accuracy

end program time_MOM_erfc
