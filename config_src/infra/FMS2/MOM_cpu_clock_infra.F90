! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Wraps the MPP cpu clock functions
!!
!! The functions and constants should be accessed via mom_cpu_clock
!!
!! Compiling with -DMOM_USE_NVTX additionally emits an NVTX range around every MOM6 cpu
!! clock, so each existing cpu_clock_id() name becomes a named range in an nsys timeline
!! with no call-site changes.  It requires nvfortran and the NVTX library
!! (-DMOM_USE_NVTX ... -cudalib=nvtx).  Undefined by default: a normal build compiles
!! exactly as before and links no extra library.
module MOM_cpu_clock_infra

! These interfaces and constants from MPP/FMS will not be directly exposed outside of this module
use fms_mod, only : clock_flag_default
use mpp_mod, only : mpp_clock_begin
use mpp_mod, only : mpp_clock_end, mpp_clock_id
use mpp_mod, only : MPP_CLOCK_COMPONENT => CLOCK_COMPONENT
use mpp_mod, only : MPP_CLOCK_SUBCOMPONENT => CLOCK_SUBCOMPONENT
use mpp_mod, only : MPP_CLOCK_MODULE_DRIVER => CLOCK_MODULE_DRIVER
use mpp_mod, only : MPP_CLOCK_MODULE => CLOCK_MODULE
use mpp_mod, only : MPP_CLOCK_ROUTINE => CLOCK_ROUTINE
use mpp_mod, only : MPP_CLOCK_LOOP => CLOCK_LOOP
use mpp_mod, only : MPP_CLOCK_INFRA => CLOCK_INFRA
#ifdef MOM_USE_NVTX
use nvtx, only : nvtxStartRange, nvtxEndRange
#endif

implicit none ; private

#ifdef MOM_USE_NVTX
!> The largest clock handle for which an NVTX range name is retained.
integer, parameter :: MAX_NVTX_CLOCKS = 4096
!> The NVTX range name for each clock handle, recorded by cpu_clock_id().  An empty entry
!! means no range is emitted for that handle.  cpu_clock_begin() and cpu_clock_end() test
!! the same condition, so starts and ends stay balanced for handles that were never named
!! or that fall outside the table.
character(len=64), dimension(MAX_NVTX_CLOCKS) :: nvtx_clock_names = ""
#endif

! Public entities
public :: cpu_clock_id, cpu_clock_begin, cpu_clock_end
public :: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER, CLOCK_MODULE
public :: CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! component, e.g. the entire MOM6 model
integer, parameter :: CLOCK_COMPONENT = MPP_CLOCK_COMPONENT

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! sub-component, e.g. dynamics or thermodynamics
integer, parameter :: CLOCK_SUBCOMPONENT = MPP_CLOCK_SUBCOMPONENT

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! module driver, e.g. a routine that calls multiple other routines
integer, parameter :: CLOCK_MODULE_DRIVER = MPP_CLOCK_MODULE_DRIVER

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! module, e.g. the main entry routine for a module
integer, parameter :: CLOCK_MODULE = MPP_CLOCK_MODULE

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! subroutine or function
integer, parameter :: CLOCK_ROUTINE = MPP_CLOCK_ROUTINE

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! section with in a routine, e.g. around a loop
integer, parameter :: CLOCK_LOOP = MPP_CLOCK_LOOP

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for an
!! infrastructure operation, e.g. a halo update
integer, parameter :: CLOCK_INFRA = MPP_CLOCK_INFRA

contains

!> Turns on clock with handle "id"
subroutine cpu_clock_begin(id)
  integer, intent(in) :: id !< Handle for clock

#ifdef MOM_USE_NVTX
  ! Opened before, and closed after, the mpp clock so the NVTX range encloses it.
  if (id > 0 .and. id <= MAX_NVTX_CLOCKS) then
    if (len_trim(nvtx_clock_names(id)) > 0) call nvtxStartRange(trim(nvtx_clock_names(id)))
  endif
#endif
  call mpp_clock_begin(id)

end subroutine cpu_clock_begin

!> Turns off clock with handle "id"
subroutine cpu_clock_end(id)
  integer, intent(in) :: id !< Handle for clock

  call mpp_clock_end(id)
#ifdef MOM_USE_NVTX
  if (id > 0 .and. id <= MAX_NVTX_CLOCKS) then
    if (len_trim(nvtx_clock_names(id)) > 0) call nvtxEndRange
  endif
#endif

end subroutine cpu_clock_end

!> Returns the integer handle for a named CPU clock.
integer function cpu_clock_id(name, sync, grain)
  character(len=*),  intent(in) :: name  !< The unique name of the CPU clock
  logical, optional, intent(in) :: sync !< A flag that controls whether the
                  !! PEs are synchronized before the cpu clocks start counting.
                  !! Synchronization occurs before the start of a clock if this
                  !! is enabled, while additional (expensive) statistics can
                  !! set for other values.
                  !! If absent, the default is taken from the settings for FMS.
  integer, optional, intent(in) :: grain !< The timing granularity for this clock, usually set to
                                       !! the values of CLOCK_COMPONENT, CLOCK_ROUTINE, CLOCK_LOOP, etc.

  integer :: clock_flags

  clock_flags = clock_flag_default
  if (present(sync)) then
    if (sync) then
      clock_flags = ibset(clock_flags, 0)
    else
      clock_flags = ibclr(clock_flags, 0)
    endif
  endif

  cpu_clock_id = mpp_clock_id(name, flags=clock_flags, grain=grain)
#ifdef MOM_USE_NVTX
  if (cpu_clock_id > 0 .and. cpu_clock_id <= MAX_NVTX_CLOCKS) &
    nvtx_clock_names(cpu_clock_id) = name
#endif
end function cpu_clock_id

end module MOM_cpu_clock_infra
