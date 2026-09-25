! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains routines that implement physical fluxes of tracers (e.g. due
!! to surface fluxes or mixing). These are intended to be called from call_tracer_column_fns
!! in the MOM_tracer_flow_control module.
module MOM_tracer_diabatic

use MOM_grid,             only : ocean_grid_type
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_forcing_type,     only : forcing
use MOM_error_handler,    only : MOM_error, FATAL, WARNING

implicit none ; private

#include <MOM_memory.h>
#include "do_concurrent_compat.h"
public tracer_vertdiff, tracer_vertdiff_Eulerian
public applyTracerBoundaryFluxesInOut

contains

!> This subroutine solves a tridiagonal equation for the final tracer concentrations after the
!! dual-entrainments, and possibly sinking or surface and bottom sources, are applied.  The sinking
!! is implemented with an fully implicit upwind advection scheme.  Alternate time units can be
!! used for the timestep, surface and bottom fluxes and sink_rate provided they are all consistent.
subroutine tracer_vertdiff(h_old, ea, eb, dt, tr, G, GV, &
                           sfc_flux, btm_flux, btm_reservoir, sink_rate, convert_flux_in)
  type(ocean_grid_type),                     intent(in)    :: G      !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old  !< layer thickness before entrainment
                                                                     !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: ea     !< amount of fluid entrained from the layer
                                                                     !! above [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: eb     !< amount of fluid entrained from the layer
                                                                     !! below [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: tr     !< tracer concentration in concentration units [CU]
  real,                                      intent(in)    :: dt     !< amount of time covered by this call [T ~> s]
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: sfc_flux !< surface flux of the tracer in units of
                                                                     !! [CU R Z T-1 ~> CU kg m-2 s-1] or
                                                                     !! [CU H ~> CU m or CU kg m-2] if
                                                                     !! convert_flux_in is .false.
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: btm_flux !< The (negative upward) bottom flux of the
                                                                     !! tracer in [CU R Z T-1 ~> CU kg m-2 s-1] or
                                                                     !! [CU H ~> CU m or CU kg m-2] if
                                                                     !! convert_flux_in is .false.
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(inout) :: btm_reservoir !< amount of tracer in a bottom reservoir
                                                                     !! [CU R Z ~> CU kg m-2]
  real,                             optional,intent(in)    :: sink_rate !< rate at which the tracer sinks
                                                                     !! [Z T-1 ~> m s-1]
  logical,                          optional,intent(in)    :: convert_flux_in !< True if the specified sfc_flux needs
                                                                     !! to be integrated in time

  ! local variables
  real :: sink_dist !< The distance the tracer sinks in a time step [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    sfc_src, &      !< The time-integrated surface source of the tracer [CU H ~> CU m or CU kg m-2].
    btm_src         !< The time-integrated bottom source of the tracer [CU H ~> CU m or CU kg m-2].
  real, dimension(SZI_(G)) :: &
    b1, &           !< b1 is used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
    d1              !! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: c1(SZI_(G),SZK_(GV))    !< c1 is used by the tridiagonal solver [nondim].
  real :: h_minus_dsink(SZI_(G),SZK_(GV)) !< The layer thickness minus the
                    !! difference in sinking rates across the layer [H ~> m or kg m-2].
                    !! By construction, 0 <= h_minus_dsink < h_work.
  real :: sink(SZI_(G),SZK_(GV)+1) !< The tracer's sinking distances at the
                    !! interfaces, limited to prevent characteristics from
                    !! crossing within a single timestep [H ~> m or kg m-2].
  real :: b_denom_1 !< The first term in the denominator of b1 [H ~> m or kg m-2].
  real :: h_tr      !< h_tr is h at tracer points with a h_neglect added to
                    !! ensure positive definiteness [H ~> m or kg m-2].
  real :: h_neglect !< A thickness that is so small it is usually lost
                    !! in roundoff and can be neglected [H ~> m or kg m-2].
  logical :: convert_flux = .true.

  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (nz == 1) then
    call MOM_error(WARNING, "MOM_tracer_diabatic.F90, tracer_vertdiff called "//&
                            "with only one vertical level")
    return
  endif

  if (present(convert_flux_in)) convert_flux = convert_flux_in
  h_neglect = GV%H_subroundoff
  sink_dist = 0.0
  if (present(sink_rate)) sink_dist = (dt*sink_rate) * GV%Z_to_H
  !$OMP parallel default(shared) private(sink,h_minus_dsink,b_denom_1,b1,d1,h_tr,c1)
  !$OMP do
  do j=js,je ; do i=is,ie ; sfc_src(i,j) = 0.0 ; btm_src(i,j) = 0.0 ; enddo ; enddo
  if (present(sfc_flux)) then
    if (convert_flux) then
      !$OMP do
      do j=js,je ; do i=is,ie
        sfc_src(i,j) = (sfc_flux(i,j)*dt) * GV%RZ_to_H
      enddo ; enddo
    else
      !$OMP do
      do j=js,je ; do i=is,ie
        sfc_src(i,j) = sfc_flux(i,j)
      enddo ; enddo
    endif
  endif
  if (present(btm_flux)) then
    if (convert_flux) then
      !$OMP do
      do j=js,je ; do i=is,ie
        btm_src(i,j) = (btm_flux(i,j)*dt) * GV%RZ_to_H
      enddo ; enddo
    else
      !$OMP do
      do j=js,je ; do i=is,ie
        btm_src(i,j) = btm_flux(i,j)
      enddo ; enddo
    endif
  endif

  if (present(sink_rate)) then
    !$OMP do
    do j=js,je
      ! Find the sinking rates at all interfaces, limiting them if necesary
      ! so that the characteristics do not cross within a timestep.
      !   If a non-constant sinking rate were used, that would be incorprated
      ! here.
      if (present(btm_reservoir)) then
        do i=is,ie ; sink(i,nz+1) = sink_dist ; enddo
        do k=2,nz ; do i=is,ie
          sink(i,K) = sink_dist ; h_minus_dsink(i,k) = h_old(i,j,k)
        enddo ; enddo
      else
        do i=is,ie ; sink(i,nz+1) = 0.0 ; enddo
        ! Find the limited sinking distance at the interfaces.
        do k=nz,2,-1 ; do i=is,ie
          if (sink(i,K+1) >= sink_dist) then
            sink(i,K) = sink_dist
            h_minus_dsink(i,k) = h_old(i,j,k) + (sink(i,K+1) - sink(i,K))
          elseif (sink(i,K+1) + h_old(i,j,k) < sink_dist) then
            sink(i,K) = sink(i,K+1) + h_old(i,j,k)
            h_minus_dsink(i,k) = 0.0
          else
            sink(i,K) = sink_dist
            h_minus_dsink(i,k) = (h_old(i,j,k) + sink(i,K+1)) - sink(i,K)
          endif
        enddo ; enddo
      endif
      do i=is,ie
        sink(i,1) = 0.0 ; h_minus_dsink(i,1) = (h_old(i,j,1) + sink(i,2))
      enddo

      ! Now solve the tridiagonal equation for the tracer concentrations.
      do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        b_denom_1 = h_minus_dsink(i,1) + ea(i,j,1) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,1) + h_neglect
        tr(i,j,1) = (b1(i)*h_tr)*tr(i,j,1) + b1(i)*sfc_src(i,j)
      endif ; enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        c1(i,k) = eb(i,j,k-1) * b1(i)
        b_denom_1 = h_minus_dsink(i,k) + d1(i) * (ea(i,j,k) + sink(i,K)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + &
                             (ea(i,j,k) + sink(i,K)) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        c1(i,nz) = eb(i,j,nz-1) * b1(i)
        b_denom_1 = h_minus_dsink(i,nz) + d1(i) * (ea(i,j,nz) + sink(i,nz)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,nz))
        h_tr = h_old(i,j,nz) + h_neglect
        tr(i,j,nz) = b1(i) * ((h_tr * tr(i,j,nz) + btm_src(i,j)) + &
                              (ea(i,j,nz) + sink(i,nz)) * tr(i,j,nz-1))
      endif ; enddo
      if (present(btm_reservoir)) then ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        btm_reservoir(i,j) = btm_reservoir(i,j) + (sink(i,nz+1)*tr(i,j,nz)) * GV%H_to_RZ
      endif ; enddo ; endif

      do k=nz-1,1,-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  else
    !$OMP do
    do j=js,je
      do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        h_tr = h_old(i,j,1) + h_neglect
        b_denom_1 = h_tr + ea(i,j,1)
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = h_tr * b1(i)
        tr(i,j,1) = (b1(i)*h_tr)*tr(i,j,1) + b1(i)*sfc_src(i,j)
      endif ; enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        c1(i,k) = eb(i,j,k-1) * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        b_denom_1 = h_tr + d1(i) * ea(i,j,k)
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + ea(i,j,k) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        c1(i,nz) = eb(i,j,nz-1) * b1(i)
        h_tr = h_old(i,j,nz) + h_neglect
        b_denom_1 = h_tr + d1(i)*ea(i,j,nz)
        b1(i) = 1.0 / ( b_denom_1 + eb(i,j,nz))
        tr(i,j,nz) = b1(i) * (( h_tr * tr(i,j,nz) + btm_src(i,j)) + &
                              ea(i,j,nz) * tr(i,j,nz-1))
      endif ; enddo
      do k=nz-1,1,-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  endif
  !$OMP end parallel

end subroutine tracer_vertdiff


!> This subroutine solves a tridiagonal equation for the final tracer concentrations after
!! Eulerian mixing, and possibly sinking or surface and bottom sources, are applied.  The sinking
!! is implemented with an fully implicit upwind advection scheme.  Alternate time units can be
!! used for the timestep, surface and bottom fluxes and sink_rate provided they are all consistent.
subroutine tracer_vertdiff_Eulerian(h_old, ent, dt, tr, G, GV, &
                                    sfc_flux, btm_flux, btm_reservoir, sink_rate, convert_flux_in)
  type(ocean_grid_type),                     intent(in)    :: G      !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h_old  !< layer thickness before entrainment
                                                                     !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: ent    !< Amount of fluid mixed across interfaces
                                                                     !! [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: tr     !< tracer concentration in concentration units [CU]
  real,                                      intent(in)    :: dt     !< amount of time covered by this call [T ~> s]
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: sfc_flux !< surface flux of the tracer in units of
                                                                     !! [CU R Z T-1 ~> CU kg m-2 s-1] or
                                                                     !! [CU H ~> CU m or CU kg m-2] if
                                                                     !! convert_flux_in is .false.
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: btm_flux !< The (negative upward) bottom flux of the
                                                                     !! tracer in [CU kg m-2 T-1 ~> CU kg m-2 s-1] or
                                                                     !! [CU H ~> CU m or CU kg m-2] if
                                                                     !! convert_flux_in is .false.
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(inout) :: btm_reservoir !< amount of tracer in a bottom reservoir
                                                                     !! [CU R Z ~> CU kg m-2]
  real,                             optional,intent(in)    :: sink_rate !< rate at which the tracer sinks
                                                                     !! [Z T-1 ~> m s-1]
  logical,                          optional,intent(in)    :: convert_flux_in !< True if the specified sfc_flux needs
                                                                     !! to be integrated in time

  ! local variables
  real :: sink_dist !< The distance the tracer sinks in a time step [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    sfc_src, &      !< The time-integrated surface source of the tracer [CU H ~> CU m or CU kg m-2].
    btm_src         !< The time-integrated bottom source of the tracer [CU H ~> CU m or CU kg m-2].
  real, dimension(SZI_(G)) :: &
    b1, &           !< b1 is used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
    d1              !! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: c1(SZI_(G),SZK_(GV))    !< c1 is used by the tridiagonal solver [nondim].
  real :: h_minus_dsink(SZI_(G),SZK_(GV)) !< The layer thickness minus the
                    !! difference in sinking rates across the layer [H ~> m or kg m-2].
                    !! By construction, 0 <= h_minus_dsink < h_work.
  real :: sink(SZI_(G),SZK_(GV)+1) !< The tracer's sinking distances at the
                    !! interfaces, limited to prevent characteristics from
                    !! crossing within a single timestep [H ~> m or kg m-2].
  real :: b_denom_1 !< The first term in the denominator of b1 [H ~> m or kg m-2].
  real :: h_tr      !< h_tr is h at tracer points with a h_neglect added to
                    !! ensure positive definiteness [H ~> m or kg m-2].
  real :: h_neglect !< A thickness that is so small it is usually lost
                    !! in roundoff and can be neglected [H ~> m or kg m-2].
  logical :: convert_flux

  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (nz == 1) then
    call MOM_error(WARNING, "MOM_tracer_diabatic.F90, tracer_vertdiff called "//&
                            "with only one vertical level")
    return
  endif

  convert_flux = .true.
  if (present(convert_flux_in)) convert_flux = convert_flux_in
  h_neglect = GV%H_subroundoff
  sink_dist = 0.0
  if (present(sink_rate)) sink_dist = (dt*sink_rate) * GV%Z_to_H
  !$OMP parallel default(shared) private(sink,h_minus_dsink,b_denom_1,b1,d1,h_tr,c1)
  !$OMP do
  do j=js,je ; do i=is,ie ; sfc_src(i,j) = 0.0 ; btm_src(i,j) = 0.0 ; enddo ; enddo
  if (present(sfc_flux)) then
    if (convert_flux) then
      !$OMP do
      do j=js,je ; do i=is,ie
        sfc_src(i,j) = (sfc_flux(i,j)*dt) * GV%RZ_to_H
      enddo ; enddo
    else
      !$OMP do
      do j=js,je ; do i=is,ie
        sfc_src(i,j) = sfc_flux(i,j)
      enddo ; enddo
    endif
  endif
  if (present(btm_flux)) then
    if (convert_flux) then
      !$OMP do
      do j=js,je ; do i=is,ie
        btm_src(i,j) = (btm_flux(i,j)*dt) * GV%kg_m2_to_H
      enddo ; enddo
    else
      !$OMP do
      do j=js,je ; do i=is,ie
        btm_src(i,j) = btm_flux(i,j)
      enddo ; enddo
    endif
  endif

  if (present(sink_rate)) then
    !$OMP do
    do j=js,je
      ! Find the sinking rates at all interfaces, limiting them if necesary
      ! so that the characteristics do not cross within a timestep.
      !   If a non-constant sinking rate were used, that would be incorprated
      ! here.
      if (present(btm_reservoir)) then
        do i=is,ie ; sink(i,nz+1) = sink_dist ; enddo
        do k=2,nz ; do i=is,ie
          sink(i,K) = sink_dist ; h_minus_dsink(i,k) = h_old(i,j,k)
        enddo ; enddo
      else
        do i=is,ie ; sink(i,nz+1) = 0.0 ; enddo
        ! Find the limited sinking distance at the interfaces.
        do k=nz,2,-1 ; do i=is,ie
          if (sink(i,K+1) >= sink_dist) then
            sink(i,K) = sink_dist
            h_minus_dsink(i,k) = h_old(i,j,k) + (sink(i,K+1) - sink(i,K))
          elseif (sink(i,K+1) + h_old(i,j,k) < sink_dist) then
            sink(i,K) = sink(i,K+1) + h_old(i,j,k)
            h_minus_dsink(i,k) = 0.0
          else
            sink(i,K) = sink_dist
            h_minus_dsink(i,k) = (h_old(i,j,k) + sink(i,K+1)) - sink(i,K)
          endif
        enddo ; enddo
      endif
      do i=is,ie
        sink(i,1) = 0.0 ; h_minus_dsink(i,1) = (h_old(i,j,1) + sink(i,2))
      enddo

      ! Now solve the tridiagonal equation for the tracer concentrations.
      do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        b_denom_1 = h_minus_dsink(i,1) + ent(i,j,1) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + ent(i,j,2))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,1) + h_neglect
        tr(i,j,1) = (b1(i)*h_tr)*tr(i,j,1) + b1(i)*sfc_src(i,j)
      endif ; enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        c1(i,k) = ent(i,j,K) * b1(i)
        b_denom_1 = h_minus_dsink(i,k) + d1(i) * (ent(i,j,K) + sink(i,K)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + ent(i,j,K+1))
        d1(i) = b_denom_1 * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + &
                             (ent(i,j,K) + sink(i,K)) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        c1(i,nz) = ent(i,j,nz) * b1(i)
        b_denom_1 = h_minus_dsink(i,nz) + d1(i) * (ent(i,j,nz) + sink(i,nz)) + &
                    h_neglect
        b1(i) = 1.0 / (b_denom_1 + ent(i,j,nz+1))
        h_tr = h_old(i,j,nz) + h_neglect
        tr(i,j,nz) = b1(i) * ((h_tr * tr(i,j,nz) + btm_src(i,j)) + &
                              (ent(i,j,nz) + sink(i,nz)) * tr(i,j,nz-1))
      endif ; enddo
      if (present(btm_reservoir)) then ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        btm_reservoir(i,j) = btm_reservoir(i,j) + (sink(i,nz+1)*tr(i,j,nz)) * GV%H_to_RZ
      endif ; enddo ; endif

      do k=nz-1,1,-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  else
    !$OMP do
    do j=js,je
      do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        h_tr = h_old(i,j,1) + h_neglect
        b_denom_1 = h_tr + ent(i,j,1)
        b1(i) = 1.0 / (b_denom_1 + ent(i,j,2))
        d1(i) = h_tr * b1(i)
        tr(i,j,1) = (b1(i)*h_tr)*tr(i,j,1) + b1(i)*sfc_src(i,j)
      endif ; enddo
      do k=2,nz-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        c1(i,k) = ent(i,j,K) * b1(i)
        h_tr = h_old(i,j,k) + h_neglect
        b_denom_1 = h_tr + d1(i) * ent(i,j,K)
        b1(i) = 1.0 / (b_denom_1 + ent(i,j,K+1))
        d1(i) = b_denom_1 * b1(i)
        tr(i,j,k) = b1(i) * (h_tr * tr(i,j,k) + ent(i,j,K) * tr(i,j,k-1))
      endif ; enddo ; enddo
      do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        c1(i,nz) = ent(i,j,nz) * b1(i)
        h_tr = h_old(i,j,nz) + h_neglect
        b_denom_1 = h_tr + d1(i)*ent(i,j,nz)
        b1(i) = 1.0 / ( b_denom_1 + ent(i,j,nz+1))
        tr(i,j,nz) = b1(i) * (( h_tr * tr(i,j,nz) + btm_src(i,j)) + &
                              ent(i,j,nz) * tr(i,j,nz-1))
      endif ; enddo
      do k=nz-1,1,-1 ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
        tr(i,j,k) = tr(i,j,k) + c1(i,k+1)*tr(i,j,k+1)
      endif ; enddo ; enddo
    enddo
  endif
  !$OMP end parallel

end subroutine tracer_vertdiff_Eulerian


!> This routine is modeled after applyBoundaryFluxesInOut in MOM_diabatic_aux.F90
!! NOTE: Please note that in this routine sfc_flux gets set to zero to ensure that the surface
!! flux of the tracer does not get applied again during a subsequent call to tracer_vertdif
subroutine applyTracerBoundaryFluxesInOut(G, GV, Tr, dt, fluxes, h, evap_CFL_limit, minimum_forcing_depth, &
               in_flux_optional, out_flux_optional, update_h_opt)

  type(ocean_grid_type),                      intent(in   ) :: G  !< Grid structure
  type(verticalGrid_type),                    intent(in   ) :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: Tr !< Tracer concentration on T-cell [conc]
  real,                                       intent(in   ) :: dt !< Time-step over which forcing is applied [T ~> s]
  type(forcing),                              intent(in   ) :: fluxes !< Surface fluxes container
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  real,                                       intent(in   ) :: evap_CFL_limit !< Limit on the fraction of the
                                                                  !! water that can be fluxed out of the top
                                                                  !! layer in a timestep [nondim]
  real,                                       intent(in   ) :: minimum_forcing_depth !< The smallest depth over
                                                                  !! which fluxes can be applied [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in   ) :: in_flux_optional !< The total time-integrated
                                                                  !! amount of tracer that enters with freshwater
                                                                  !! [conc H ~> conc m or conc kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: out_flux_optional !< The total time-integrated
                                                                  !! amount of tracer that leaves with freshwater
                                                                  !! [conc H ~> conc m or conc kg m-2]
  logical,                          optional, intent(in) :: update_h_opt  !< Optional flag to determine whether
                                                                  !! h should be updated

  integer, parameter :: maxGroundings = 5
  integer :: numberOfGroundings, iGround(maxGroundings), jGround(maxGroundings)
  integer :: nGround           ! The running count of groundings found by the host scan [nondim]
  real :: IforcingDepthScale ! The inverse of the scale over which to apply forcing [H-1 ~> m-1 or m2 kg-1]
  real :: dThickness         ! The change in a layer's thickness [H ~> m or kg m-2]
  real :: dTracer            ! The change in the integrated tracer content of a layer [conc H ~> conc m or conc kg m-2]
  real :: fractionOfForcing  ! The fraction of the forcing to apply to a layer [nondim]
  real :: hOld               ! The layer thickness before surface forcing is applied [H ~> m or kg m-2]
  real :: Ithickness         ! The inverse of the new layer thickness [H-1 ~> m-1 or m2 kg-1]

  real :: h_k                ! A column-local copy of the thickness of layer k [H ~> m or kg m-2]
  real :: Tr_k               ! A column-local copy of the tracer concentration in layer k [conc]
  real :: in_flux(SZI_(G),SZJ_(G))  ! The total time-integrated amount of tracer that
                                    ! enters with freshwater [conc H ~> conc m or conc kg m-2]
  real :: out_flux(SZI_(G),SZJ_(G)) ! The total time-integrated amount of tracer that
                                    ! leaves with freshwater [conc H ~> conc m or conc kg m-2]
  real :: hGround(SZI_(G),SZJ_(G))  ! The amount of tracer in each column that could not be
                                    ! applied because the column grounded out [conc H ~> conc m or conc kg m-2]
  real :: netMassOut_col  ! The remaining mass leaving the ocean surface in a column [H ~> m or kg m-2]
  real :: in_flux_col     ! The remaining amount of tracer that enters with the freshwater in a
                          ! column [conc H ~> conc m or conc kg m-2]
  real :: out_flux_col    ! The remaining amount of tracer that leaves with the freshwater in a
                          ! column [conc H ~> conc m or conc kg m-2]
  real :: hGrounding(maxGroundings) ! The remaining fresh water flux that was not able to be
                               ! supplied from a column that grounded out [H ~> m or kg m-2]
  logical :: update_h
  integer :: i, j, is, ie, js, je, k, nz
  character(len=45) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ! If no freshwater fluxes, nothing needs to be done in this routine
  if ( (.not. associated(fluxes%netMassIn)) .or. (.not. associated(fluxes%netMassOut)) ) return

  !   Tr, h and the fluxes type with its netMassIn and netMassOut arrays are mapped to the device
  ! by the caller.  The optional tracer fluxes are mapped here, because a caller might pass them
  ! as expressions that it can not map by name.
  !$omp target enter data map(alloc: in_flux, out_flux, hGround)
  !$omp target enter data map(to: in_flux_optional) if(present(in_flux_optional))
  !$omp target enter data map(to: out_flux_optional) if(present(out_flux_optional))

  do concurrent (j=js:je, i=is:ie)
    in_flux(i,j) = 0.0 ; out_flux(i,j) = 0.0
  enddo
  if (present(in_flux_optional)) then
    do concurrent (j=js:je, i=is:ie)
      in_flux(i,j) = in_flux_optional(i,j)
    enddo
  endif
  if (present(out_flux_optional)) then
    do concurrent (j=js:je, i=is:ie)
      out_flux(i,j) = out_flux_optional(i,j)
    enddo
  endif

  if (present(update_h_opt)) then
    update_h = update_h_opt
  else
    update_h = .true.
  endif

  numberOfGroundings = 0

  !   Each column is handled by one thread, working down the column with copies of the layer
  ! thickness and tracer concentration.  The surface forcing is contained in the fluxes type.
  ! We aggregate the thermodynamic forcing for a time step into the following:
  ! These should have been set and stored during a call to applyBoundaryFluxesInOut
  ! netMassIn    = net mass entering at ocean surface over a timestep
  ! netMassOut   = net mass leaving ocean surface [H ~> m or kg m-2] over a time step.
  !                netMassOut < 0 means mass leaves ocean.
  ! Note here that the aggregateFW flag has already been taken care of in the call to
  ! applyBoundaryFluxesInOut
  do concurrent (j=js:je, i=is:ie) DO_LOCALITY(local(k, h_k, Tr_k, hOld, Ithickness, dThickness, dTracer)) &
                                 & DO_LOCALITY(local(fractionOfForcing, IforcingDepthScale)) &
                                 & DO_LOCALITY(local(netMassOut_col, in_flux_col, out_flux_col)) &
                                 & DO_LOCALITY(reduce(+:numberOfGroundings))
    in_flux_col = in_flux(i,j)
    out_flux_col = out_flux(i,j)
    if (G%mask2dT(i,j) > 0.) then
      netMassOut_col = fluxes%netMassOut(i,j)

      ! Apply the surface boundary fluxes in two steps:
      ! A/ update concentration from mass entering the ocean, which only changes the top layer
      ! B/ update concentration from mass leaving ocean.
      !   Step A does not change netMassOut_col or out_flux_col, so it is done here at the top of
      ! the k-loop of step B, which recovers the answers from doing the two steps in turn.
      do k=1,nz
        h_k = h(i,j,k)
        Tr_k = Tr(i,j,k)

        if (k == 1) then
          ! A/ Update tracer due to incoming mass flux.

          ! Change in state due to forcing
          dThickness = fluxes%netMassIn(i,j) ! Since we are adding mass, we can use all of it
          dTracer = 0.

          ! Update the forcing by the part to be consumed within the present k-layer.
          dTracer = dTracer + in_flux_col
          in_flux_col = 0.0

          ! Update state
          hOld = h_k                ! Keep original thickness in hand
          h_k  = h_k + dThickness   ! New thickness
          if (h_k > 0.0) then
            Ithickness  = 1.0/h_k   ! Inverse new thickness
            ! The "if"s below avoid changing T/S by roundoff unnecessarily
            if (dThickness /= 0. .or. dTracer /= 0.) Tr_k = (hOld*Tr_k + dTracer)*Ithickness
          endif
        endif

        ! B/ Update tracer from mass leaving ocean

        ! Place forcing into this layer if this layer has nontrivial thickness.
        ! For layers thin relative to 1/IforcingDepthScale, then distribute
        ! forcing into deeper layers.
        IforcingDepthScale = 1. / max(GV%H_subroundoff, minimum_forcing_depth - netMassOut_col )
        ! fractionOfForcing = 1.0, unless h_k is less than IforcingDepthScale.
        fractionOfForcing = min(1.0, h_k*IforcingDepthScale)

        ! In the case with (-1)*netMassOut*fractionOfForcing greater than cfl*h, we
        ! limit the forcing applied to this cell, leaving the remaining forcing to
        ! be distributed downwards.
        if (-fractionOfForcing*netMassOut_col > evap_CFL_limit*h_k) then
          fractionOfForcing = -evap_CFL_limit*h_k/netMassOut_col
        endif

        ! Change in state due to forcing
        dThickness = max( fractionOfForcing*netMassOut_col, -h_k )
        ! Note this is slightly different to how salt is currently treated
        dTracer =  fractionOfForcing*out_flux_col

        ! Update the forcing by the part to be consumed within the present k-layer.
        ! If fractionOfForcing = 1, then new netMassOut vanishes.
        netMassOut_col = netMassOut_col - dThickness
        out_flux_col = out_flux_col - dTracer

        ! Update state by the appropriate increment.
        hOld = h_k                ! Keep original thickness in hand
        h_k  = h_k + dThickness   ! New thickness
        if (h_k > 0.) then
          Ithickness = 1.0/h_k    ! Inverse of new thickness
          Tr_k = (hOld*Tr_k + dTracer)*Ithickness
        endif

        ! C/ Copy the updated tracer concentration and thickness back into the model state.
        Tr(i,j,k) = Tr_k
        if (update_h) h(i,j,k) = h_k
      enddo ! k
    endif

    !   If anything remains after the k-loop, then we have grounded out, which is a problem.  This
    ! is recorded for land points too, and the locations are gathered on the host below.
    hGround(i,j) = abs(in_flux_col) + abs(out_flux_col)
    if (hGround(i,j) /= 0.0) numberOfGroundings = numberOfGroundings + 1
  enddo

  if (numberOfGroundings>0) then
    !   Record the locations of the first few groundings in the order in which a serial i-j loop
    ! would find them, for the warning messages.
    !$omp target update from(hGround)
    nGround = 0
    do j=js,je ; do i=is,ie ; if (hGround(i,j) /= 0.0) then
      nGround = nGround + 1
      if (nGround<=maxGroundings) then
        iGround(nGround) = i ! Record i,j location of event for
        jGround(nGround) = j ! warning message
        hGrounding(nGround) = hGround(i,j)
      endif
    endif ; enddo ; enddo

    do i = 1, min(numberOfGroundings, maxGroundings)
      write(mesg(1:45),'(3es15.3)') G%geoLonT( iGround(i), jGround(i) ), &
                             G%geoLatT( iGround(i), jGround(i)) , hGrounding(i)
      call MOM_error(WARNING, "MOM_tracer_diabatic.F90, applyTracerBoundaryFluxesInOut(): "//&
                              "Tracer created. x,y,dh= "//trim(mesg), all_print=.true.)
    enddo

    if (numberOfGroundings - maxGroundings > 0) then
      write(mesg, '(I0)') numberOfGroundings - maxGroundings
      call MOM_error(WARNING, "MOM_tracer_vertical.F90, applyTracerBoundaryFluxesInOut(): "//&
                              trim(mesg) // " groundings remaining", all_print=.true.)
    endif
  endif

  !$omp target exit data map(release: out_flux_optional) if(present(out_flux_optional))
  !$omp target exit data map(release: in_flux_optional) if(present(in_flux_optional))
  !$omp target exit data map(delete: in_flux, out_flux, hGround)

end subroutine applyTracerBoundaryFluxesInOut
end module MOM_tracer_diabatic
