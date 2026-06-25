! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Calculations of isoneutral slopes and stratification.
module MOM_isopycnal_slopes

use MOM_debugging,     only : hchksum, uvchksum
use MOM_error_handler, only : MOM_error, FATAL
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : calculate_density_derivs, calculate_density_second_derivs, EOS_domain
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S

implicit none ; private

#include <MOM_memory.h>
#include "do_concurrent_compat.h"

public calc_isoneutral_slopes, vert_fill_TS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Calculate isopycnal slopes, and optionally return other stratification dependent functions such as N^2
!! and dz*S^2*g-prime used, or calculable from factors used, during the calculation.
subroutine calc_isoneutral_slopes(G, GV, US, h, e, tv, dt_kappa_smooth, use_stanley, slope_x, slope_y, &
                                  niblock, njblock, nkblock, N2_u, N2_v, dzu, dzv, dzSxN, dzSyN, halo, &
                                  OBC, OBC_N2)
  integer,                                     intent(in)    :: niblock !< Blocksize in i direction
  integer,                                     intent(in)    :: njblock !< Blocksize in j direction
  integer,                                     intent(in)    :: nkblock !< Blocksize in k direction
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: e    !< Interface heights [Z ~> m]
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to various
                                                                     !! thermodynamic variables
  real,                                        intent(in)    :: dt_kappa_smooth !< A smoothing vertical
                                                                     !! diffusivity times a smoothing
                                                                     !! timescale [H Z ~> m2 or kg m-1]
  logical,                                     intent(in)    :: use_stanley !< turn on stanley param in slope
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: slope_x !< Isopycnal slope in i-dir [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: slope_y !< Isopycnal slope in j-dir [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: N2_u !< Brunt-Vaisala frequency squared at
                                                                     !! interfaces between u-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: N2_v !< Brunt-Vaisala frequency squared at
                                                                     !! interfaces between v-points [L2 Z-2 T-2 ~> s-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzu  !< Z-thickness at u-points [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzv  !< Z-thickness at v-points [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzSxN !< Z-thickness times zonal slope contribution to
                                                                     !! Eady growth rate at u-points. [Z T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                                     optional, intent(inout) :: dzSyN !< Z-thickness times meridional slope contrib. to
                                                                     !! Eady growth rate at v-points. [Z T-1 ~> m s-1]
  integer,                           optional, intent(in)    :: halo !< Halo width over which to compute
  type(ocean_OBC_type),              optional, pointer       :: OBC  !< Open boundaries control structure.
  logical,                           optional, intent(in)    :: OBC_N2 !< If present and true, use interior data
                                                                     !! to calculate stratification at open boundary
                                                                     !! condition faces.

  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: &
    T, &          ! The temperature [C ~> degC], with the values in
                  ! in massless layers filled vertically by diffusion.
    S             ! The filled salinity [S ~> ppt], with the values in
                  ! in massless layers filled vertically by diffusion.
  real, dimension(SZI_(G), SZJ_(G),SZK_(GV)+1) :: &
    pres          ! The pressure at an interface [R L2 T-2 ~> Pa].

  ! Tiled work arrays. First and second dimensions use niblock+1 / njblock+1 so that
  ! the Stanley h-point fills (which need one extra column/row beyond the u/v tile extent)
  ! fit without a separate buffer. drho_dT and drho_dS only ever need niblock x njblock,
  ! but are dimensioned +1 for alignment with the rest.
  real, dimension(niblock+1, njblock+1, nkblock) :: &
    drho_dT, &      ! The derivative of density with temperature [R C-1 ~> kg m-3 degC-1].
    drho_dS, &      ! The derivative of density with salinity [R S-1 ~> kg m-3 ppt-1].
    drho_dT_dT_h, & ! The second derivative of density with temperature at h points [R C-2 ~> kg m-3 degC-2]
    T_uvh, &        ! Temperature at the u, v or h-points for derivative calculations [C ~> degC].
    S_uvh, &        ! Salinity at the u, v or h-points for derivative calculations [S ~> ppt].
    GxSpV_uvh, &    ! g * specific volume at u/v-point interfaces [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]
    pres_uvh        ! Pressure at the u, v or h-points [R L2 T-2 ~> Pa].
  real, dimension(niblock+1) :: scrap ! Ignored output from calculate_density_second_derivs()

  real :: drdiA, drdiB  ! Along layer zonal potential density  gradients in the layers above (A)
                        ! and below (B) the interface times the grid spacing [R ~> kg m-3].
  real :: drdjA, drdjB  ! Along layer meridional potential density  gradients in the layers above (A)
                        ! and below (B) the interface times the grid spacing [R ~> kg m-3].
  real :: drdkL, drdkR  ! Vertical density differences across an interface [R ~> kg m-3].
  real :: hg2A, hg2B    ! Squares of geometric mean thicknesses [H2 ~> m2 or kg2 m-4].
  real :: hg2L, hg2R    ! Squares of geometric mean thicknesses [H2 ~> m2 or kg2 m-4].
  real :: haA, haB, haL, haR  ! Arithmetic mean thicknesses [H ~> m or kg m-2].
  real :: dzaL, dzaR    ! Temporary thicknesses in eta units [Z ~> m].
  real :: wtA, wtB      ! Unnormalized weights of the slopes above and below [H3 ~> m3 or kg3 m-6]
  real :: wtL, wtR      ! Unnormalized weights of the slopes to the left and right [H3 Z ~> m4 or kg3 m-5]
  real :: drdx, drdy    ! Zonal and meridional density gradients [R L-1 ~> kg m-4].
  real :: drdz          ! Vertical density gradient [R Z-1 ~> kg m-4].
  real :: slope         ! The slope of density surfaces, calculated in a way
                        ! that is always between -1 and 1. [Z L-1 ~> nondim]
  real :: mag_grad2     ! The squared magnitude of the 3-d density gradient [R2 Z-2 ~> kg2 m-8].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: h_neglect2    ! h_neglect^2 [H2 ~> m2 or kg2 m-4].
  real :: dz_neglect    ! A change in interface heights that is so small it is usually lost
                        ! in roundoff and can be neglected [Z ~> m].
  logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.
  real :: G_Rho0        ! The gravitational acceleration divided by density [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]

  logical :: present_N2_u, present_N2_v
  logical :: local_open_u_BC, local_open_v_BC ! True if u- or v-face OBCs exist anywhere in the global domain.
  logical :: OBC_friendly  ! If true, open boundary conditions are in use and only interior data should
                        ! be used to calculate N2 at OBC faces.
  integer, dimension(2) :: EOSdom_h1 ! The shifted i-computational domain to use for equation of
                                     ! state calculations at h points with 1 extra halo point
  integer, dimension(2) :: EOSdom_tile   !< 1-based EOS domain for the current tile [nondim].
  integer, dimension(2) :: EOSdom_tile_h1 !< 1-based EOS domain for h-point fills (one extra column) [nondim].
  integer :: is, ie, js, je, nz, IsdB
  integer :: i, j, k
  integer :: istart, iend !< First and last global i (or I) indices of the current tile [nondim].
  integer :: jstart, jend !< First and last global j (or J) indices of the current tile [nondim].
  integer :: kstart, kend !< First and last global K (or K) indices of the current tile [nondim].
  integer :: ii, jj, kk   !< Tile-local 1-based i, j, and k indices [nondim].

  ! Allocate locals on device
  !$omp target enter data map(alloc: T, S, pres, T_uvh, S_uvh, pres_uvh, GxSpV_uvh, &
  !$omp & scrap, drho_dT, drho_dS, drho_dT_dT_h)

  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
    EOSdom_h1(:) = EOS_domain(G%HI, halo=halo+1)
  else
    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
    EOSdom_h1(:) = EOS_domain(G%HI, halo=1)
  endif

  nz = GV%ke ; IsdB = G%IsdB

  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect**2
  dz_neglect = GV%dZ_subroundoff

  local_open_u_BC = .false.
  local_open_v_BC = .false.
  OBC_friendly = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    local_open_u_BC = OBC%open_u_BCs_exist_globally
    local_open_v_BC = OBC%open_v_BCs_exist_globally
    if (present(OBC_N2)) OBC_friendly = OBC_N2
  endif ; endif

  use_EOS = associated(tv%eqn_of_state)

  present_N2_u = PRESENT(N2_u)
  present_N2_v = PRESENT(N2_v)
  G_Rho0 = GV%g_Earth / GV%Rho0
  if (present_N2_u) then
    do concurrent (j=js:je, I=is-1:ie)
      N2_u(I,j,1) = 0.
      N2_u(I,j,nz+1) = 0.
    enddo
  endif
  if (present_N2_v) then
    do concurrent (J=js-1:je, i=is:ie)
      N2_v(i,J,1) = 0.
      N2_v(i,J,nz+1) = 0.
    enddo
  endif
  if (present(dzu)) then
    do concurrent (j=js:je, I=is-1:ie)
      dzu(I,j,1) = 0.
      dzu(I,j,nz+1) = 0.
    enddo
  endif
  if (present(dzv)) then
    do concurrent (J=js-1:je, i=is:ie)
      dzv(i,J,1) = 0.
      dzv(i,J,nz+1) = 0.
    enddo
  endif
  if (present(dzSxN)) then
    do concurrent (j=js:je, I=is-1:ie)
      dzSxN(I,j,1) = 0.
      dzSxN(I,j,nz+1) = 0.
    enddo
  endif
  if (present(dzSyN)) then
    do concurrent (J=js-1:je, i=is:ie)
      dzSyN(i,J,1) = 0.
      dzSyN(i,J,nz+1) = 0.
    enddo
  endif

  if (use_EOS) then
    if (present(halo)) then
      call vert_fill_TS(h, tv%T, tv%S, dt_kappa_smooth, T, S, G, GV, US, halo+1)
    else
      call vert_fill_TS(h, tv%T, tv%S, dt_kappa_smooth, T, S, G, GV, US, 1)
    endif
  endif

  if ((use_EOS .and. allocated(tv%SpV_avg) .and. (tv%valid_SpV_halo < 1)) .and. &
      (present_N2_u .or. present(dzSxN) .or. present_N2_v .or. present(dzSyN))) then
    if (tv%valid_SpV_halo < 0) then
      call MOM_error(FATAL, "calc_isoneutral_slopes called in fully non-Boussinesq mode "//&
                            "with invalid values of SpV_avg.")
    else
      call MOM_error(FATAL, "calc_isoneutral_slopes called in fully non-Boussinesq mode "//&
                            "with insufficiently large SpV_avg halos of width 0 but 1 is needed.")
    endif
  endif

  ! Find the maximum and minimum permitted streamfunction.
  if (associated(tv%p_surf)) then
    do concurrent (j=js-1:je+1, i=is-1:ie+1)
      pres(i,j,1) = tv%p_surf(i,j)
    enddo
  else
    do concurrent (j=js-1:je+1, i=is-1:ie+1)
      pres(i,j,1) = 0.0
    enddo
  endif
  ! UMW NOTE: K must be serial
  do concurrent( j=js-1:je+1 )
    do k=1,nz
      do concurrent( i=is-1:ie+1 )
        pres(i,j,K+1) = pres(i,j,K) + GV%g_Earth * GV%H_to_RZ * h(i,j,k)
      enddo
    enddo
  enddo

  ! ============================================================
  ! Zonal tiling loop: compute zonal isopycnal slopes
  ! ============================================================
  ! Outer loop tiles the (I=is-1:ie, j=js:je) domain into niblock x njblock patches.
  ! Within each tile:
  !   1. Fill tile-sized T_uvh/S_uvh/pres_uvh at u-points on device (do concurrent).
  !   2. Transfer tile arrays to host (target update from) for CPU-only EOS calls.
  !   3. Apply OBC overrides (CPU, host-side OBC derived types).
  !   4. Pre-fill GxSpV_uvh at u-points.
  !   5. Call calculate_density_derivs over tile.
  !   6. If use_stanley: refill at h-points, call calculate_density_second_derivs.
  !   7. Transfer derivatives to device (target update to).
  !   8. Compute slopes in do concurrent using tile-local indices.
  do jstart=js, je, njblock ; do istart=is-1, ie, niblock ; do kstart=2,nz, nkblock
    iend = min(istart + niblock - 1, ie)
    jend = min(jstart + njblock - 1, je)
    kend = min(kstart + nkblock - 1, nz)
    EOSdom_tile(1) = 1 ; EOSdom_tile(2) = iend - istart + 1

    ! Initialise GxSpV tile to the Boussinesq default G/rho0; overwritten below if
    ! non-Boussinesq SpV_avg is available.
    do concurrent( jj=1:jend-jstart+1, kk=1:kend-kstart+1, ii=1:iend-istart+1 )
      GxSpV_uvh(ii,jj,kk) = G_Rho0
    enddo

    if (use_EOS) then
      ! Fill tile T_uvh/S_uvh/pres_uvh at u-points (u-face = average of i and i+1)
      do concurrent(jj=1:jend-jstart+1, kk=1:kend-kstart+1, ii=1:iend-istart+1)
        i = istart + ii - 1  ! global I (u-face index)
        j = jstart + jj - 1  ! global j
        k = kstart + kk - 1  ! global k
        pres_uvh(ii,jj,kk) = 0.5*(pres(i,j,K) + pres(i+1,j,K))
        T_uvh(ii,jj,kk) = 0.25*((T(i,j,k) + T(i+1,j,k)) + (T(i,j,k-1) + T(i+1,j,k-1)))
        S_uvh(ii,jj,kk) = 0.25*((S(i,j,k) + S(i+1,j,k)) + (S(i,j,k-1) + S(i+1,j,k-1)))
      enddo

      ! Apply OBC overrides: at open boundary faces use only the interior cell's T/S/P
      ! rather than the two-cell average computed above. Runs on CPU because it requires
      ! access to OBC%segnum_u, a derived-type component not mapped to the device.
      ! UMW TODO: fix the indices in these two loops
      if (OBC_friendly) then
        ! East-facing open boundaries: interior cell is at i (=istart+ii-1), use pres/T/S(i,j,K)
        if (OBC%u_E_OBCs_on_PE) then
          do j = max(jstart, OBC%js_u_E_obc), min(jend, OBC%je_u_E_obc)
            jj = j - jstart + 1
            do K = nz, 2, -1
              do i = max(istart, OBC%Is_u_E_obc), min(iend, OBC%Ie_u_E_obc)
                ii = i - istart + 1
                if (OBC%segnum_u(i,j) > 0) then ! OBC_DIRECTION_E
                  pres_uvh(ii,jj,K) = pres(i,j,K)
                  T_uvh(ii,jj,K) = 0.5*(T(i,j,K) + T(i,j,K-1))
                  S_uvh(ii,jj,K) = 0.5*(S(i,j,K) + S(i,j,K-1))
                endif
              enddo
            enddo
          enddo
        endif
        ! West-facing open boundaries: interior cell is at i+1, use pres/T/S(i+1,j,K)
        if (OBC%u_W_OBCs_on_PE) then
          do j = max(jstart, OBC%js_u_W_obc), min(jend, OBC%je_u_W_obc)
            jj = j - jstart + 1
            do K = nz, 2, -1
              do i = max(istart, OBC%Is_u_W_obc), min(iend, OBC%Ie_u_W_obc)
                ii = i - istart + 1
                if (OBC%segnum_u(i,j) < 0) then ! OBC_DIRECTION_W
                  pres_uvh(ii,jj,K) = pres(i+1,j,K)
                  T_uvh(ii,jj,K) = 0.5*(T(i+1,j,K) + T(i+1,j,K-1))
                  S_uvh(ii,jj,K) = 0.5*(S(i+1,j,K) + S(i+1,j,K-1))
                endif
              enddo
            enddo
          enddo
        endif
      endif

      ! Pre-fill GxSpV at u-points: four-cell SpV average over the u-face and the two
      ! adjacent layers. Individual OBC faces may override this inside the slope loop.
      if (present_N2_u .or. present(dzSxN)) then
        if (allocated(tv%SpV_avg)) then
          do jj=1, jend-jstart+1 ; do kk=1,kend-kstart+1 ; do ii=1, iend-istart+1
            i = istart + ii - 1
            j = jstart + jj - 1
            k = kstart + kk - 1
            GxSpV_uvh(ii,jj,kk) = GV%g_Earth * 0.25 * ((tv%SpV_avg(i,j,K) + tv%SpV_avg(i+1,j,K)) + &
                                                        (tv%SpV_avg(i,j,K-1) + tv%SpV_avg(i+1,j,K-1)))
          enddo ; enddo ; enddo
        endif
      endif

      ! Density derivatives at u-points. EOSdom_tile is 1-based so no IsdB offset is needed.
      do jj=1,jend-jstart+1 ; do kk=1,kend-kstart+1
        call calculate_density_derivs(T_uvh(:,jj,kk), S_uvh(:,jj,kk), pres_uvh(:,jj,kk), &
                                      drho_dT(:,jj,kk), drho_dS(:,jj,kk), &
                                      tv%eqn_of_state, EOSdom_tile)
      enddo ; enddo

      if (use_stanley) then
        ! Refill T_uvh/S_uvh/pres_uvh at h-points for the Stanley second-derivative calculation.
        ! The fill extends one extra column (ii up to iend-istart+2) because the slope compute
        ! accesses drho_dT_dT_h(ii,jj,K) and drho_dT_dT_h(ii+1,jj,K) simultaneously.
        EOSdom_tile_h1(1) = 1 ; EOSdom_tile_h1(2) = iend - istart + 2
        do concurrent(jj=1:jend-jstart+1, kk=1:kend-kstart+1, ii=1:iend-istart+2)
          i = istart + ii - 1
          j = jstart + jj - 1
          k = kstart + kk - 1
          pres_uvh(ii,jj,kk) = pres(i,j,K)
          T_uvh(ii,jj,kk) = 0.5*(T(i,j,K) + T(i,j,K-1))
          S_uvh(ii,jj,kk) = 0.5*(S(i,j,K) + S(i,j,K-1))
        enddo

        do jj=1,jend-jstart+1 ; do kk=1,kend-kstart+1
          call calculate_density_second_derivs(T_uvh(:,jj,kk), S_uvh(:,jj,kk), pres_uvh(:,jj,kk), &
                     scrap, scrap, drho_dT_dT_h(:,jj,kk), scrap, scrap, &
                     tv%eqn_of_state, dom=EOSdom_tile_h1)
        enddo ; enddo
      endif ! end use_stanley

    endif ! end use_EOS for zonal tile

    ! Zonal slope compute over the tile
    do concurrent( jj=1:jend-jstart+1, kk=1:kend-kstart+1, ii=1:iend-istart+1 ) &
        local(drdkL, drdkR, drdiA, drdiB, I, j)
      I = istart + ii - 1  ! global I (u-face)
      j = jstart + jj - 1  ! global j (h-point)
      k = kstart + kk - 1  ! global k

      if (use_EOS) then
        drdiA = drho_dT(ii,jj,kk) * (T(i+1,j,k-1)-T(i,j,k-1)) + &
                drho_dS(ii,jj,kk) * (S(i+1,j,k-1)-S(i,j,k-1))
        drdiB = drho_dT(ii,jj,kk) * (T(i+1,j,k)-T(i,j,k)) + &
                drho_dS(ii,jj,kk) * (S(i+1,j,k)-S(i,j,k))
        drdkL = (drho_dT(ii,jj,kk) * (T(i,j,k)-T(i,j,k-1)) + &
                  drho_dS(ii,jj,kk) * (S(i,j,k)-S(i,j,k-1)))
        drdkR = (drho_dT(ii,jj,kk) * (T(i+1,j,k)-T(i+1,j,k-1)) + &
                  drho_dS(ii,jj,kk) * (S(i+1,j,k)-S(i+1,j,k-1)))
        if (use_stanley) then
          drdiA = drdiA + 0.5 * ((drho_dT_dT_h(ii+1,jj,kk) * tv%varT(i+1,j,K-1)) - &
                                (drho_dT_dT_h(ii,jj,kk) * tv%varT(i,j,K-1)) )
          drdiB = drdiB + 0.5 * ((drho_dT_dT_h(ii+1,jj,kk) * tv%varT(i+1,j,K)) - &
                                (drho_dT_dT_h(ii,jj,kk) * tv%varT(i,j,K)) )
        endif
      else ! .not. use_EOS: layers are constant density
        drdiA = 0.0 ; drdiB = 0.0
        drdkL = GV%Rlay(K)-GV%Rlay(K-1) ; drdkR = drdkL
      endif

      hg2A = h(i,j,k-1)*h(i+1,j,k-1) + h_neglect2
      hg2B = h(i,j,k)*h(i+1,j,k) + h_neglect2
      hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
      hg2R = h(i+1,j,k-1)*h(i+1,j,k) + h_neglect2
      haA = 0.5*(h(i,j,k-1) + h(i+1,j,k-1)) + h_neglect
      haB = 0.5*(h(i,j,k) + h(i+1,j,k)) + h_neglect
      haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
      haR = 0.5*(h(i+1,j,k-1) + h(i+1,j,k)) + h_neglect
      if (GV%Boussinesq) then
        dzaL = haL * GV%H_to_Z ; dzaR = haR * GV%H_to_Z
      else
        dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
        dzaR = 0.5*(e(i+1,j,K-1) - e(i+1,j,K+1)) + dz_neglect
      endif
      if (present(dzu)) dzu(I,j,K) = 0.5*( dzaL + dzaR )
      wtA = hg2A*haB ; wtB = hg2B*haA
      wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

      drdz = ((wtL * drdkL) + (wtR * drdkR)) / ((dzaL*wtL) + (dzaR*wtR))
      if (present_N2_u) then
        if (OBC_friendly) then ; if (OBC%segnum_u(I,j) /= 0) then
          if (OBC%segnum_u(I,j) > 0) then !  OBC_DIRECTION_E
            drdz = drdkL / dzaL
            if (use_EOS .and. allocated(tv%SpV_avg)) &
              GxSpV_uvh(ii,jj,k) = GV%g_Earth * 0.5 * (tv%SpV_avg(i,j,k) + tv%SpV_avg(i,j,k-1))
          elseif (OBC%segnum_u(I,j) < 0) then !  OBC_DIRECTION_W
            drdz = drdkR / dzaR
            if (use_EOS .and. allocated(tv%SpV_avg)) &
              GxSpV_uvh(ii,jj,k) = GV%g_Earth * 0.5 * (tv%SpV_avg(i+1,j,k) + tv%SpV_avg(i+1,j,k-1))
          endif
        endif ; endif
        N2_u(I,j,K) = GxSpV_uvh(ii,jj,kk) * drdz * G%mask2dCu(I,j)
      endif

      if (use_EOS) then
        drdx = ((wtA * drdiA + wtB * drdiB) / (wtA + wtB) - &
                drdz * (e(i,j,K)-e(i+1,j,K))) * G%IdxCu(I,j)
        mag_grad2 = (US%Z_to_L*drdx)**2 + drdz**2
        if (mag_grad2 > 0.0) then
          slope = drdx / sqrt(mag_grad2)
        else
          slope = 0.0
        endif
      else
        slope = (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
      endif

      if (local_open_u_BC) then
        if (OBC%segnum_u(I,j) /= 0) then
          if (OBC%segment(abs(OBC%segnum_u(I,j)))%open) slope = 0.
        endif
        slope = slope * max(G%mask2dT(i,j), G%mask2dT(i+1,j))
      endif

      slope_x(I,j,K) = slope
      if (present(dzSxN)) &
        dzSxN(I,j,K) = sqrt( GxSpV_uvh(ii,jj,kk) * max(0., (wtL * ( dzaL * drdkL )) &
                                                + (wtR * ( dzaR * drdkR ))) / (wtL + wtR) ) &
                        * abs(slope) * G%mask2dCu(I,j)
    enddo ! end zonal tile do concurrent
    enddo ; enddo ; enddo ! end zonal outer tile loops

  ! ============================================================
  ! Meridional tiling loop: compute meridional isopycnal slopes
  ! ============================================================
  ! Structure mirrors the zonal loop above, but tiling over (i=is:ie, J=js-1:je).
  ! The Stanley second-derivative fill extends to jend+1 in J because the slope
  ! compute accesses drho_dT_dT_h(ii,jj,K) and drho_dT_dT_h(ii,jj+1,K) simultaneously.
  do jstart=js-1, je, njblock ; do istart=is, ie, niblock ; do kstart=2,nz, nkblock
    iend = min(istart + niblock - 1, ie)
    jend = min(jstart + njblock - 1, je)
    kend = min(kstart + nkblock - 1, nz)
    EOSdom_tile(1) = 1 ; EOSdom_tile(2) = iend - istart + 1

    ! Re-initialise GxSpV tile: the zonal pass may have set entries here, so reset to
    ! G_Rho0 before the meridional fills.
    do concurrent( jj=1:jend-jstart+1, kk=1:kend-kstart+1, ii=1:iend-istart+1 )
      GxSpV_uvh(ii,jj,kk) = G_Rho0
    enddo

    if (use_EOS) then
      ! Fill tile T_uvh/S_uvh/pres_uvh at v-points (v-face = average of j and j+1)
      do concurrent(jj=1:jend-jstart+1, kk=1:kend-kstart+1, ii=1:iend-istart+1)
        i = istart + ii - 1  ! global i (h-point)
        j = jstart + jj - 1  ! global J (v-face index)
        k = kstart + kk - 1  ! global k
        pres_uvh(ii,jj,kk) = 0.5*(pres(i,j,K) + pres(i,j+1,K))
        T_uvh(ii,jj,kk) = 0.25*((T(i,j,K) + T(i,j+1,K)) + (T(i,j,K-1) + T(i,j+1,K-1)))
        S_uvh(ii,jj,kk) = 0.25*((S(i,j,K) + S(i,j+1,K)) + (S(i,j,K-1) + S(i,j+1,K-1)))
      enddo

      ! Apply OBC overrides at v-points: north-facing uses interior (j), south-facing uses j+1.
      ! UMW TODO: Fix indices here
      if (OBC_friendly) then
        if (OBC%v_N_OBCs_on_PE) then
          do j = max(jstart, OBC%Js_v_N_obc), min(jend, OBC%Je_v_N_obc)
            jj = j - jstart + 1
            do K = nz, 2, -1
              do i = max(istart, OBC%is_v_N_obc), min(iend, OBC%ie_v_N_obc)
                ii = i - istart + 1
                if (OBC%segnum_v(i,j) > 0) then ! OBC_DIRECTION_N
                  pres_uvh(ii,jj,K) = pres(i,j,K)
                  T_uvh(ii,jj,K) = 0.5*(T(i,j,K) + T(i,j,K-1))
                  S_uvh(ii,jj,K) = 0.5*(S(i,j,K) + S(i,j,K-1))
                endif
              enddo
            enddo
          enddo
        endif
        if (OBC%v_S_OBCs_on_PE) then
          do j = max(jstart, OBC%Js_v_S_obc), min(jend, OBC%Je_v_S_obc)
            jj = j - jstart + 1
            do K = nz, 2, -1
              do i = max(istart, OBC%is_v_S_obc), min(iend, OBC%ie_v_S_obc)
                ii = i - istart + 1
                if (OBC%segnum_v(i,j) < 0) then ! OBC_DIRECTION_S
                  pres_uvh(ii,jj,K) = pres(i,j+1,K)
                  T_uvh(ii,jj,K) = 0.5*(T(i,j+1,K) + T(i,j+1,K-1))
                  S_uvh(ii,jj,K) = 0.5*(S(i,j+1,K) + S(i,j+1,K-1))
                endif
              enddo
            enddo
          enddo
        endif
      endif

      ! Pre-fill GxSpV at v-points.
      if (present_N2_v .or. present(dzSyN)) then
        if (allocated(tv%SpV_avg)) then
          do jj=1, jend-jstart+1 ; do kk=1, kend-kstart+1 ; do ii=1, iend-istart+1
            i = istart + ii - 1
            j = jstart + jj - 1
            k = kstart + kk - 1
            GxSpV_uvh(ii,jj,kk) = GV%g_Earth * 0.25 * ((tv%SpV_avg(i,j,K) + tv%SpV_avg(i,j+1,K)) + &
                                                        (tv%SpV_avg(i,j,K-1) + tv%SpV_avg(i,j+1,K-1)))
          enddo ; enddo ; enddo
        endif
      endif

      ! Density derivatives at v-points.
      do jj=1,jend-jstart+1 ; do kk=1,kend-kstart+1
        call calculate_density_derivs(T_uvh(:,jj,kk), S_uvh(:,jj,kk), pres_uvh(:,jj,kk), &
                                      drho_dT(:,jj,kk), drho_dS(:,jj,kk), &
                                      tv%eqn_of_state, EOSdom_tile)
      enddo ; enddo

      if (use_stanley) then
        ! Refill at h-points for Stanley second derivatives. Extend to jend+1 in J because
        ! the slope compute reads drho_dT_dT_h(ii,jj+1,K) (the J+1 neighbour).
        EOSdom_tile_h1(1) = 1 ; EOSdom_tile_h1(2) = iend - istart + 1
        do concurrent(jj=1:jend-jstart+2, kk=1:kend-kstart, ii=1:iend-istart+1)
          i = istart + ii - 1
          j = jstart + jj - 1
          k = kstart + kk - 1
          pres_uvh(ii,jj,kk) = pres(i,j,K)
          T_uvh(ii,jj,kk) = 0.5*(T(i,j,K) + T(i,j,K-1))
          S_uvh(ii,jj,kk) = 0.5*(S(i,j,K) + S(i,j,K-1))
        enddo

        do jj=1,jend-jstart+2 ; do kk=1,kend-kstart+1
          call calculate_density_second_derivs(T_uvh(:,jj,kk), S_uvh(:,jj,kk), pres_uvh(:,jj,kk), &
                     scrap, scrap, drho_dT_dT_h(:,jj,kk), scrap, scrap, &
                     tv%eqn_of_state, dom=EOSdom_tile_h1)
        enddo ; enddo
      endif ! end use_stanley

    endif ! end use_EOS for meridional tile

    ! Meridional slope compute over the tile
    do concurrent( jj=1:jend-jstart+1, kk=1:kend-kstart, ii=1:iend-istart+1 ) &
        local(drdkL, drdkR, drdjA, drdjB, i, J)
      i = istart + ii - 1  ! global i (h-point)
      J = jstart + jj - 1  ! global J (v-face)
      K = kstart + kk - 1  ! global k

      if (use_EOS) then
        drdjA = drho_dT(ii,jj,kk) * (T(i,j+1,k-1)-T(i,j,k-1)) + &
                drho_dS(ii,jj,kk) * (S(i,j+1,k-1)-S(i,j,k-1))
        drdjB = drho_dT(ii,jj,kk) * (T(i,j+1,k)-T(i,j,k)) + &
                drho_dS(ii,jj,kk) * (S(i,j+1,k)-S(i,j,k))
        drdkL = (drho_dT(ii,jj,kk) * (T(i,j,K)-T(i,j,K-1)) + &
                  drho_dS(ii,jj,kk) * (S(i,j,K)-S(i,j,K-1)))
        drdkR = (drho_dT(ii,jj,kk) * (T(i,j+1,K)-T(i,j+1,K-1)) + &
                  drho_dS(ii,jj,kk) * (S(i,j+1,K)-S(i,j+1,K-1)))
        if (use_stanley) then
          drdjA = drdjA + 0.5 * ((drho_dT_dT_h(ii,jj+1,kk) * tv%varT(i,j+1,K-1)) - &
                                (drho_dT_dT_h(ii,jj,kk) * tv%varT(i,j,K-1)) )
          drdjB = drdjB + 0.5 * ((drho_dT_dT_h(ii,jj+1,kk) * tv%varT(i,j+1,K)) - &
                                (drho_dT_dT_h(ii,jj,kk) * tv%varT(i,j,K)) )
        endif
      else
        drdjA = 0.0 ; drdjB = 0.0
        drdkL = GV%Rlay(K)-GV%Rlay(K-1) ; drdkR = drdkL
      endif

      hg2A = h(i,j,k-1)*h(i,j+1,k-1) + h_neglect2
      hg2B = h(i,j,k)*h(i,j+1,k) + h_neglect2
      hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
      hg2R = h(i,j+1,k-1)*h(i,j+1,k) + h_neglect2
      haA = 0.5*(h(i,j,k-1) + h(i,j+1,k-1)) + h_neglect
      haB = 0.5*(h(i,j,k) + h(i,j+1,k)) + h_neglect
      haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
      haR = 0.5*(h(i,j+1,k-1) + h(i,j+1,k)) + h_neglect
      if (GV%Boussinesq) then
        dzaL = haL * GV%H_to_Z ; dzaR = haR * GV%H_to_Z
      else
        dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
        dzaR = 0.5*(e(i,j+1,K-1) - e(i,j+1,K+1)) + dz_neglect
      endif
      if (present(dzv)) dzv(i,J,K) = 0.5*( dzaL + dzaR )
      wtA = hg2A*haB ; wtB = hg2B*haA
      wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

      drdz = ((wtL * drdkL) + (wtR * drdkR)) / ((dzaL*wtL) + (dzaR*wtR))
      if (present_N2_v) then
        if (OBC_friendly) then ; if (OBC%segnum_v(i,J) /= 0) then
          if (OBC%segnum_v(i,J) > 0) then !  OBC_DIRECTION_N
            drdz = drdkL / dzaL
            if (use_EOS .and. allocated(tv%SpV_avg)) &
              GxSpV_uvh(ii,jj,kk) = GV%g_Earth * 0.5 * (tv%SpV_avg(i,j,k) + tv%SpV_avg(i,j,k-1))
          elseif (OBC%segnum_v(i,J) < 0) then !  OBC_DIRECTION_S
            drdz = drdkL / dzaL
            if (use_EOS .and. allocated(tv%SpV_avg)) &
              GxSpV_uvh(ii,jj,kk) = GV%g_Earth * 0.5 * (tv%SpV_avg(i,j+1,k) + tv%SpV_avg(i,j+1,k-1))
          endif
        endif ; endif
        N2_v(i,J,K) = GxSpV_uvh(ii,jj,kk) * drdz * G%mask2dCv(i,J)
      endif

      if (use_EOS) then
        drdy = ((wtA * drdjA + wtB * drdjB) / (wtA + wtB) - &
                drdz * (e(i,j,K)-e(i,j+1,K))) * G%IdyCv(i,J)
        mag_grad2 = (US%Z_to_L*drdy)**2 + drdz**2
        if (mag_grad2 > 0.0) then
          slope = drdy / sqrt(mag_grad2)
        else
          slope = 0.0
        endif
      else
        slope = (e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J)
      endif

      if (local_open_v_BC) then
        if (OBC%segnum_v(i,J) /= 0) then
          if (OBC%segment(abs(OBC%segnum_v(i,J)))%open) slope = 0.
        endif
        slope = slope * max(G%mask2dT(i,j), G%mask2dT(i,j+1))
      endif
      slope_y(i,J,K) = slope
      if (present(dzSyN)) &
        dzSyN(i,J,K) = sqrt( GxSpV_uvh(ii,jj,kk) * max(0., (wtL * ( dzaL * drdkL )) &
                                                + (wtR * ( dzaR * drdkR ))) / (wtL + wtR) ) &
                        * abs(slope) * G%mask2dCv(i,J)
    enddo ! end meridional tile do concurrent
  enddo ; enddo ; enddo! end meridional outer tile loops

  ! Delete all tile and field arrays from device
  !$omp target exit data map(delete: T, S, pres, T_uvh, S_uvh, pres_uvh, GxSpV_uvh, &
  !$omp   scrap, drho_dT, drho_dS, drho_dT_dT_h)

end subroutine calc_isoneutral_slopes

!> Returns tracer arrays (nominally T and S) with massless layers filled with
!! sensible values, by diffusing vertically with a small but constant diffusivity.
subroutine vert_fill_TS(h, T_in, S_in, kappa_dt, T_f, S_f, G, GV, US, halo_here, larger_h_denom)
  type(ocean_grid_type),                     intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                     intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: T_in !< Input temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: S_in !< Input salinity [S ~> ppt]
  real,                                      intent(in)  :: kappa_dt !< A vertical diffusivity to use for smoothing
                                                                 !! times a smoothing timescale [H Z ~> m2 or kg m-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T_f  !< Filled temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S_f  !< Filled salinity [S ~> ppt]
  integer,                         optional, intent(in)  :: halo_here !< Number of halo points to work on,
                                                                 !! 0 by default
  logical,                         optional, intent(in)  :: larger_h_denom !< Present and true, add a large
                                                                 !! enough minimal thickness in the denominator of
                                                                 !! the flux calculations so that the fluxes are
                                                                 !! never so large as eliminate the transmission
                                                                 !! of information across groups of massless layers.
  ! Local variables
  real :: ent(SZI_(G),SZK_(GV)+1)  ! The diffusive entrainment (kappa*dt)/dz
                                   ! between layers in a timestep [H ~> m or kg m-2].
  real :: b1(SZI_(G))              ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  real :: d1(SZI_(G))              ! A variable used by the tridiagonal solver [nondim], d1 = 1 - c1.
  real :: c1(SZI_(G),SZK_(GV))     ! A variable used by the tridiagonal solver [nondim].
  real :: kap_dt_x2                ! The 2*kappa_dt converted to H units [H2 ~> m2 or kg2 m-4].
  real :: h_neglect                ! A negligible thickness [H ~> m or kg m-2], to allow for zero thicknesses.
  real :: h0                       ! A negligible thickness to allow for zero thickness layers without
                                   ! completely decoupling groups of layers [H ~> m or kg m-2].
                                   ! Often 0 < h_neglect << h0.
  real :: h_tr                     ! h_tr is h at tracer points with a tiny thickness
                                   ! added to ensure positive definiteness [H ~> m or kg m-2].
  integer :: i, j, k, is, ie, js, je, nz, halo

  halo=0 ; if (present(halo_here)) halo = max(halo_here,0)

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke

  h_neglect = GV%H_subroundoff
  ! The use of the fixed rescaling factor in the next line avoids an extra call to thickness_to_dz()
  ! and the use of an extra 3-d array of vertical distnaces across layers (dz).  This would be more
  ! physically consistent, but it would also be more expensive, and given that this routine applies
  ! a small (but arbitrary) amount of mixing to clean up the properties of nearly massless layers,
  ! the added expense is hard to justify.
  kap_dt_x2 = (2.0*kappa_dt) * (US%Z_to_m*GV%m_to_H) ! Usually the latter term is GV%Z_to_H.
  h0 = h_neglect
  if (present(larger_h_denom)) then
    if (larger_h_denom) h0 = 1.0e-16*sqrt(0.5*kap_dt_x2)
  endif

  !$omp target enter data map(alloc: ent, b1, d1, c1)

  if (kap_dt_x2 <= 0.0) then
    do concurrent (k=1:nz , j=js:je , i=is:ie)
      T_f(i,j,k) = T_in(i,j,k) ; S_f(i,j,k) = S_in(i,j,k)
    enddo
  else
    !$omp target teams loop private( ent, b1, d1, c1, h_tr )
    do j=js,je
      do concurrent( i=is:ie )
        ent(i,2) = kap_dt_x2 / ((h(i,j,1)+h(i,j,2)) + h0)
        h_tr = h(i,j,1) + h_neglect
        b1(i) = 1.0 / (h_tr + ent(i,2))
        d1(i) = b1(i) * h_tr
        T_f(i,j,1) = (b1(i)*h_tr)*T_in(i,j,1)
        S_f(i,j,1) = (b1(i)*h_tr)*S_in(i,j,1)
      enddo
      do k=2,nz-1 ; do concurrent( i=is:ie )
        ent(i,K+1) = kap_dt_x2 / ((h(i,j,k)+h(i,j,k+1)) + h0)
        h_tr = h(i,j,k) + h_neglect
        c1(i,k) = ent(i,K) * b1(i)
        b1(i) = 1.0 / ((h_tr + d1(i)*ent(i,K)) + ent(i,K+1))
        d1(i) = b1(i) * (h_tr + d1(i)*ent(i,K))
        T_f(i,j,k) = b1(i) * (h_tr*T_in(i,j,k) + ent(i,K)*T_f(i,j,k-1))
        S_f(i,j,k) = b1(i) * (h_tr*S_in(i,j,k) + ent(i,K)*S_f(i,j,k-1))
      enddo ; enddo
      do concurrent( i=is:ie )
        c1(i,nz) = ent(i,nz) * b1(i)
        h_tr = h(i,j,nz) + h_neglect
        b1(i) = 1.0 / (h_tr + d1(i)*ent(i,nz))
        T_f(i,j,nz) = b1(i) * (h_tr*T_in(i,j,nz) + ent(i,nz)*T_f(i,j,nz-1))
        S_f(i,j,nz) = b1(i) * (h_tr*S_in(i,j,nz) + ent(i,nz)*S_f(i,j,nz-1))
      enddo
      do k=nz-1,1,-1 ; do concurrent( i=is:ie )
        T_f(i,j,k) = T_f(i,j,k) + c1(i,k+1)*T_f(i,j,k+1)
        S_f(i,j,k) = S_f(i,j,k) + c1(i,k+1)*S_f(i,j,k+1)
      enddo ; enddo
    enddo
  endif

  !$omp target exit data map(delete: ent, b1, d1, c1)

end subroutine vert_fill_TS

end module MOM_isopycnal_slopes
