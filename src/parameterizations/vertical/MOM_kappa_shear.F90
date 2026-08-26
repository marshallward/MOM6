! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Shear-dependent mixing following Jackson et al. 2008.
module MOM_kappa_shear

use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,         only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,     only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,     only : diag_ctrl, time_type
use MOM_debugging,         only : hchksum, Bchksum
use MOM_error_handler,     only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_EOS,               only : calculate_density_derivs
use MOM_EOS,               only : calculate_density, calculate_specific_vol_derivs

implicit none ; private

#include <MOM_memory.h>
#include "do_concurrent_compat.h"

public Calculate_kappa_shear, Calc_kappa_shear_vertex, kappa_shear_init
public kappa_shear_is_used, kappa_shear_at_vertex

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds the parameters that regulate shear mixing
type, public :: Kappa_shear_CS ; private
  real    :: RiNo_crit       !< The critical shear Richardson number for
                             !! shear-entrainment [nondim]. The theoretical value is 0.25.
                             !! The values found by Jackson et al. are 0.25-0.35.
  real    :: Shearmix_rate   !< A nondimensional rate scale for shear-driven
                             !! entrainment [nondim].  The value given by Jackson et al.
                             !! is 0.085-0.089.
  real    :: FRi_curvature   !<   A constant giving the curvature of the function
                             !! of the Richardson number that relates shear to
                             !! sources in the kappa equation [nondim].
                             !! The values found by Jackson et al. are -0.97 - -0.89.
  real    :: C_N             !<   The coefficient for the decay of TKE due to
                             !! stratification (i.e. proportional to N*tke) [nondim].
                             !! The values found by Jackson et al. are 0.24-0.28.
  real    :: C_S             !<   The coefficient for the decay of TKE due to
                             !! shear (i.e. proportional to |S|*tke) [nondim].
                             !! The values found by Jackson et al. are 0.14-0.12.
  real    :: lambda          !<   The coefficient for the buoyancy length scale
                             !! in the kappa equation [nondim].
                             !! The values found by Jackson et al. are 0.82-0.81.
  real    :: lambda2_N_S     !<   The square of the ratio of the coefficients of
                             !! the buoyancy and shear scales in the diffusivity
                             !! equation, 0 to eliminate the shear scale [nondim].
  real    :: lz_rescale      !<   A coefficient to rescale the distance to the nearest
                             !! solid boundary. This adjustment is to account for
                             !! regions where 3 dimensional turbulence prevents the
                             !! growth of shear instabilities [nondim].
  real    :: TKE_bg          !<   The background level of TKE [Z2 T-2 ~> m2 s-2].
  real    :: kappa_0         !<   The background diapycnal diffusivity [H Z T-1 ~> m2 s-1 or Pa s]
  real    :: kappa_seed      !<   A moderately large seed value of diapycnal diffusivity that
                             !! is used as a starting turbulent diffusivity in the iterations
                             !! to finding an energetically constrained solution for the
                             !! shear-driven diffusivity [H Z T-1 ~> m2 s-1 or Pa s]
  real    :: kappa_trunc     !< Diffusivities smaller than this are rounded to 0 [H Z T-1 ~> m2 s-1 or Pa s]
  real    :: kappa_tol_err   !<   The fractional error in kappa that is tolerated [nondim].
  real    :: Prandtl_turb    !< Prandtl number used to convert Kd_shear into viscosity [nondim].
  integer :: nkml            !<   The number of layers in the mixed layer, as
                             !! treated in this routine.  If the pieces of the
                             !! mixed layer are not to be treated collectively,
                             !! nkml is set to 1.
  integer :: max_RiNo_it     !< The maximum number of iterations that may be used
                             !! to estimate the instantaneous shear-driven mixing.
  integer :: max_KS_it       !< The maximum number of iterations that may be used
                             !! to estimate the time-averaged diffusivity.
  logical :: dKdQ_iteration_bug !< If true. use an older, dimensionally inconsistent estimate of
                             !! the derivative of diffusivity with energy in the Newton's method
                             !! iteration.  The bug causes under-corrections when dz > 1m.
  logical :: KS_at_vertex    !< If true, do the calculations of the shear-driven mixing
                             !! at the cell vertices (i.e., the vorticity points).
  logical :: eliminate_massless !< If true, massless layers are merged with neighboring
                             !! massive layers in this calculation.
                             !  I can think of no good reason why this should be false. - RWH
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  real    :: kappa_src_max_chg !< The maximum permitted increase in the kappa source within an
                             !! iteration relative to the local source [nondim].  This must be
                             !! greater than 1.  The lower limit for the permitted fractional
                             !! decrease is (1 - 0.5/kappa_src_max_chg).  These limits could
                             !! perhaps be made dynamic with an improved iterative solver.
  real    :: VS_GeoMean_Kdmin !< A minimum diffusivity for computing the horizontal averages
                             !! when using the geometric mean with VERTEX_SHEAR=True.  The model
                             !! is sensitive to this value, which is a drawback of using the
                             !! geometric average as currently implemented.
  logical :: psurf_bug       !< If true, do a simple average of the cell surface pressures to get a
                             !! surface pressure at the corner if VERTEX_SHEAR=True.  Otherwise mask
                             !! out any land points in the average.
  logical :: all_layer_TKE_bug !< If true, report back the latest estimate of TKE instead of the
                             !! time average TKE when there is mass in all layers.  Otherwise always
                             !! report the time-averaged TKE, as is currently done when there
                             !! are some massless layers.
  logical :: VS_viscosity_bug !< If true, use a bug in the calculation of the viscosity that sets
                             !! it to zero for all vertices that are on a coastline.
  logical :: vertex_shear_OBC_bug !< If false, use extra masking when interpolating thicknesses to velocity
                             !! points for setting up the shear velocities at vertices to avoid using
                             !! external thicknesses at open boundaries.  When OBCs are not in use,
                             !! this parameter does not change answers, but true is more efficient.
  logical :: VS_GeometricMean !< If true use geometric averaging for Kd from vertices to tracer points
  logical :: VS_ThicknessMean !< If true use thickness weighting when averaging Kd from vertices to
                             !! tracer points
  logical :: restrictive_tolerance_check !< If false, uses the less restrictive tolerance check to
                             !! determine if a timestep is acceptable for the KS_it outer iteration
                             !! loop, as the code was originally written.  True uses the more
                             !! restrictive check.
!  logical :: layer_stagger = .false. ! If true, do the calculations centered at
                             !  layers, rather than the interfaces.
  integer :: niblock         !< The i block size used in the kappa shear calculations [nondim].
  integer :: njblock         !< The j block size used in the kappa shear calculations [nondim].
  logical :: debug = .false. !< If true, write verbose debugging messages.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_Kd_shear = -1, id_TKE = -1, id_Kd_vertex = -1, &
             id_S2_init = -1, id_N2_init = -1, id_S2_mean = -1, id_N2_mean = -1
  !>@}
end type Kappa_shear_CS

! integer :: id_clock_project, id_clock_KQ, id_clock_avg, id_clock_setup

!> A compile-time ceiling on the number of layers in GPU builds, used to give the
!! device-executed column routines fixed-size (stack) local arrays instead of
!! per-call device-heap automatic allocations, which exhaust the default device
!! heap and serialize on the device allocator.  Checked against GV%ke in
!! kappa_shear_init.  Unused in CPU builds, where the declarations keep their
!! exact runtime sizes.
integer, parameter :: GPU_nk_max = 128

contains

!> Subroutine for calculating shear-driven diffusivity and TKE in tracer columns
subroutine Calculate_kappa_shear(u_in, v_in, h, tv, p_surf, kappa_io, tke_io, &
                                 kv_io, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: u_in   !< Initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: v_in   !< Initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa] (or NULL).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(inout) :: kappa_io !< The diapycnal diffusivity at each interface
                                                   !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].  Initially this
                                                   !! is the value from the previous timestep, which may
                                                   !! accelerate the iteration toward convergence.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out) :: tke_io   !< The turbulent kinetic energy per unit mass at
                                                   !! each interface (not layer!) [Z2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(inout) :: kv_io  !< The vertical viscosity at each interface
                                                   !! (not layer!) [H Z T-1 ~> m2 s-1 or Pa s]. This discards any
                                                   !! previous value (i.e. it is intent out) and
                                                   !! simply sets Kv = Prandtl * Kd_shear
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  type(Kappa_shear_CS),    pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to kappa_shear_init.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    diag_N2_init, & ! Diagnostic of N2 as provided to this routine [T-2 ~> s-2]
    diag_S2_init, & ! Diagnostic of S2 as provided to this routine [T-2 ~> s-2]
    diag_N2_mean, & ! Diagnostic of N2 averaged over the timestep applied [T-2 ~> s-2]
    diag_S2_mean    ! Diagnostic of S2 averaged over the timestep applied [T-2 ~> s-2]

  real, dimension(merge(G%iec-G%isc+1, CS%niblock, CS%niblock==0), &
                  merge(G%jec-G%jsc+1, CS%njblock, CS%njblock==0), SZK_(GV)+1) :: &
    kappa, &        ! The shear-driven diapycnal diffusivity at an interface [H Z T-1 ~> m2 s-1 or Pa s]
    tke, &          ! The Turbulent Kinetic Energy per unit mass at an interface [Z2 T-2 ~> m2 s-2].
    kappa_avg, &    ! The time-weighted average of kappa [H Z T-1 ~> m2 s-1 or Pa s]
    tke_avg, &      ! The time-weighted average of TKE [Z2 T-2 ~> m2 s-2]
    N2_init, &      ! N2 as provided to this routine on the vertically reduced grid,
                    ! valid only for K=1:nzc+1 [T-2 ~> s-2].
    S2_init, &      ! S2 as provided to this routine on the vertically reduced grid,
                    ! valid only for K=1:nzc+1 [T-2 ~> s-2].
    N2_mean, &      ! The time-weighted average of N2 on the vertically reduced grid,
                    ! valid only for K=1:nzc+1 [T-2 ~> s-2].
    S2_mean, &      ! The time-weighted average of S2 on the vertically reduced grid,
                    ! valid only for K=1:nzc+1 [T-2 ~> s-2].
    kf, &           ! The fractional weight of interface kc+1 for
                    ! interpolating back to the original index space [nondim].
    u, &            ! The zonal velocity after a timestep of mixing [L T-1 ~> m s-1].
    v, &            ! The meridional velocity after a timestep of mixing [L T-1 ~> m s-1].
    T, &            ! The potential temperature after a timestep of mixing [C ~> degC].
    Sal, &          ! The salinity after a timestep of mixing [S ~> ppt].
    Pressure, &     ! The pressure after a timestep of mixing [R L2 T-2 ~> Pa].
    I_dz_int, &     ! The inverse of the distance between velocity & density points
                    ! above and below an interface [Z-1 ~> m-1].  This is used to
                    ! calculate N2, shear and fluxes
    T_int, &        ! The temperature interpolated to an interface [C ~> degC].
    Sal_int, &      ! The salinity interpolated to an interface [S ~> ppt].
    dbuoy_dT, &     ! The partial derivative of buoyancy with changes in temperature [Z T-2 C-1 ~> m s-2 degC-1]
    dbuoy_dS, &     ! The partial derivative of buoyancy with changes in salinity [Z T-2 S-1 ~> m s-2 ppt-1]
    dSpV_dT, &      ! The partial derivative of specific volume with changes in temperature [R-1 C-1 ~> m3 kg-1 degC-1]
    dSpV_dS, &      ! The partial derivative of specific volume with changes in salinity [R-1 S-1 ~> m3 kg-1 ppt-1]
    rho_int, &      ! The in situ density interpolated to an interface [R ~> kg m-3]
    kappa_full, &   ! kappa mapped back to the original interfaces [H Z T-1 ~> m2 s-1 or Pa s]
    tke_full        ! tke mapped back to the original interfaces [Z2 T-2 ~> m2 s-2].

  integer, dimension(merge(G%iec-G%isc+1, CS%niblock, CS%niblock==0), &
                     merge(G%jec-G%jsc+1, CS%njblock, CS%njblock==0), SZK_(GV)+1) :: &
    kc              ! The index map between the original
                    ! interfaces and the interfaces with massless layers
                    ! merged into nearby massive layers.
  integer, dimension(merge(G%iec-G%isc+1, CS%niblock, CS%niblock==0), &
                     merge(G%jec-G%jsc+1, CS%njblock, CS%njblock==0)) :: &
    nzc_2d             ! The number of interfaces in the column after massless layers
                    ! have been merged into nearby massive layers.

  real, dimension(SZK_(GV)+1) :: &
    a1_col, &   ! a1 is the coupling between adjacent interfaces in the TKE,
                ! velocity, and density equations [H ~> m or kg m-2]
    c1_col      ! c1 is used in the tridiagonal (and similar) solvers [nondim].
  real, dimension(merge(G%iec-G%isc+1, CS%niblock, CS%niblock==0), &
                  merge(G%jec-G%jsc+1, CS%njblock, CS%njblock==0), SZK_(GV)) :: &
    dz_3d, &    ! Vertical distance between interface heights [Z ~> m], block-local.
    Idz, &      ! The inverse of the thickness of the merged layers [H-1 ~> m2 kg-1].
    h_lay, &    ! The layer thickness [H ~> m or kg m-2]
    dz_lay, &   ! The geometric layer thickness in height units [Z ~> m]
    u0xdz, &    ! The initial zonal velocity times dz [H L T-1 ~> m2 s-1 or kg m-1 s-1]
    v0xdz, &    ! The initial meridional velocity times dz [H L T-1 ~> m2 s-1 or kg m-1 s-1]
    T0xdz, &    ! The initial temperature times thickness [C H ~> degC m or degC kg m-2] or if
                ! temperature is not a state variable, the density times thickness [R H ~> kg m-2 or kg2 m-5]
    S0xdz       ! The initial salinity times dz [S H ~> ppt m or ppt kg m-2].

  real :: f2    ! The squared Coriolis parameter of each column [T-2 ~> s-2].
  real :: surface_pres  ! The top surface pressure [R L2 T-2 ~> Pa].
  real :: b1            ! The inverse of the pivot in the tridiagonal equations [H-1 ~> m-1 or m2 kg-1].
  real :: d1            ! 1 - c1 in the tridiagonal equations [nondim]
  real :: bd1           ! A term in the denominator of b1 [H ~> m or kg m-2].
  real :: dz_in_lay     !   The running sum of the thickness in a layer [H ~> m or kg m-2]
  real :: k0dt          ! The background diffusivity times the timestep [H Z ~> m2 or kg m-1]
  real :: gR0           ! A conversion factor from H to pressure, Rho_0 times g in Boussinesq
                        ! mode, or just g when non-Boussinesq [R L2 T-2 H-1 ~> kg m-2 s-2 or m s-2].
  real :: g_R0          ! g_R0 is a rescaled version of g/Rho [Z R-1 T-2 ~> m4 kg-1 s-2].
  real :: dz_massless   ! A layer thickness that is considered massless [H ~> m or kg m-2]
  logical :: use_temperature  !  If true, temperature and salinity have been
                        ! allocated and are being used as state variables.
  integer, dimension(3,2) :: EOSdom    !< 1-based EOS domain density deriv calculations
  integer :: is, ie, js, je, i, j, k, nz, nzc
  integer :: nii, njj   ! The block sizes actually used, with 0 resolved to the full domain.
  integer :: isb, ieb   ! I-block bounds of the current block.
  integer :: jsb, jeb   ! J-block bounds of the current block.
  integer :: iie, jje   ! Valid block-local end bounds for blocked loops.
  integer :: ii, jj     ! Block-local 1-based i and j indices.
  integer :: nzc_max    ! Maximum value of nzc across domain

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  nii = merge(ie-is+1, CS%niblock, CS%niblock==0)
  njj = merge(je-js+1, CS%njblock, CS%njblock==0)

  use_temperature = associated(tv%T)

  k0dt = dt*CS%kappa_0
  gR0 = GV%H_to_RZ * GV%g_Earth
  g_R0 = GV%g_Earth_Z_T2 / GV%Rho0
  dz_massless = 0.1*sqrt((US%Z_to_m*GV%m_to_H)*k0dt)
  nzc_max = 0

  if ((CS%id_N2_init>0) .or. CS%debug) diag_N2_init(:,:,:) = 0.0
  if ((CS%id_S2_init>0) .or. CS%debug) diag_S2_init(:,:,:) = 0.0
  if (CS%id_N2_mean>0) diag_N2_mean(:,:,:) = 0.0
  if (CS%id_S2_mean>0) diag_S2_mean(:,:,:) = 0.0

  !---------------------------------------
  ! Work on each column.
  !---------------------------------------

  ! Diagnostics to allocated to avoid large transfers to and from the device.
  !$omp target enter data map(alloc: diag_N2_init, diag_S2_init, diag_N2_mean, diag_S2_mean)

  !$omp target enter data map(alloc: h_lay, dz_lay, u0xdz, v0xdz, T0xdz, S0xdz, kf, kc, kappa, Idz, &
  !$omp &                            tke, kappa_avg, tke_avg, N2_init, S2_init, N2_mean, S2_mean, &
  !$omp &                            dbuoy_dT, dbuoy_dS, dSpV_dT, dSpV_dS, rho_int, T_int, Sal_int, &
  !$omp &                            pressure, I_dz_int, u, v, T, Sal, nzc_2d, a1_col, c1_col,&
  !$omp &                            kappa_full, tke_full, dz_3d )

  do jsb=js,je,njj ; do isb=is,ie,nii
    ieb = min(isb + nii - 1, ie)
    jeb = min(jsb + njj - 1, je)
    iie = ieb - isb + 1 ; jje = jeb - jsb + 1

    ! Convert layer thicknesses into geometric thickness in height units.
    call thickness_to_dz(h, tv, dz_3d, G, GV, US, niB=nii, njB=njj, do_offload=.true., &
                         i_lo=isb, i_hi=ieb, j_lo=jsb, j_hi=jeb)

    ! Zero out arrays passed to density derivative calculation to prevent uninitialized values
    ! from entering call. Keeping this out of the per-column k-loop so it vectorizes on the CPU
    do concurrent( k=1:nz+1, jj=1:jje, ii=1:iie )
      T_int(ii,jj,k) = 0.0 ; Sal_int(ii,jj,k) = 0.0 ; pressure(ii,jj,k) = 0.0
      kappa(ii,jj,K) = CS%kappa_seed
    enddo

    nzc_max = 0

    !$omp target teams loop collapse(2) &
    !$omp & private(nzc, surface_pres, dz_in_lay, b1, d1, bd1, i, j, a1_col, c1_col) &
    !$omp & reduction(max:nzc_max)
    do  jj=1,jje ; do ii=1,iie ; if (G%mask2dT(isb+ii-1, jsb+jj-1) > 0.0 ) then
      i = isb + ii - 1 ; j = jsb + jj - 1

      ! Store a transposed version of the initial arrays.
      ! Any elimination of massless layers would occur here.
      if (CS%eliminate_massless) then
        nzc = 1
        ! This section just finds massive layers and sets the values of
        ! various variables in those layers.
        do k=1,nz
          ! Zero out the thicknesses of all layers, even if they are unused.
          h_lay(ii,jj,k) = 0.0 ; dz_lay(ii,jj,k) = 0.0 ; u0xdz(ii,jj,k) = 0.0 ; v0xdz(ii,jj,k) = 0.0
          T0xdz(ii,jj,k) = 0.0 ; S0xdz(ii,jj,k) = 0.0

          ! Add a new layer if this one has mass.
          ! if ((h_lay(nzc) > 0.0) .and. (h_1d(k) > dz_massless)) nzc = nzc+1
          if ((k>CS%nkml) .and. (h_lay(ii,jj,nzc) > 0.0) .and. (h(i,j,k) > dz_massless)) then
            nzc = nzc+1
          endif

          kc(ii,jj,k) = nzc
          h_lay(ii,jj,nzc) = h_lay(ii,jj,nzc) + h(i,j,k)
          dz_lay(ii,jj,nzc) = dz_lay(ii,jj,nzc) + dz_3d(ii,jj,k)
          u0xdz(ii,jj,nzc) = u0xdz(ii,jj,nzc) + u_in(i,j,k)*h(i,j,k)
          v0xdz(ii,jj,nzc) = v0xdz(ii,jj,nzc) + v_in(i,j,k)*h(i,j,k)
          if (use_temperature) then
            T0xdz(ii,jj,nzc) = T0xdz(ii,jj,nzc) + tv%T(i,j,k)*h(i,j,k)
            S0xdz(ii,jj,nzc) = S0xdz(ii,jj,nzc) + tv%S(i,j,k)*h(i,j,k)
          else
            T0xdz(ii,jj,nzc) = T0xdz(ii,jj,nzc) + GV%Rlay(k)*h(i,j,k)
            S0xdz(ii,jj,nzc) = S0xdz(ii,jj,nzc) + GV%Rlay(k)*h(i,j,k)
          endif
        enddo
        kc(ii,jj,nz+1) = nzc+1

        ! Set up Idz as the inverse of layer thicknesses.
        do k=1,nzc ; Idz(ii,jj,k) = 1.0 / h_lay(ii,jj,k) ; enddo

        !   Now determine kf, the fractional weight of interface kc when
        ! interpolating between interfaces kc and kc+1.
        kf(ii,jj,1) = 0.0 ; dz_in_lay = h(i,j,1)
        do k=2,nz
          if (kc(ii,jj,k) > kc(ii,jj,k-1)) then
            kf(ii,jj,k) = 0.0
            dz_in_lay = h(i,j,k)
          else
            kf(ii,jj,k) = dz_in_lay*Idz(ii,jj,kc(ii,jj,k))
            dz_in_lay = dz_in_lay + h(i,j,k)
          endif
        enddo
        kf(ii,jj,nz+1) = 0.0
      else
        do k=1,nz
          h_lay(ii,jj,k) = h(i,j,k)
          dz_lay(ii,jj,k) = dz_3d(ii,jj,k)
          u0xdz(ii,jj,k) = u_in(i,j,k)*h_lay(ii,jj,k) ; v0xdz(ii,jj,k) = v_in(i,j,k)*h_lay(ii,jj,k)
        enddo
        if (use_temperature) then
          do k=1,nz
            T0xdz(ii,jj,k) = tv%T(i,j,k)*h_lay(ii,jj,k) ; S0xdz(ii,jj,k) = tv%S(i,j,k)*h_lay(ii,jj,k)
          enddo
        else
          do k=1,nz
            T0xdz(ii,jj,k) = GV%Rlay(k)*h_lay(ii,jj,k) ; S0xdz(ii,jj,k) = GV%Rlay(k)*h_lay(ii,jj,k)
          enddo
        endif
        nzc = nz
        do k=1,nzc+1 ; kc(ii,jj,k) = k ; kf(ii,jj,k) = 0.0 ; enddo
      endif

      ! Get inverse of layer thicknesses for applying background diffisivity
      !   Set up I_dz_int as the inverse of the distance between
      ! adjacent layer centers.
      I_dz_int(ii,jj,1) = 2.0 / dz_lay(ii,jj,1)
      do K=2,nzc
        I_dz_int(ii,jj,K) = 2.0 / (dz_lay(ii,jj,K-1) + dz_lay(ii,jj,K))
      enddo
      I_dz_int(ii,jj,nzc+1) = 2.0 / dz_lay(ii,jj,nzc)

      ! Get velocities, thickness, and thermodynamic tracers after background diffusion step
      ! Determine the velocities and thicknesses after eliminating massless
      ! layers and applying a time-step of background diffusion.
      if (nzc > 1) then
        a1_col(2) = k0dt*I_dz_int(ii,jj,2)
        b1 = 1.0 / (h_lay(ii,jj,1) + a1_col(2))
        u(ii,jj,1) = b1 * u0xdz(ii,jj,1) ; v(ii,jj,1) = b1 * v0xdz(ii,jj,1)
        T(ii,jj,1) = b1 * T0xdz(ii,jj,1) ; Sal(ii,jj,1) = b1 * S0xdz(ii,jj,1)
        c1_col(2) = a1_col(2) * b1 ; d1 = h_lay(ii,jj,1) * b1 ! = 1 - c1
        do k=2,nzc-1
          bd1 = h_lay(ii,jj,k) + d1*a1_col(k)
          a1_col(k+1) = k0dt*I_dz_int(ii,jj,k+1)
          b1 = 1.0 / (bd1 + a1_col(k+1))
          u(ii,jj,k) = b1 * (u0xdz(ii,jj,k) + a1_col(k)*u(ii,jj,k-1))
          v(ii,jj,k) = b1 * (v0xdz(ii,jj,k) + a1_col(k)*v(ii,jj,k-1))
          T(ii,jj,k) = b1 * (T0xdz(ii,jj,k) + a1_col(k)*T(ii,jj,k-1))
          Sal(ii,jj,k) = b1 * (S0xdz(ii,jj,k) + a1_col(k)*Sal(ii,jj,k-1))
          c1_col(k+1) = a1_col(k+1) * b1 ; d1 = bd1 * b1 ! d1 = 1 - c1
        enddo
        ! rho or T and S have insulating boundary conditions, u & v use no-slip
        ! bottom boundary conditions (if kappa0 > 0).
        ! For no-slip bottom boundary conditions
        b1 = 1.0 / ((h_lay(ii,jj,nzc) + d1*a1_col(nzc)) + k0dt*I_dz_int(ii,jj,nzc+1))
        u(ii,jj,nzc) = b1 * (u0xdz(ii,jj,nzc) + a1_col(nzc)*u(ii,jj,nzc-1))
        v(ii,jj,nzc) = b1 * (v0xdz(ii,jj,nzc) + a1_col(nzc)*v(ii,jj,nzc-1))
        ! For insulating boundary conditions
        b1 = 1.0 / (h_lay(ii,jj,nzc) + d1*a1_col(nzc))
        T(ii,jj,nzc) = b1 * (T0xdz(ii,jj,nzc) + a1_col(nzc)*T(ii,jj,nzc-1))
        Sal(ii,jj,nzc) = b1 * (S0xdz(ii,jj,nzc) + a1_col(nzc)*Sal(ii,jj,nzc-1))
        do k=nzc-1,1,-1
          u(ii,jj,k) = u(ii,jj,k) + c1_col(k+1)*u(ii,jj,k+1) ; v(ii,jj,k) = v(ii,jj,k) + c1_col(k+1)*v(ii,jj,k+1)
          T(ii,jj,k) = T(ii,jj,k) + c1_col(k+1)*T(ii,jj,k+1) ; Sal(ii,jj,k) = Sal(ii,jj,k) + c1_col(k+1)*Sal(ii,jj,k+1)
        enddo
      else
        ! This is correct, but probably unnecessary.
        b1 = 1.0 / (h_lay(ii,jj,1) + k0dt*I_dz_int(ii,jj,2))
        u(ii,jj,1) = b1 * u0xdz(ii,jj,1) ; v(ii,jj,1) = b1 * v0xdz(ii,jj,1)
        b1 = 1.0 / h_lay(ii,jj,1)
        T(ii,jj,1) = b1 * T0xdz(ii,jj,1) ; Sal(ii,jj,1) = b1 * S0xdz(ii,jj,1)
      endif

      ! Get T,S and pressure at interfaces for calculating dbuoy_dT and dbuoy_dS
      surface_pres = 0.0 ; if (associated(p_surf)) surface_pres = p_surf(i,j)
      if (use_temperature) then
        pressure(ii,jj,1) = surface_pres
        do k=2,nzc
          pressure(ii,jj,k) = pressure(ii,jj,k-1) + gR0*h_lay(ii,jj,k-1)
          T_int(ii,jj,k) = 0.5*(T(ii,jj,k-1) + T(ii,jj,k))
          Sal_int(ii,jj,k) = 0.5*(Sal(ii,jj,k-1) + Sal(ii,jj,k))
        enddo
      endif

      ! Save the number of collapsed layers for each column
      nzc_2d(ii,jj) = nzc
      nzc_max = max( nzc_max, nzc)

    endif ; enddo ; enddo ! end of j-loop, ! end of i-loop, !end of if (G%mask2dT(i,j) > 0.0)


    ! Get dbouy_Dt and dbouy_DS
    ! Calculate thermodynamic coefficients and an initial estimate of N2.
    if (use_temperature) then
      if (GV%Boussinesq .or. GV%semi_Boussinesq) then
        EOSdom(1,1) = 1 ; EOSdom(1,2) = iie
        EOSdom(2,1) = 1 ; EOSdom(2,2) = jje
        EOSdom(3,1) = 2 ; EOSdom(3,2) = nzc_max
        call calculate_density_derivs(T_int, Sal_int, pressure, dbuoy_dT, dbuoy_dS, &
                                      tv%eqn_of_state, EOSdom, scale=-g_R0 )
      else
        ! These should perhaps be combined into a single call to calculate the thermal expansion
        ! and haline contraction coefficients?
        !$omp target update from(T_int, Sal_int, pressure, nzc_2d)
        do jj=1,jje ; do ii=1,iie
          if (G%mask2dT(isb+ii-1, jsb+jj-1) > 0.0) then
            nzc = nzc_2d(ii,jj)
            call calculate_specific_vol_derivs(T_int(ii,jj,:), Sal_int(ii,jj,:), pressure(ii,jj,:), &
                                          dSpV_dT(ii,jj,:), dSpV_dS(ii,jj,:), tv%eqn_of_state, (/2,nzc/) )
            call calculate_density(T_int(ii,jj,:), Sal_int(ii,jj,:), pressure(ii,jj,:), &
                                   rho_int(ii,jj,:), tv%eqn_of_state, (/2,nzc/) )
            do K=2,nzc
              dbuoy_dT(ii,jj,K) = GV%g_Earth_Z_T2 * (rho_int(ii,jj,K) * dSpV_dT(ii,jj,K))
              dbuoy_dS(ii,jj,K) = GV%g_Earth_Z_T2 * (rho_int(ii,jj,K) * dSpV_dS(ii,jj,K))
            enddo
          endif
        enddo ; enddo
        !$omp target update to(dbuoy_dT, dbuoy_dS)
      endif
    elseif (GV%Boussinesq .or. GV%semi_Boussinesq) then
      !$omp target update from(nzc_2d)
      do jj=1,jje ; do ii=1,iie
        if (G%mask2dT(isb+ii-1, jsb+jj-1) > 0.0) then
          nzc = nzc_2d(ii,jj)
          do K=1,nzc+1 ; dbuoy_dT(ii,jj,K) = -g_R0 ; dbuoy_dS(ii,jj,K) = 0.0 ; enddo
        endif
      enddo ; enddo
      !$omp target update to(dbuoy_dT, dbuoy_dS)
    else
      !$omp target update from(nzc_2d)
      do jj=1,jje ; do ii=1,iie
        if (G%mask2dT(isb+ii-1, jsb+jj-1) > 0.0) then
          nzc = nzc_2d(ii,jj)
          do K=1,nzc+1 ; dbuoy_dS(ii,jj,K) = 0.0 ; enddo
          dbuoy_dT(ii,jj,1) = -GV%g_Earth_Z_T2 / GV%Rlay(1)
          do K=2,nzc
            dbuoy_dT(ii,jj,K) = -GV%g_Earth_Z_T2 / (0.5*(GV%Rlay(K-1) + GV%Rlay(K)))
          enddo
          dbuoy_dT(ii,jj,nzc+1) = -GV%g_Earth_Z_T2 / GV%Rlay(nzc)
        endif
      enddo ; enddo
      !$omp target update to(dbuoy_dT, dbuoy_dS)
    endif

    do concurrent( jj=1:jje, ii=1:iie ) DO_LOCALITY(local( f2, k, i, j ))
      i = isb + ii - 1 ; j = jsb + jj - 1
      if (G%mask2dT(i,j) > 0.0) then

      f2 = 0.25 * ((G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J-1)) + &
                    (G%Coriolis2Bu(I,J-1) + G%Coriolis2Bu(I-1,J)))

      call kappa_shear_column(kappa, tke, dt, nzc_2d(ii,jj), f2, dbuoy_dT, dbuoy_dS, &
                             h_lay, dz_lay, I_dz_int, kappa_avg, u, v, T, Sal, &
                             tke_avg, N2_init, S2_init, N2_mean, S2_mean, &
                             CS, GV, US, nii, njj, ii, jj)

      ! call cpu_clock_begin(id_clock_setup)
      ! Extrapolate from the vertically reduced grid back to the original layers.
      if (nz == nzc_2d(ii,jj)) then
        do K=1,nz+1
          kappa_full(ii,jj,K) = kappa_avg(ii,jj,K)
          if (CS%all_layer_TKE_bug) then
            tke_full(ii,jj,K) = tke(ii,jj,K)
          else
            tke_full(ii,jj,K) = tke_avg(ii,jj,K)
          endif
        enddo
        if (CS%id_N2_mean>0) then ; do K=1,nz+1
          diag_N2_mean(i,j,K) = N2_mean(ii,jj,K)
        enddo ; endif
        if (CS%id_S2_mean>0) then ; do K=1,nz+1
          diag_S2_mean(i,j,K) = S2_mean(ii,jj,K)
        enddo ; endif
        if ((CS%id_N2_init>0) .or. CS%debug) then ; do K=1,nz+1
          diag_N2_init(i,j,K) = N2_init(ii,jj,K)
        enddo ; endif
        if ((CS%id_S2_init>0) .or. CS%debug) then ; do K=1,nz+1
          diag_S2_init(i,j,K) = S2_init(ii,jj,K)
        enddo ; endif
      else
        do K=1,nz+1
          if (kf(ii,jj,K) == 0.0) then
            kappa_full(ii,jj,K) = kappa_avg(ii,jj,kc(ii,jj,K))
            tke_full(ii,jj,K) = tke_avg(ii,jj,kc(ii,jj,K))
          else
            kappa_full(ii,jj,K) = (1.0-kf(ii,jj,K)) * kappa_avg(ii,jj,kc(ii,jj,K)) + kf(ii,jj,K) &
                                  * kappa_avg(ii,jj,kc(ii,jj,K)+1)
            tke_full(ii,jj,K) = (1.0-kf(ii,jj,K)) * tke_avg(ii,jj,kc(ii,jj,K)) + kf(ii,jj,K) &
                                * tke_avg(ii,jj,kc(ii,jj,K)+1)
          endif
        enddo
        do K=1,nz+1
          if (kf(ii,jj,K) == 0.0) then
            if (CS%id_N2_mean>0) diag_N2_mean(i,j,K) = N2_mean(ii,jj,kc(ii,jj,K))
            if (CS%id_S2_mean>0) diag_S2_mean(i,j,K) = S2_mean(ii,jj,kc(ii,jj,K))
            if ((CS%id_N2_init>0) .or. CS%debug) diag_N2_init(i,j,K) = N2_init(ii,jj,kc(ii,jj,K))
            if ((CS%id_S2_init>0) .or. CS%debug) diag_S2_init(i,j,K) = S2_init(ii,jj,kc(ii,jj,K))
          else
            if (CS%id_N2_mean>0) &
              diag_N2_mean(i,j,K) = (1.0-kf(ii,jj,K)) * N2_mean(ii,jj,kc(ii,jj,K)) + kf(ii,jj,K) &
                                    * N2_mean(ii,jj,kc(ii,jj,K)+1)
            if (CS%id_S2_mean>0) &
              diag_S2_mean(i,j,K) = (1.0-kf(ii,jj,K)) * S2_mean(ii,jj,kc(ii,jj,K)) + kf(ii,jj,K) &
                                    * S2_mean(ii,jj,kc(ii,jj,K)+1)
            if ((CS%id_N2_init>0) .or. CS%debug) &
              diag_N2_init(i,j,K) = (1.0-kf(ii,jj,K)) * N2_init(ii,jj,kc(ii,jj,K)) + kf(ii,jj,K) &
                                    * N2_init(ii,jj,kc(ii,jj,K)+1)
            if ((CS%id_S2_init>0) .or. CS%debug) &
              diag_S2_init(i,j,K) = (1.0-kf(ii,jj,K)) * S2_init(ii,jj,kc(ii,jj,K)) + kf(ii,jj,K) &
                                    * S2_init(ii,jj,kc(ii,jj,K)+1)
          endif
        enddo
      endif ! end of if (nz == nzc)

      ! call cpu_clock_end(id_clock_setup)
    else  ! Land points, still inside the i,j-loop.
      do K=1,nz+1
        kappa_full(ii,jj,K) = 0.0 ; tke_full(ii,jj,K) = 0.0
      enddo
      if (CS%id_N2_mean>0) then ; do K=1,nz+1
        diag_N2_mean(i,j,K) = 0.0
      enddo ; endif
      if (CS%id_S2_mean>0) then ; do K=1,nz+1
        diag_S2_mean(i,j,K) = 0.0
      enddo ; endif
      if ((CS%id_N2_init>0) .or. CS%debug) then ; do K=1,nz+1
        diag_N2_init(i,j,K) = 0.0
      enddo ; endif
      if ((CS%id_S2_init>0) .or. CS%debug) then ; do K=1,nz+1
        diag_S2_init(i,j,K) = 0.0
      enddo ; endif
    endif ; enddo ! end of j-loop, ! end of i-loop, !end of if (G%mask2dT(i,j) > 0.0)

    ! Write the results to full domain output arrays outside of ji loop to vectorize on CPU
    ! This mirrors what Calc_kappa_shear_vertex already does with kappa_full and tke_full.
    do concurrent( K=1:nz+1, jj=1:jje, ii=1:iie ) DO_LOCALITY(local(i, j))
      i = isb + ii - 1 ; j = jsb + jj - 1
      kappa_io(i,j,K) = kappa_full(ii,jj,K)
      tke_io(i,j,K) = tke_full(ii,jj,K)
      kv_io(i,j,K) = kappa_full(ii,jj,K) * CS%Prandtl_turb
    enddo

  enddo ; enddo ! end of the i- and j-block loops

  if (CS%debug) then
    !$omp target update from( diag_N2_init, diag_S2_init, kappa_io, tke_io )
    call hchksum(diag_N2_init, "kappa_shear N2_init", G%HI, unscale=US%s_to_T**2)
    call hchksum(diag_S2_init, "kappa_shear S2_init", G%HI, unscale=US%s_to_T**2)
    call hchksum(kappa_io, "kappa", G%HI, unscale=GV%HZ_T_to_m2_s)
    call hchksum(tke_io, "tke", G%HI, unscale=US%Z_to_m**2*US%s_to_T**2)
  endif

  !$omp target exit data map(delete: h_lay, dz_3d, dz_lay, u0xdz, v0xdz, T0xdz, S0xdz, kf, kc, kappa, Idz, &
  !$omp &                            tke, kappa_avg, tke_avg, N2_init, S2_init, N2_mean, S2_mean, &
  !$omp &                            diag_N2_init, diag_S2_init, diag_N2_mean, diag_S2_mean, &
  !$omp &                            dbuoy_dT, dbuoy_dS, dSpV_dT, dSpV_dS, rho_int, T_int, Sal_int, &
  !$omp &                            pressure, I_dz_int, u, v, T, Sal, a1_col, c1_col, nzc_2d, &
  !$omp &                            kappa_full, tke_full )

  if (CS%id_Kd_shear > 0) call post_data(CS%id_Kd_shear, kappa_io, CS%diag)
  if (CS%id_TKE > 0) call post_data(CS%id_TKE, tke_io, CS%diag)
  if (CS%id_N2_init > 0) call post_data(CS%id_N2_init, diag_N2_init, CS%diag)
  if (CS%id_S2_init > 0) call post_data(CS%id_S2_init, diag_S2_init, CS%diag)
  if (CS%id_N2_mean > 0) call post_data(CS%id_N2_mean, diag_N2_mean, CS%diag)
  if (CS%id_S2_mean > 0) call post_data(CS%id_S2_mean, diag_S2_mean, CS%diag)

end subroutine Calculate_kappa_shear


!> Subroutine for calculating shear-driven diffusivity and TKE in corner columns
subroutine Calc_kappa_shear_vertex(u_in, v_in, h, T_in, S_in, tv, p_surf, kappa_io, tke_io, &
                                   kv_io, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)   :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: u_in   !< Initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   &
                           intent(in)    :: v_in   !< Initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: T_in   !< Layer potential temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: S_in   !< Layer salinities [S ~> ppt]
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
                                                   !! (or NULL).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out)   :: kappa_io !< The diapycnal diffusivity at each interface
                                                   !! (not layer!) [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                           intent(out)   :: tke_io !< The turbulent kinetic energy per unit mass at
                                                   !! each interface (not layer!) [Z2 T-2 ~> m2 s-2].
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                           intent(inout) :: kv_io  !< The vertical viscosity at each interface
                                                   !! [H Z T-1 ~> m2 s-1 or Pa s].
                                                   !! The previous value is used to initialize kappa
                                                   !! in the vertex columns as Kappa = Kv/Prandtl
                                                   !! to accelerate the iteration toward convergence.
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  type(Kappa_shear_CS),    pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to kappa_shear_init.

  ! Local variables
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1) :: &
    diag_N2_init, & ! Diagnostic of N2 as provided to this routine [T-2 ~> s-2]
    diag_S2_init, & ! Diagnostic of S2 as provided to this routine [T-2 ~> s-2]
    diag_N2_mean, & ! Diagnostic of N2 averaged over the timestep applied [T-2 ~> s-2]
    diag_S2_mean, & ! Diagnostic of S2 averaged over the timestep applied [T-2 ~> s-2]
    kappa_vertex    ! Diffusivity at interfaces and vertices [H Z T-1 ~> m2 s-1 or Pa s]

  real, dimension(merge(G%iecB-(G%isc-1)+2, CS%niblock, CS%niblock==0), &
                  merge(G%jecB-(G%jsc-1)+2, CS%njblock, CS%njblock==0), SZK_(GV)) :: &
    dz_3d, &        ! Vertical distance between interface heights [Z ~> m], block-local, with the
                    ! same one-column/row overlap as h_at_u and h_at_v to support the corner reads.
    h_at_u, &       ! A mask-weighted thickness interpolated to u-points [H ~> m or kg m-2]
    h_at_v, &       ! A mask-weighted thickness interpolated to v-points [H ~> m or kg m-2]
    h_vrt, &        ! h interpolated to vertices [H ~> m or kg m-2].
    dz_vrt, &       ! Vertical distance between interface heights [Z ~> m].
    u_vrt, v_vrt, & ! u_in and v_in interpolated to vertices [L T-1 ~> m s-1].
    T_vrt, S_vrt, rho_vrt, & ! T [C ~> degC], S [S ~> ppt], and rho [R ~> kg m-3] at vertices.
    Idz, &      ! The inverse of the thickness of the merged layers [H-1 ~> m2 kg-1].
    h_lay, &    ! The layer thickness [H ~> m or kg m-2]
    dz_lay, &   ! The geometric layer thickness in height units [Z ~> m]
    u0xdz, &    ! The initial zonal velocity times dz [L H T-1 ~> m2 s-1 or kg m-1 s-1].
    v0xdz, &    ! The initial meridional velocity times dz [H L T-1 ~> m2 s-1 or kg m-1 s-1]
    T0xdz, &    ! The initial temperature times dz [C H ~> degC m or degC kg m-2]
    S0xdz       ! The initial salinity times dz [S H ~> ppt m or ppt kg m-2]
  real, dimension(merge(G%iecB-(G%isc-1)+2, CS%niblock, CS%niblock==0), &
                  merge(G%jecB-(G%jsc-1)+2, CS%njblock, CS%njblock==0), SZK_(GV)+1) :: &
    kappa, &    ! The shear-driven diapycnal diffusivity at an interface [H Z T-1 ~> m2 s-1 or Pa s]
    tke, &      ! The Turbulent Kinetic Energy per unit mass at an interface [Z2 T-2 ~> m2 s-2].
    kappa_avg, & ! The time-weighted average of kappa [H Z T-1 ~> m2 s-1 or Pa s]
    tke_avg, &  ! The time-weighted average of TKE [Z2 T-2 ~> m2 s-2].
    N2_init, &  ! N2 as provided to this routine [T-2 ~> s-2].
    S2_init, &  ! S2 as provided to this routine [T-2 ~> s-2].
    N2_mean, &  ! The time-weighted average of N2 [T-2 ~> s-2].
    S2_mean, &  ! The time-weighted average of S2 [T-2 ~> s-2].
    I_dz_int, & ! The inverse of the distance between velocity & density points
                ! above and below an interface [Z-1 ~> m-1].
    u, v, &     ! The velocities after a timestep of mixing [L T-1 ~> m s-1].
    T, &        ! The potential temperature after a timestep of mixing [C ~> degC].
    Sal, &      ! The salinity after a timestep of mixing [S ~> ppt].
    pressure, & ! The pressure at an interface [R L2 T-2 ~> Pa].
    T_int, &    ! The temperature interpolated to an interface [C ~> degC].
    Sal_int, &  ! The salinity interpolated to an interface [S ~> ppt].
    dbuoy_dT, & ! The partial derivative of buoyancy with changes in temperature [Z T-2 C-1 ~> m s-2 degC-1]
    dbuoy_dS, & ! The partial derivative of buoyancy with changes in salinity [Z T-2 S-1 ~> m s-2 ppt-1]
    dSpV_dT, &  ! The partial derivative of specific volume with changes in temperature [R-1 C-1 ~> m3 kg-1 degC-1]
    dSpV_dS, &  ! The partial derivative of specific volume with changes in salinity [R-1 S-1 ~> m3 kg-1 ppt-1]
    rho_int, &  ! The in situ density interpolated to an interface [R ~> kg m-3]
    kappa_full, & ! kappa mapped back to the original interfaces [H Z T-1 ~> m2 s-1 or Pa s]
    tke_full, & ! tke mapped back to the original interfaces [Z2 T-2 ~> m2 s-2].
    kf          ! The fractional weight of interface kc+1 for
                ! interpolating back to the original index space [nondim].
  real, dimension(SZK_(GV)+1) :: &
    a1_col, &   ! a1 is the coupling between adjacent interfaces in the TKE,
                ! velocity, and density equations [H ~> m or kg m-2]
    c1_col      ! c1 is used in the tridiagonal (and similar) solvers [nondim].

  real :: f2    ! The squared Coriolis parameter of each column [T-2 ~> s-2].
  real :: surface_pres  ! The top surface pressure [R L2 T-2 ~> Pa].

  real :: dz_in_lay     !   The running sum of the thickness in a layer [H ~> m or kg m-2]
  real :: k0dt          ! The background diffusivity times the timestep [H Z ~> m2 or kg m-1]
  real :: dz_massless   ! A layer thickness that is considered massless [H ~> m or kg m-2]
  real :: I_hwt         ! The inverse of the sum of the adjacent masked thickness weights [H-1 ~> m-1 or m2 kg-1]
  real :: I_htot        ! The inverse of the sum of the thicknesses at adjacent vertices [H-1 ~> m-1 or m2 kg-1]
  real :: I_Prandtl     ! The inverse of the turbulent Prandtl number [nondim].
  real :: Prandtl_turb  ! A local copy of the turbulent Prandtl number [nondim].
  real :: b1            ! The inverse of the pivot in the tridiagonal equations [H-1 ~> m-1 or m2 kg-1].
  real :: d1            ! 1 - c1 in the tridiagonal equations [nondim]
  real :: bd1           ! A term in the denominator of b1 [H ~> m or kg m-2].
  real :: gR0           ! A conversion factor from H to pressure, Rho_0 times g in Boussinesq
                        ! mode, or just g when non-Boussinesq [R L2 T-2 H-1 ~> kg m-2 s-2 or m s-2].
  real :: g_R0          ! g_R0 is a rescaled version of g/Rho [Z R-1 T-2 ~> m4 kg-1 s-2].
  logical :: use_temperature  !  If true, temperature and salinity have been
                        ! allocated and are being used as state variables.

  integer, dimension(merge(G%iecB-(G%isc-1)+2, CS%niblock, CS%niblock==0), &
                     merge(G%jecB-(G%jsc-1)+2, CS%njblock, CS%njblock==0), SZK_(GV)+1) :: &
    kc                  ! The index map between the original
                        ! interfaces and the interfaces with massless layers
                        ! merged into nearby massive layers.
  integer, dimension(merge(G%iecB-(G%isc-1)+2, CS%niblock, CS%niblock==0), &
                     merge(G%jecB-(G%jsc-1)+2, CS%njblock, CS%njblock==0)) :: &
    nzc_2d              ! The number of layers in each column after massless layers
                        ! have been merged into nearby massive layers.
  integer, dimension(3,2) :: EOSdom  ! The 1-based domain of indices for the density derivatives
  real :: h_SW, h_SE, h_NW, h_NE ! Thicknesses at adjacent vertices [H ~> m or kg m-2]
  real :: mks_to_HZ_T   ! A factor used to restore dimensional scaling after the geometric mean
                        ! diffusivity is taken using thickness weighted powers [H Z s m-2 T-1 ~> 1]
                        ! or [H Z m s kg-1 T-1 ~> 1]
  real :: H_tiny        ! A sub-roundoff thickness to use in the denominator when calculating
                        ! thickness-weighted averages [H ~> m or kg m-2]
  integer :: IsB, IeB, JsB, JeB, i, j, k, nz, nzc
  integer :: nIIB, nJJB ! The block sizes actually used, with 0 resolved to the full domain.
  integer :: delta_i, delta_j  ! The block strides, one less than the block sizes because
                               ! consecutive blocks overlap by a column and a row.
  integer :: IsbB, IebB        ! I-block bounds written by the current block.
  integer :: JsbB, JebB        ! J-block bounds written by the current block.
  integer :: IebB_read         ! Last I index filled into the block, one beyond IebB where there
                               ! is room, to supply the II+1 that the corner interpolation reads.
  integer :: JebB_read         ! Last J index filled into the block, as IebB_read but in J.
  integer :: IIe, JJe          ! Valid block-local end bounds for face-indexed loops.
  integer :: ii, jj            ! Block-local 1-based i and j indices
  integer :: nzc_max           ! Max value of nzc across all columns

  ! Diagnostics that should be deleted?
  isB = G%isc-1 ; ieB = G%iecB ; jsB = G%jsc-1 ; jeB = G%jecB ; nz = GV%ke

  !   Consecutive blocks overlap by one column and row, so only nIIB-1 by nJJB-1
  ! vertices are written per block.
  nIIB = merge(IeB-IsB+2, CS%niblock, CS%niblock==0)
  nJJB = merge(JeB-JsB+2, CS%njblock, CS%njblock==0)
  delta_i = nIIB - 1 ; delta_j = nJJB - 1

  use_temperature = associated(tv%T)

  k0dt =  dt*CS%kappa_0
  gR0 = GV%H_to_RZ * GV%g_Earth
  g_R0 = GV%g_Earth_Z_T2 / GV%Rho0
  dz_massless = 0.1*sqrt((US%Z_to_m*GV%m_to_H)*k0dt)
  I_Prandtl = 0.0 ; if (CS%Prandtl_turb > 0.0) I_Prandtl = 1.0 / CS%Prandtl_turb
  Prandtl_turb = CS%Prandtl_turb
  H_tiny = 0.5 * GV%H_subroundoff
  nzc_max = 0

  !$omp target enter data map(alloc: diag_N2_init, diag_S2_init, diag_N2_mean, diag_S2_mean, &
  !$omp &                            dz_3d, h_at_u, h_at_v, kappa_vertex, h_vrt, dz_vrt, &
  !$omp &                            u_vrt, v_vrt, T_vrt, S_vrt, rho_vrt, kappa_full, tke_full, &
  !$omp &                            Idz, h_lay, dz_lay, u0xdz, v0xdz, T0xdz, S0xdz, &
  !$omp &                            kappa, tke, kappa_avg, tke_avg, N2_init, S2_init, N2_mean, S2_mean, &
  !$omp &                            I_dz_int, u, v, T, Sal, pressure, T_int, Sal_int, dbuoy_dT, dbuoy_dS, &
  !$omp &                            dSpV_dT, dSpV_dS, rho_int, a1_col, c1_col, kc, kf, nzc_2d)

  if ((CS%id_N2_init>0) .or. CS%debug) diag_N2_init(:,:,:) = 0.0
  if ((CS%id_S2_init>0) .or. CS%debug) diag_S2_init(:,:,:) = 0.0
  if (CS%id_N2_mean>0) diag_N2_mean(:,:,:) = 0.0
  if (CS%id_S2_mean>0) diag_S2_mean(:,:,:) = 0.0
  kappa_vertex(:,:,:) = 0.0

  ! The blocks stride by delta_i and delta_j rather than by the block size, because
  ! the corner interpolation reads the ii+1 and jj+1 neighbours of h_at_v and h_at_u;
  do JsbB=JsB,JeB,delta_j ; do IsbB=IsB,IeB,delta_i
    IebB = min(IsbB + delta_i - 1, IeB)
    JebB = min(JsbB + delta_j - 1, JeB)
    IebB_read = min(IsbB + nIIB - 1, IeB+1)
    JebB_read = min(JsbB + nJJB - 1, JeB+1)
    IIe = IebB - IsbB + 1 ; JJe = JebB - JsbB + 1

    ! Convert layer thicknesses into geometric thickness in height units.
    call thickness_to_dz(h, tv, dz_3d, G, GV, US, niB=nIIB, njB=nJJB, do_offload=.true., &
                         i_lo=IsbB, i_hi=IebB_read, j_lo=JsbB, j_hi=JebB_read)

    if (CS%vertex_shear_OBC_bug) then
      do concurrent( k=1:nz, JJ=1:JebB_read-JsbB+1, II=1:IIe ) DO_LOCALITY(local(I, j))
        I = IsbB + II - 1 ; j = JsbB + JJ - 1
        h_at_u(II,JJ,k) = G%mask2dCu(I,j) * (h(i,j,k) + h(i+1,j,k)) * 0.5
      enddo
      do concurrent( k=1:nz, JJ=1:JJe, II=1:IebB_read-IsbB+1 ) DO_LOCALITY(local(i, J))
        i = IsbB + II - 1 ; J = JsbB + JJ - 1
        h_at_v(II,JJ,k) = G%mask2dCv(i,J) * (h(i,j,k) + h(i,j+1,k)) * 0.5
      enddo
    else
      ! Because G%mask2dCu(I,j) is zero if either G%mask2dT(i,j) or G%mask2dT(i+1,j) except at OBC
      ! faces, the following form give equivalent answers to those above unless OBCs are in use,
      ! although the former is clearly less complicated and costly.
      do concurrent( k=1:nz, JJ=1:JebB_read-JsbB+1, II=1:IIe ) DO_LOCALITY(local(I, j))
        I = IsbB + II - 1 ; j = JsbB + JJ - 1
        h_at_u(II,JJ,k) = G%mask2dCu(I,j) * (G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i+1,j) * h(i+1,j,k)) / &
                                        (G%mask2dT(i,j) + G%mask2dT(i+1,j) + 1.0e-36)
      enddo
      do concurrent( k=1:nz, JJ=1:JJe, II=1:IebB_read-IsbB+1 ) DO_LOCALITY(local(i, J))
        i = IsbB + II - 1 ; J = JsbB + JJ - 1
        h_at_v(II,JJ,k) = G%mask2dCv(i,J) * (G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i,j+1) * h(i,j+1,k)) / &
                                        (G%mask2dT(i,j) + G%mask2dT(i,j+1) + 1.0e-36)
      enddo
    endif


    ! Interpolate the various quantities to the corners, using masks.
    do concurrent( k=1:nz, JJ=1:JJe, II=1:IIe ) DO_LOCALITY(local(I_hwt, I, J))
      I = IsbB + II - 1 ; J = JsbB + JJ - 1
      u_vrt(II,JJ,k) = ( (u_in(I,j,k) * h_at_u(II,JJ,k)) + (u_in(I,j+1,k) * h_at_u(II,JJ+1,k)) ) / &
                  ( (h_at_u(II,JJ,k) + h_at_u(II,JJ+1,k)) + H_tiny )
      v_vrt(II,JJ,k) = ( (v_in(i,J,k) * h_at_v(II,JJ,k)) + (v_in(i+1,J,k) * h_at_v(II+1,JJ,k)) ) / &
                  ( (h_at_v(II,JJ,k) + h_at_v(II+1,JJ,k)) + H_tiny )
      I_hwt = 1.0 / (((G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) + &
                      (G%mask2dT(i+1,j) * h(i+1,j,k) + G%mask2dT(i,j+1) * h(i,j+1,k))) + &
                     GV%H_subroundoff)
      if (use_temperature) then
        T_vrt(II,JJ,k) = ( (G%mask2dT(i,j) * (h(i,j,k) * T_in(i,j,k)) + &
                       G%mask2dT(i+1,j+1) * (h(i+1,j+1,k) * T_in(i+1,j+1,k))) + &
                      (G%mask2dT(i+1,j) * (h(i+1,j,k) * T_in(i+1,j,k)) + &
                       G%mask2dT(i,j+1) * (h(i,j+1,k) * T_in(i,j+1,k))) ) * I_hwt
        S_vrt(II,JJ,k) = ( (G%mask2dT(i,j) * (h(i,j,k) * S_in(i,j,k)) + &
                       G%mask2dT(i+1,j+1) * (h(i+1,j+1,k) * S_in(i+1,j+1,k))) + &
                      (G%mask2dT(i+1,j) * (h(i+1,j,k) * S_in(i+1,j,k)) + &
                       G%mask2dT(i,j+1) * (h(i,j+1,k) * S_in(i,j+1,k))) ) * I_hwt
      endif
      h_vrt(II,JJ,k) = ((G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) + &
                   (G%mask2dT(i+1,j) * h(i+1,j,k) + G%mask2dT(i,j+1) * h(i,j+1,k)) ) / &
                  ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + &
                   (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)) + 1.0e-36 )
      dz_vrt(II,JJ,k) = ((G%mask2dT(i,j) * dz_3d(II,JJ,k) + G%mask2dT(i+1,j+1) * dz_3d(II+1,JJ+1,k)) + &
                    (G%mask2dT(i+1,j) * dz_3d(II+1,JJ,k) + G%mask2dT(i,j+1) * dz_3d(II,JJ+1,k)) ) / &
                   ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + &
                    (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)) + 1.0e-36 )
!       h_vrt(II,JJ,k) = 0.25*((h(i,j,k) + h(i+1,j+1,k)) + (h(i+1,j,k) + h(i,j+1,k)))
!       h_vrt(II,JJ,k) = (((h(i,j,k)**2) + (h(i+1,j+1,k)**2)) + &
!                    ((h(i+1,j,k)**2) + (h(i,j+1,k)**2))) * I_hwt
    enddo

    if (.not.use_temperature) then
      do concurrent( k=1:nz, JJ=1:JJe, II=1:IIe )
        rho_vrt(II,JJ,k) = GV%Rlay(k)
      enddo
    endif

!---------------------------------------
! Work on each column.
!---------------------------------------
    ! Zero out blocked arrays so that the calculate_density_derivs call is defined everywhere it
    ! evaluates, and seed kappa.  Doing this outside of the column loop so that it vectorizes on CPU.
    do concurrent( K=1:nz+1, JJ=1:JJe, II=1:IIe )
      pressure(II,JJ,K) = 0.0 ; T_int(II,JJ,K) = 0.0 ; Sal_int(II,JJ,K) = 0.0
      kappa(II,JJ,K) = CS%kappa_seed
    enddo
    do concurrent( JJ=1:JJe, II=1:IIe )
      nzc_2d(II,JJ) = 0
    enddo

    !$omp target teams loop collapse(2) &
    !$omp & private(nzc, surface_pres, dz_in_lay, b1, d1, bd1, I, J, a1_col, c1_col) &
    !$omp & reduction(max:nzc_max)
    do JJ=1,JJe ; do II=1,IIe
      I = IsbB + II - 1 ; J = JsbB + JJ - 1

      if ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
          (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) > 0.0) then
        ! Store a transposed version of the initial arrays.
        ! Any elimination of massless layers would occur here.
        if (CS%eliminate_massless) then
          nzc = 1
          do k=1,nz
            ! Zero out the thicknesses of all layers, even if they are unused.
            h_lay(II,JJ,k) = 0.0 ; dz_lay(II,JJ,k) = 0.0 ; u0xdz(II,JJ,k) = 0.0 ; v0xdz(II,JJ,k) = 0.0
            T0xdz(II,JJ,k) = 0.0 ; S0xdz(II,JJ,k) = 0.0

            ! Add a new layer if this one has mass.
  !          if ((h_lay(II,JJ,nzc) > 0.0) .and. (h_vrt(II,JJ,k) > dz_massless)) nzc = nzc+1
            if ((k>CS%nkml) .and. (h_lay(II,JJ,nzc) > 0.0) .and. &
                (h_vrt(II,JJ,k) > dz_massless)) nzc = nzc+1

            ! Only merge clusters of massless layers.
  !         if ((h_lay(II,JJ,nzc) > dz_massless) .or. &
  !             ((h_lay(II,JJ,nzc) > 0.0) .and. (h_vrt(II,JJ,k) > dz_massless))) nzc = nzc+1

            kc(II,JJ,k) = nzc
            h_lay(II,JJ,nzc) = h_lay(II,JJ,nzc) + h_vrt(II,JJ,k)
            dz_lay(II,JJ,nzc) = dz_lay(II,JJ,nzc) + dz_vrt(II,JJ,k)
            u0xdz(II,JJ,nzc) = u0xdz(II,JJ,nzc) + u_vrt(II,JJ,k)*h_vrt(II,JJ,k)
            v0xdz(II,JJ,nzc) = v0xdz(II,JJ,nzc) + v_vrt(II,JJ,k)*h_vrt(II,JJ,k)
            if (use_temperature) then
              T0xdz(II,JJ,nzc) = T0xdz(II,JJ,nzc) + T_vrt(II,JJ,k)*h_vrt(II,JJ,k)
              S0xdz(II,JJ,nzc) = S0xdz(II,JJ,nzc) + S_vrt(II,JJ,k)*h_vrt(II,JJ,k)
            else
              T0xdz(II,JJ,nzc) = T0xdz(II,JJ,nzc) + rho_vrt(II,JJ,k)*h_vrt(II,JJ,k)
              S0xdz(II,JJ,nzc) = S0xdz(II,JJ,nzc) + rho_vrt(II,JJ,k)*h_vrt(II,JJ,k)
            endif
          enddo
          kc(II,JJ,nz+1) = nzc+1

          ! Set up Idz as the inverse of layer thicknesses.
          do k=1,nzc
            Idz(II,JJ,k) = 1.0 / h_lay(II,JJ,k)
          enddo

          !   Now determine kf, the fractional weight of interface kc when
          ! interpolating between interfaces kc and kc+1.
          kf(II,JJ,1) = 0.0 ; dz_in_lay = h_vrt(II,JJ,1)
          do k=2,nz
            if (kc(II,JJ,k) > kc(II,JJ,k-1)) then
              kf(II,JJ,k) = 0.0 ; dz_in_lay = h_vrt(II,JJ,k)
            else
              kf(II,JJ,k) = dz_in_lay*Idz(II,JJ,kc(II,JJ,k)) ; dz_in_lay = dz_in_lay + h_vrt(II,JJ,k)
            endif
          enddo
          kf(II,JJ,nz+1) = 0.0
        else
          do k=1,nz
            h_lay(II,JJ,k) = h_vrt(II,JJ,k)
            dz_lay(II,JJ,k) = dz_vrt(II,JJ,k)
            u0xdz(II,JJ,k) = u_vrt(II,JJ,k)*h_lay(II,JJ,k) ; v0xdz(II,JJ,k) = v_vrt(II,JJ,k)*h_lay(II,JJ,k)
          enddo
          if (use_temperature) then
            do k=1,nz
              T0xdz(II,JJ,k) = T_vrt(II,JJ,k)*h_lay(II,JJ,k) ; S0xdz(II,JJ,k) = S_vrt(II,JJ,k)*h_lay(II,JJ,k)
            enddo
          else
            do k=1,nz
              T0xdz(II,JJ,k) = rho_vrt(II,JJ,k)*h_lay(II,JJ,k) ; S0xdz(II,JJ,k) = rho_vrt(II,JJ,k)*h_lay(II,JJ,k)
            enddo
          endif
          nzc = nz
          do k=1,nzc+1
            kc(II,JJ,k) = k
            kf(II,JJ,k) = 0.0
          enddo
        endif

        !   Set up I_dz_int as the inverse of the distance between adjacent layer centers, for
        ! applying the background diffusivity.
        I_dz_int(II,JJ,1) = 2.0 / dz_lay(II,JJ,1)
        do K=2,nzc
          I_dz_int(II,JJ,K) = 2.0 / (dz_lay(II,JJ,K-1) + dz_lay(II,JJ,K))
        enddo
        I_dz_int(II,JJ,nzc+1) = 2.0 / dz_lay(II,JJ,nzc)

        !   Determine the velocities and thermodynamic tracers after eliminating massless layers
        ! and applying a time-step of background diffusion.
        if (nzc > 1) then
          a1_col(2) = k0dt*I_dz_int(II,JJ,2)
          b1 = 1.0 / (h_lay(II,JJ,1) + a1_col(2))
          u(II,JJ,1) = b1 * u0xdz(II,JJ,1) ; v(II,JJ,1) = b1 * v0xdz(II,JJ,1)
          T(II,JJ,1) = b1 * T0xdz(II,JJ,1) ; Sal(II,JJ,1) = b1 * S0xdz(II,JJ,1)
          c1_col(2) = a1_col(2) * b1 ; d1 = h_lay(II,JJ,1) * b1 ! = 1 - c1
          do k=2,nzc-1
            bd1 = h_lay(II,JJ,k) + d1*a1_col(k)
            a1_col(k+1) = k0dt*I_dz_int(II,JJ,k+1)
            b1 = 1.0 / (bd1 + a1_col(k+1))
            u(II,JJ,k) = b1 * (u0xdz(II,JJ,k) + a1_col(k)*u(II,JJ,k-1))
            v(II,JJ,k) = b1 * (v0xdz(II,JJ,k) + a1_col(k)*v(II,JJ,k-1))
            T(II,JJ,k) = b1 * (T0xdz(II,JJ,k) + a1_col(k)*T(II,JJ,k-1))
            Sal(II,JJ,k) = b1 * (S0xdz(II,JJ,k) + a1_col(k)*Sal(II,JJ,k-1))
            c1_col(k+1) = a1_col(k+1) * b1 ; d1 = bd1 * b1 ! d1 = 1 - c1
          enddo
          ! rho or T and S have insulating boundary conditions, u & v use no-slip
          ! bottom boundary conditions (if kappa0 > 0).
          ! For no-slip bottom boundary conditions
          b1 = 1.0 / ((h_lay(II,JJ,nzc) + d1*a1_col(nzc)) + k0dt*I_dz_int(II,JJ,nzc+1))
          u(II,JJ,nzc) = b1 * (u0xdz(II,JJ,nzc) + a1_col(nzc)*u(II,JJ,nzc-1))
          v(II,JJ,nzc) = b1 * (v0xdz(II,JJ,nzc) + a1_col(nzc)*v(II,JJ,nzc-1))
          ! For insulating boundary conditions
          b1 = 1.0 / (h_lay(II,JJ,nzc) + d1*a1_col(nzc))
          T(II,JJ,nzc) = b1 * (T0xdz(II,JJ,nzc) + a1_col(nzc)*T(II,JJ,nzc-1))
          Sal(II,JJ,nzc) = b1 * (S0xdz(II,JJ,nzc) + a1_col(nzc)*Sal(II,JJ,nzc-1))
          do k=nzc-1,1,-1
            u(II,JJ,k) = u(II,JJ,k) + c1_col(k+1)*u(II,JJ,k+1)
            v(II,JJ,k) = v(II,JJ,k) + c1_col(k+1)*v(II,JJ,k+1)
            T(II,JJ,k) = T(II,JJ,k) + c1_col(k+1)*T(II,JJ,k+1)
            Sal(II,JJ,k) = Sal(II,JJ,k) + c1_col(k+1)*Sal(II,JJ,k+1)
          enddo
        else
          ! This is correct, but probably unnecessary.
          b1 = 1.0 / (h_lay(II,JJ,1) + k0dt*I_dz_int(II,JJ,2))
          u(II,JJ,1) = b1 * u0xdz(II,JJ,1) ; v(II,JJ,1) = b1 * v0xdz(II,JJ,1)
          b1 = 1.0 / h_lay(II,JJ,1)
          T(II,JJ,1) = b1 * T0xdz(II,JJ,1) ; Sal(II,JJ,1) = b1 * S0xdz(II,JJ,1)
        endif

        ! Get T, S and pressure at the interfaces for calculating dbuoy_dT and dbuoy_dS.
        surface_pres = 0.0
        if (associated(p_surf)) then
          if (CS%psurf_bug) then
            ! This is wrong because it is averaging values from land in some places.
            surface_pres = 0.25 * ((p_surf(i,j) + p_surf(i+1,j+1)) + &
                                   (p_surf(i+1,j) + p_surf(i,j+1)))
          else
            surface_pres = ((G%mask2dT(i,j) * p_surf(i,j) + G%mask2dT(i+1,j+1) * p_surf(i+1,j+1)) + &
                            (G%mask2dT(i+1,j) * p_surf(i+1,j) + G%mask2dT(i,j+1) * p_surf(i,j+1)) ) / &
                           ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + &
                            (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)) + 1.0e-36 )
          endif
        endif
        if (use_temperature) then
          pressure(II,JJ,1) = surface_pres
          do k=2,nzc
            pressure(II,JJ,k) = pressure(II,JJ,k-1) + gR0*h_lay(II,JJ,k-1)
            T_int(II,JJ,k) = 0.5*(T(II,JJ,k-1) + T(II,JJ,k))
            Sal_int(II,JJ,k) = 0.5*(Sal(II,JJ,k-1) + Sal(II,JJ,k))
          enddo
        endif

        ! Save the number of merged layers for the stages below.
        nzc_2d(II,JJ) = nzc
        nzc_max = max( nzc_max, nzc )
      endif ! end if ((G%mask2dCu(I,j)
    enddo ; enddo ! end of the setup I- and J-loops

    ! Calculate the thermodynamic coefficients for all of the vertex columns at once.
    if (use_temperature) then
      if (GV%Boussinesq .or. GV%semi_Boussinesq) then
        EOSdom(1,1) = 1 ; EOSdom(1,2) = IIe
        EOSdom(2,1) = 1 ; EOSdom(2,2) = JJe
        EOSdom(3,1) = 2 ; EOSdom(3,2) = nzc_max
        call calculate_density_derivs(T_int, Sal_int, pressure, dbuoy_dT, dbuoy_dS, &
                                      tv%eqn_of_state, EOSdom, scale=-g_R0 )
      else
        ! These should perhaps be combined into a single call to calculate the thermal expansion
        ! and haline contraction coefficients?
        !$omp target update from(T_int, Sal_int, pressure, nzc_2d)
        do JJ=1,JJe ; do II=1,IIe
          if (nzc_2d(II,JJ) > 0) then
            nzc = nzc_2d(II,JJ)
            call calculate_specific_vol_derivs(T_int(II,JJ,:), Sal_int(II,JJ,:), pressure(II,JJ,:), &
                                          dSpV_dT(II,JJ,:), dSpV_dS(II,JJ,:), tv%eqn_of_state, (/2,nzc/) )
            call calculate_density(T_int(II,JJ,:), Sal_int(II,JJ,:), pressure(II,JJ,:), &
                                   rho_int(II,JJ,:), tv%eqn_of_state, (/2,nzc/) )
            do K=2,nzc
              dbuoy_dT(II,JJ,K) = GV%g_Earth_Z_T2 * (rho_int(II,JJ,K) * dSpV_dT(II,JJ,K))
              dbuoy_dS(II,JJ,K) = GV%g_Earth_Z_T2 * (rho_int(II,JJ,K) * dSpV_dS(II,JJ,K))
            enddo
          endif
        enddo ; enddo
        !$omp target update to(dbuoy_dT, dbuoy_dS)
      endif
    elseif (GV%Boussinesq .or. GV%semi_Boussinesq) then
      !$omp target update from(nzc_2d)
      do JJ=1,JJe ; do II=1,IIe
        if (nzc_2d(II,JJ) > 0) then
          nzc = nzc_2d(II,JJ)
          do K=1,nzc+1 ; dbuoy_dT(II,JJ,K) = -g_R0 ; dbuoy_dS(II,JJ,K) = 0.0 ; enddo
        endif
      enddo ; enddo
      !$omp target update to(dbuoy_dT, dbuoy_dS)
    else
      !$omp target update from(nzc_2d)
      do JJ=1,JJe ; do II=1,IIe
        if (nzc_2d(II,JJ) > 0) then
          nzc = nzc_2d(II,JJ)
          do K=1,nzc+1 ; dbuoy_dS(II,JJ,K) = 0.0 ; enddo
          dbuoy_dT(II,JJ,1) = -GV%g_Earth_Z_T2 / GV%Rlay(1)
          do K=2,nzc
            dbuoy_dT(II,JJ,K) = -GV%g_Earth_Z_T2 / (0.5*(GV%Rlay(k-1) + GV%Rlay(k)))
          enddo
          dbuoy_dT(II,JJ,nzc+1) = -GV%g_Earth_Z_T2 / GV%Rlay(nzc)
        endif
      enddo ; enddo
      !$omp target update to(dbuoy_dT, dbuoy_dS)
    endif

  !---------------------------------------
  ! Work on each column.
  !---------------------------------------
    do concurrent( JJ=1:JJe, II=1:IIe ) DO_LOCALITY(local(nzc, f2, I, J))
      I = IsbB + II - 1 ; J = JsbB + JJ - 1
      if ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
          (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) > 0.0) then
        nzc = nzc_2d(II,JJ)
        f2 = G%Coriolis2Bu(I,J)

        call kappa_shear_column(kappa, tke, dt, nzc, f2, dbuoy_dT, dbuoy_dS, &
                                h_lay, dz_lay, I_dz_int, kappa_avg, u, v, T, Sal, &
                                tke_avg, N2_init, S2_init, N2_mean, S2_mean, &
                                CS, GV, US, nIIB, nJJB, II, JJ)

      ! Extrapolate from the vertically reduced grid back to the original layers.
        if (nz == nzc) then
          do K=1,nz+1
            kappa_full(II,JJ,K) = kappa_avg(II,JJ,K)
            if (CS%all_layer_TKE_bug) then
              tke_full(II,JJ,K) = tke(II,JJ,K)
            else
              tke_full(II,JJ,K) = tke_avg(II,JJ,K)
            endif
          enddo
          if (CS%id_N2_mean>0) then ; do K=1,nz+1
            diag_N2_mean(I,J,K) = N2_mean(II,JJ,K)
          enddo ; endif
          if (CS%id_S2_mean>0) then ; do K=1,nz+1
            diag_S2_mean(I,J,K) = S2_mean(II,JJ,K)
          enddo ; endif
          if ((CS%id_N2_init>0) .or. CS%debug) then ; do K=1,nz+1
            diag_N2_init(I,J,K) = N2_init(II,JJ,K)
          enddo ; endif
          if ((CS%id_S2_init>0) .or. CS%debug) then ; do K=1,nz+1
            diag_S2_init(I,J,K) = S2_init(II,JJ,K)
          enddo ; endif
        else
          do K=1,nz+1
            if (kf(II,JJ,K) == 0.0) then
              kappa_full(II,JJ,K) = kappa_avg(II,JJ,kc(II,JJ,K))
              tke_full(II,JJ,K) = tke_avg(II,JJ,kc(II,JJ,K))
            else
              kappa_full(II,JJ,K) = (1.0-kf(II,JJ,K)) * kappa_avg(II,JJ,kc(II,JJ,K)) + &
                                  kf(II,JJ,K) * kappa_avg(II,JJ,kc(II,JJ,K)+1)
              tke_full(II,JJ,K) = (1.0-kf(II,JJ,K)) * tke_avg(II,JJ,kc(II,JJ,K)) + &
                                kf(II,JJ,K) * tke_avg(II,JJ,kc(II,JJ,K)+1)
            endif
          enddo
          do K=1,nz+1
            if (kf(II,JJ,K) == 0.0) then
              if (CS%id_N2_mean>0) diag_N2_mean(I,J,K) = N2_mean(II,JJ,kc(II,JJ,K))
              if (CS%id_S2_mean>0) diag_S2_mean(I,J,K) = S2_mean(II,JJ,kc(II,JJ,K))
              if ((CS%id_N2_init>0) .or. CS%debug) diag_N2_init(I,J,K) = N2_init(II,JJ,kc(II,JJ,K))
              if ((CS%id_S2_init>0) .or. CS%debug) diag_S2_init(I,J,K) = S2_init(II,JJ,kc(II,JJ,K))
            else
              if (CS%id_N2_mean>0) &
                diag_N2_mean(I,J,K) = (1.0-kf(II,JJ,K)) * N2_mean(II,JJ,kc(II,JJ,K)) + &
                                      kf(II,JJ,K) * N2_mean(II,JJ,kc(II,JJ,K)+1)
              if (CS%id_S2_mean>0) &
                diag_S2_mean(I,J,K) = (1.0-kf(II,JJ,K)) * S2_mean(II,JJ,kc(II,JJ,K)) + &
                                      kf(II,JJ,K) * S2_mean(II,JJ,kc(II,JJ,K)+1)
              if ((CS%id_N2_init>0) .or. CS%debug) &
                diag_N2_init(I,J,K) = (1.0-kf(II,JJ,K)) * N2_init(II,JJ,kc(II,JJ,K)) + &
                                      kf(II,JJ,K) * N2_init(II,JJ,kc(II,JJ,K)+1)
              if ((CS%id_S2_init>0) .or. CS%debug) &
                diag_S2_init(I,J,K) = (1.0-kf(II,JJ,K)) * S2_init(II,JJ,kc(II,JJ,K)) + &
                                      kf(II,JJ,K) * S2_init(II,JJ,kc(II,JJ,K)+1)
            endif
          enddo
        endif
      else  ! Land points.
        do K=1,nz+1
          kappa_full(II,JJ,K) = 0.0 ; tke_full(II,JJ,K) = 0.0
        enddo
      endif
    enddo ! end of the column I- and J-loops

    ! Store the results for restarts or interpolation back to tracer points.
    if (CS%VS_viscosity_bug) then
      do concurrent( K=1:nz+1, JJ=1:JJe, II=1:IIe ) DO_LOCALITY(local(I, J))
        I = IsbB + II - 1 ; J = JsbB + JJ - 1
        kappa_vertex(I,J,K) = kappa_full(II,JJ,K)
        tke_io(I,J,K) = G%mask2dBu(I,J) * tke_full(II,JJ,K)
        kv_io(I,J,K) = ( G%mask2dBu(I,J) * kappa_vertex(I,J,K) ) * Prandtl_turb
      enddo
    else
      do concurrent( K=1:nz+1, JJ=1:JJe, II=1:IIe ) DO_LOCALITY(local(I, J))
        I = IsbB + II - 1 ; J = JsbB + JJ - 1
        kappa_vertex(I,J,K) = kappa_full(II,JJ,K)
        tke_io(I,J,K) = tke_full(II,JJ,K)
        kv_io(I,J,K) = kappa_vertex(I,J,K) * Prandtl_turb
      enddo
    endif

  enddo ; enddo ! end of the i- and j-block loops

  ! Set the diffusivities in tracer columns from the values at vertices.

  do concurrent( j=G%jsc:G%jec, i=G%isc:G%iec )
    ! The turbulent length scales (and hence turbulent diffusivity) should always go to 0 at the top and bottom.
    kappa_io(i,j,1) = 0.0
    kappa_io(i,j,nz+1) = 0.0
  enddo
  if (CS%VS_ThicknessMean) then
    ! This conversion factor is required to allow for arbitrary fractional powers of the diffusivities.
    mks_to_HZ_T = 1.0 /  GV%HZ_T_to_MKS

    ! Stride by delta_i/j to account for stencils when calculating h at adjacent vertices.
    do jsbB=G%jsc,G%jec,delta_j ; do isbB=G%isc,G%iec,delta_i
      iebB = min(isbB + delta_i - 1, G%iec)
      jebB = min(jsbB + delta_j - 1, G%jec)

      ! +2 so that we write h_vrt over the entire nIIB x nJJB extent
      do concurrent( k=1:nz, jj=1:jebB-jsbB+2, ii=1:iebB-isbB+2 ) DO_LOCALITY(local(I, J))
        I = isbB + ii - 2 ; J = jsbB + jj - 2
        h_vrt(ii,jj,k) = ((G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) + &
                     (G%mask2dT(i+1,j) * h(i+1,j,k) + G%mask2dT(i,j+1) * h(i,j+1,k)) ) / &
                    ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + &
                     (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)) + 1.0e-36 )
      enddo

      if (CS%VS_GeometricMean) then
        ! +1 so that we only read h_vrt over delta_i/j
        do concurrent( K=2:nz, jj=1:jebB-jsbB+1, ii=1:iebB-isbB+1 ) &
                     DO_LOCALITY(local(h_SW, h_NE, h_NW, h_SE, I_htot, I, J))
          I = isbB + ii - 1 ; J = jsbB + jj - 1
          h_SW = 0.5 * (h_vrt(ii,jj,k) + h_vrt(ii,jj,k-1))
          h_NE = 0.5 * (h_vrt(ii+1,jj+1,k) + h_vrt(ii+1,jj+1,k-1))
          h_NW = 0.5 * (h_vrt(ii,jj+1,k) + h_vrt(ii,jj+1,k-1))
          h_SE = 0.5 * (h_vrt(ii+1,jj,k) + h_vrt(ii+1,jj,k-1))
          if ((h_SW + h_NE) + (h_NW + h_SE) > 0.0) then
            !  The geometric mean is zero if any component is zero, hence the need to use a floor
            !  on the value of kappa_trunc in regions on boundaries of shear zones.
            I_htot = 1.0 / ((h_SW + h_NE) + (h_NW + h_SE))
            kappa_io(i,j,K) = G%mask2dT(i,j) * mks_to_HZ_T * &
                            ( ((GV%HZ_T_to_MKS * max(kappa_vertex(I-1,J-1,K), CS%VS_GeoMean_Kdmin))**(h_SW*I_htot) * &
                               (GV%HZ_T_to_MKS * max(kappa_vertex(I,J,K), CS%VS_GeoMean_Kdmin))**(h_NE*I_htot)) * &
                              ((GV%HZ_T_to_MKS * max(kappa_vertex(I-1,J,K), CS%VS_GeoMean_Kdmin))**(h_NW*I_htot) * &
                               (GV%HZ_T_to_MKS * max(kappa_vertex(I,J-1,K), CS%VS_GeoMean_Kdmin))**(h_SE*I_htot)) )
          else
            ! If all points have zero thickness, the thickness-weighted geometric mean is undefined, so use
            ! the non-thickness weighted geometric mean instead.
            kappa_io(i,j,K) = G%mask2dT(i,j) * sqrt(sqrt( &
              (max(kappa_vertex(I-1,J-1,K),CS%VS_GeoMean_Kdmin) * max(kappa_vertex(I,J,K),CS%VS_GeoMean_Kdmin)) * &
              (max(kappa_vertex(I-1,J,K),CS%VS_GeoMean_Kdmin) * max(kappa_vertex(I,J-1,K),CS%VS_GeoMean_Kdmin)) ))
          endif
        enddo
      else   ! Use thickness-weighted arithmetic mean diffusivities.
        ! +1 so that we only read h_vrt over delta_i/j
        do concurrent( K=2:nz, jj=1:jebB-jsbB+1, ii=1:iebB-isbB+1 ) &
                     DO_LOCALITY(local(h_SW, h_NE, h_NW, h_SE, I_htot, I, J))
          I = isbB + ii - 1 ; J = jsbB + jj - 1
          h_SW = 0.5 * (h_vrt(ii,jj,k) + h_vrt(ii,jj,k-1))
          h_NE = 0.5 * (h_vrt(ii+1,jj+1,k) + h_vrt(ii+1,jj+1,k-1))
          h_NW = 0.5 * (h_vrt(ii,jj+1,k) + h_vrt(ii,jj+1,k-1))
          h_SE = 0.5 * (h_vrt(ii+1,jj,k) + h_vrt(ii+1,jj,k-1))
          ! The following expression is a thickness weighted arithmetic mean at tracer points:
          I_htot = 1.0 / (((h_SW + h_NE) + (h_NW + h_SE)) + GV%H_subroundoff)
          kappa_io(i,j,K) = G%mask2dT(i,j) * &
            (((kappa_vertex(I-1,J-1,K)*h_SW) + (kappa_vertex(I,J,K)*h_NE)) + &
             ((kappa_vertex(I-1,J,K)*h_NW) + (kappa_vertex(I,J-1,K)*h_SE))) * I_htot
        enddo
      endif
    enddo ; enddo ! end of the i- and j-block loops
  elseif (CS%VS_GeometricMean) then   ! The geometic mean diffusivities are not thickness weighted.
    do concurrent( K=2:nz, j=G%jsc:G%jec, i=G%isc:G%iec )
      kappa_io(i,j,K) = G%mask2dT(i,j) * sqrt(sqrt( &
            (max(kappa_vertex(I-1,J-1,K),CS%VS_GeoMean_Kdmin) * max(kappa_vertex(I,J,K),CS%VS_GeoMean_Kdmin)) * &
            (max(kappa_vertex(I-1,J,K),CS%VS_GeoMean_Kdmin) * max(kappa_vertex(I,J-1,K),CS%VS_GeoMean_Kdmin)) ))
    enddo
  else   ! Use a non-thickness weighted arithmetic mean.
    do concurrent( K=2:nz, j=G%jsc:G%jec, i=G%isc:G%iec )
      kappa_io(i,j,K) = G%mask2dT(i,j) * 0.25 * &
                       ((kappa_vertex(I-1,J-1,K) + kappa_vertex(I,J,K)) +&
                        (kappa_vertex(I-1,J,K) + kappa_vertex(I,J-1,K)))
    enddo
  endif

  !   The checksums and the diagnostics below all read host memory, so anything they touch has to
  ! be brought back first.  kappa_io, tke_io and kv_io are the caller's, and are still device
  ! resident here; the caller copies them back after this routine returns.
  if (CS%debug) then
    !$omp target update from( diag_N2_init, diag_S2_init, kappa_io, tke_io )
    call Bchksum(diag_N2_init, "shear_vertex N2_init", G%HI, unscale=US%s_to_T**2)
    call Bchksum(diag_S2_init, "shear_vertex S2_init", G%HI, unscale=US%s_to_T**2)
    call hchksum(kappa_io, "kappa", G%HI, unscale=GV%HZ_T_to_m2_s)
    call Bchksum(tke_io, "tke", G%HI, unscale=US%Z_to_m**2*US%s_to_T**2)
  endif

  if (CS%id_Kd_shear > 0) then
    !$omp target update from( kappa_io )
    call post_data(CS%id_Kd_shear, kappa_io, CS%diag)
  endif
  if (CS%id_TKE > 0) then
    !$omp target update from( tke_io )
    call post_data(CS%id_TKE, tke_io, CS%diag)
  endif
  if (CS%id_Kd_vertex > 0) then
    !$omp target update from( kappa_vertex )
    call post_data(CS%id_Kd_vertex, kappa_vertex, CS%diag)
  endif
  if (CS%id_N2_init > 0) then
    !$omp target update from( diag_N2_init )
    call post_data(CS%id_N2_init, diag_N2_init, CS%diag)
  endif
  if (CS%id_S2_init > 0) then
    !$omp target update from( diag_S2_init )
    call post_data(CS%id_S2_init, diag_S2_init, CS%diag)
  endif
  if (CS%id_N2_mean > 0) then
    !$omp target update from( diag_N2_mean )
    call post_data(CS%id_N2_mean, diag_N2_mean, CS%diag)
  endif
  if (CS%id_S2_mean > 0) then
    !$omp target update from( diag_S2_mean )
    call post_data(CS%id_S2_mean, diag_S2_mean, CS%diag)
  endif

  !$omp target exit data map(delete: diag_N2_init, diag_S2_init, diag_N2_mean, diag_S2_mean, &
  !$omp &                            dz_3d, h_at_u, h_at_v, kappa_vertex, h_vrt, dz_vrt, &
  !$omp &                            u_vrt, v_vrt, T_vrt, S_vrt, rho_vrt, kappa_full, tke_full, &
  !$omp &                            Idz, h_lay, dz_lay, u0xdz, v0xdz, T0xdz, S0xdz, &
  !$omp &                            kappa, tke, kappa_avg, tke_avg, N2_init, S2_init, N2_mean, S2_mean, &
  !$omp &                            I_dz_int, u, v, T, Sal, pressure, T_int, Sal_int, dbuoy_dT, dbuoy_dS, &
  !$omp &                            dSpV_dT, dSpV_dS, rho_int, a1_col, c1_col, kc, kf, nzc_2d)

end subroutine Calc_kappa_shear_vertex


!> This subroutine calculates shear-driven diffusivity and TKE in a single column
pure subroutine kappa_shear_column(kappa, tke, dt, nzc, f2, dbuoy_dT, dbuoy_dS, hlay, dz_lay, I_dz_int, &
                              kappa_avg, u, v, T, Sal,   tke_avg, N2_init, S2_init, &
                              N2_mean, S2_mean, CS, GV, US, nib, njb, ii, jj )
  !$omp declare target
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure.
  integer,           intent(in)    :: nib  !< The block size of the input arrays in the i-direction
  integer,           intent(in)    :: njb  !< The block size of the input arrays in the j-direction.
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(inout) :: kappa    !< The time-weighted average of kappa [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(out)   :: tke      !< The Turbulent Kinetic Energy per unit mass at
                                               !! an interface [Z2 T-2 ~> m2 s-2].
  integer,           intent(in)    :: nzc      !< The number of active layers in the column.
  real,              intent(in)    :: f2       !< The square of the Coriolis parameter [T-2 ~> s-2].
  real, dimension(1:nib,1:njb,SZK_(GV)), &
                     intent(in)    :: hlay     !< The layer thickness [H ~> m or kg m-2]
  real, dimension(1:nib,1:njb,SZK_(GV)), &
                     intent(in)    :: dz_lay   !< The geometric layer thickness in height units [Z ~> m]
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(in)    :: I_dz_int !< The inverse of the distance between velocity & density
                                               !! points above and below an interface [Z-1 ~> m-1].
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(in)    :: dbuoy_dT !< The derivative of buoyancy with respect to temperature
                                               !! [H L T-2 ~> m s-2 degC-1 or m s-2 kg-1 m2]
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(in)    :: dbuoy_dS !< The derivative of buoyancy with respect to salinity
                                               !! [H L T-2 ~> m s-2 ppt-1 or m s-2 kg-1 m2]
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(out)   :: kappa_avg !< The time-weighted average of kappa [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(inout) :: u !< On entry, the zonal velocity after eliminating massless layers
                                        !! and applying background diffusion; on exit, the projected
                                        !! zonal velocity [H Z T-1 ~> m s-1]
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(inout) :: v !< On entry, the meridional velocity after eliminating massless
                                        !! layers and applying background diffusion; on exit, the
                                        !! projected meridional velocity [H Z T-1 ~> m s-1]
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(inout) :: T !< On entry, the temperature after eliminating massless layers
                                        !! and applying background diffusion; on exit, the projected
                                        !! temperature [H Z T-1 ~> degC]
  real, dimension(1:nib,1:njb ,SZK_(GV)+1), &
                     intent(inout) :: Sal !< On entry, the salinity after eliminating massless layers
                                          !! and applying background diffusion; on exit, the projected
                                          !! salinity [H Z T-1 ~> ppt]
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(out)   :: tke_avg      !< The time-weighted average of TKE [Z2 T-2 ~> m2 s-2].
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(out)   :: N2_mean      !< The time-weighted average of N2 [Z2 T-2 ~> m2 s-2].
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(out)   :: S2_mean      !< The time-weighted average of S2 [Z2 T-2 ~> m2 s-2].
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(out)   :: N2_init      !< The initial value of N2 [Z2 T-2 ~> m2 s-2].
  real, dimension(1:nib,1:njb,SZK_(GV)+1), &
                     intent(out)   :: S2_init      !< The initial value of S2 [Z2 T-2 ~> m2 s-2].
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  type(Kappa_shear_CS),    intent(in)    :: CS     !< The control structure returned by a previous
                                                   !! call to kappa_shear_init.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  integer,                 intent(in)    :: ii      !< The block-local i index of the column.
  integer,                 intent(in)    :: jj      !< The block-local j index of the column.

  ! Local variables
  ! In GPU builds these locals have compile-time-constant sizes so each device call uses
  ! stack ("local memory") arrays; runtime-sized automatics are device-heap allocated per
  ! call, which exhausts the default heap and serializes on the device allocator.
#ifdef __NVCOMPILER_OPENMP_GPU
  real, dimension(GPU_nk_max) :: &
#else
  real, dimension(nzc) :: &
#endif
    Idz, &      ! The inverse of the distance between TKE points [Z-1 ~> m-1].
    u_test, v_test, & ! Temporary velocities [L T-1 ~> m s-1].
    T_test, S_test, & ! Temporary temperatures [C ~> degC] and salinities [S ~> ppt].
    u_col, v_col, & ! Contiguous per-column copies of u(ii,jj,:) and v(ii,jj,:) used as
                    ! input/output for calculate_projected_state, to avoid passing
                    ! non-contiguous array sections or mismatched-rank whole arrays.
    T_col, Sal_col, & ! Contiguous per-column copies of T(ii,jj,:) and Sal(ii,jj,:), as above.
    hlay_col, & ! A contiguous per-column copy of hlay(ii,jj,:), as above.
    dz_lay_col  ! A contiguous per-column copy of dz_lay(ii,jj,:).  The setup below walks
                ! the layer thicknesses several times over, and reading them from the
                ! caller's blocked arrays gives every one of those passes a stride of
                ! nib*njb elements.

#ifdef __NVCOMPILER_OPENMP_GPU
  real, dimension(GPU_nk_max+1) :: &
#else
  real, dimension(nzc+1) :: &
#endif
    N2, &       ! The squared buoyancy frequency at an interface [T-2 ~> s-2].
    h_Int, &    ! The extent of a finite-volume space surrounding an interface,
                ! as used in calculating kappa and TKE [H ~> m or kg m-2]
    dz_Int, &   ! The vertical distance with the space surrounding an interface,
                ! as used in calculating kappa and TKE [Z ~> m]
    dz_h_Int, & ! The ratio of the vertical distances to the thickness around an
                ! interface [Z H-1 ~> nondim or m3 kg-1].  In non-Boussinesq mode
                ! this is the specific volume, otherwise it is a scaling factor.
    I_dz_int_col, & ! A contiguous per-column copy of I_dz_int(ii,jj,:), the inverse of
                ! the distance between velocity & density points above and below an
                ! interface [Z-1 ~> m-1], used as above for calculate_projected_state.
    dbuoy_dT_col, & ! A contiguous per-column copy of dbuoy_dT(ii,jj,:), as above.
    dbuoy_dS_col, & ! A contiguous per-column copy of dbuoy_dS(ii,jj,:), as above.
    S2, &       ! The squared shear at an interface [T-2 ~> s-2].
    k_src, &    ! The shear-dependent source term in the kappa equation [T-1 ~> s-1]
    kappa_src, & ! The shear-dependent source term in the kappa equation [T-1 ~> s-1]
    kappa_out, & ! The kappa that results from the kappa equation [H Z T-1 ~> m2 s-1 or Pa s]
    kappa_mid, & ! The average of the initial and predictor estimates of kappa [H Z T-1 ~> m2 s-1 or Pa s]
    tke_pred, & ! The value of TKE from a predictor step [Z2 T-2 ~> m2 s-2].
    kappa_pred, & ! The value of kappa from a predictor step [H Z T-1 ~> m2 s-1 or Pa s]
    kappa_col, & ! A contiguous per-column copy of kappa(ii,jj,:) used as an input to
                 ! find_kappa_tke, to avoid passing a non-contiguous array section.
    tke_col, &   ! A contiguous per-column copy of tke(ii,jj,:) used as an output from
                 ! find_kappa_tke, to avoid passing a non-contiguous array section.
    kappa_avg_col, & ! Contiguous per-column accumulators for the time-weighted averages.
    tke_avg_col, &   ! The iteration below updates these once per interface per iteration;
    N2_mean_col, &   ! keeping them here rather than in the caller's blocked arrays keeps
    S2_mean_col, &   ! the whole iteration working on a few contiguous kilobytes that stay
                     ! resident in L1, instead of striding through the blocked arrays.
                     ! They are written back to the blocked arrays once, after the loop.
    pressure, & ! The pressure at an interface [R L2 T-2 ~> Pa].
    I_L2_bdry, &   ! The inverse of the square of twice the harmonic mean
                   ! distance to the top and bottom boundaries [H-1 Z-1 ~> m-2 or m kg-1].
    K_Q, &         ! Diffusivity divided by TKE [H T Z-1 ~> s or kg s m-3]
    K_Q_tmp, &     ! A temporary copy of diffusivity divided by TKE [H T Z-1 ~> s or kg s m-3]
    local_src_avg, & ! The time-integral of the local source [nondim]
    tol_min, & ! Minimum tolerated ksrc for the corrector step [T-1 ~> s-1].
    tol_max, & ! Maximum tolerated ksrc for the corrector step [T-1 ~> s-1].
    tol_chg, & ! The tolerated kappa change integrated over a timestep [nondim].
    dist_from_top, &  ! The distance from the top surface [Z ~> m].
    h_from_top, & ! The total thickness above an interface [H ~> m or kg m-2]
    local_src     ! The sum of all sources of kappa, including kappa_src and
                  ! sources from the elliptic term [T-1 ~> s-1]

  real :: dist_from_bot ! The distance from the bottom surface [Z ~> m].
  real :: h_from_bot    ! The total thickness below and interface [H ~> m or kg m-2]
  real :: gR0           ! A conversion factor from H to pressure, Rho_0 times g in Boussinesq
                        ! mode, or just g when non-Boussinesq [R L2 T-2 H-1 ~> kg m-2 s-2 or m s-2].
  real :: g_R0          ! g_R0 is a rescaled version of g/Rho [Z R-1 T-2 ~> m4 kg-1 s-2].
  real :: Norm          ! A factor that normalizes two weights to 1 [H-2 ~> m-2 or m4 kg-2].
  real :: tol_dksrc     ! Tolerance for the change in the kappa source within an iteration
                        ! relative to the local source [nondim].  This must be greater than 1.
  real :: tol2          ! The tolerance for the change in the kappa source within an iteration
                        ! relative to the average local source over previous iterations [nondim].
  real :: tol_dksrc_low ! The tolerance for the fractional decrease in ksrc
                        ! within an iteration [nondim].  0 < tol_dksrc_low < 1.
  real :: Ri_crit       !   The critical shear Richardson number for shear-
                        ! driven mixing [nondim]. The theoretical value is 0.25.
  real :: dt_rem        !   The remaining time to advance the solution [T ~> s].
  real :: dt_now        !   The time step used in the current iteration [T ~> s].
  real :: dt_wt         !   The fractional weight of the current iteration [nondim].
  real :: dt_test       !   A time-step that is being tested for whether it
                        ! gives acceptably small changes in k_src [T ~> s].
  real :: Idtt          !   Idtt = 1 / dt_test [T-1 ~> s-1].
  real :: dt_inc        !   An increment to dt_test that is being tested [T ~> s].
  real :: wt_a          ! The fraction of a layer thickness identified with the interface
                        ! above a layer [nondim]
  real :: wt_b          ! The fraction of a layer thickness identified with the interface
                        ! below a layer [nondim]
  real :: k0dt          ! The background diffusivity times the timestep [H Z ~> m2 or kg m-1].
  real :: I_lz_rescale_sqr ! The inverse of a rescaling factor for L2_bdry (Lz) squared [nondim].
  logical :: valid_dt   ! If true, all levels so far exhibit acceptably small changes in k_src.
  integer :: ks_kappa, ke_kappa  ! The k-range with nonzero kappas.
  integer :: dt_refinements ! The number of 2-fold refinements that will be used
                           ! to estimate the maximum permitted time step.  I.e.,
                           ! the resolution is 1/2^dt_refinements.
  integer :: k, itt, itt_dt

  ! This calculation of N2 is for debugging only.
  ! real, dimension(SZK_(GV)+1) :: &
  !   N2_debug, & ! A version of N2 for debugging [T-2 ~> s-2]

  Ri_crit = CS%Rino_crit
  gR0 = GV%H_to_RZ * GV%g_Earth
  g_R0 = GV%g_Earth_Z_T2 / GV%Rho0
  k0dt = dt*CS%kappa_0

  I_lz_rescale_sqr = 1.0 ; if (CS%lz_rescale > 0) I_lz_rescale_sqr = 1/(CS%lz_rescale*CS%lz_rescale)

  tol_dksrc = CS%kappa_src_max_chg
  if (tol_dksrc == 10.0) then
    ! This is equivalent to the expression below, but avoids changes at roundoff for the default value.
    tol_dksrc_low = 0.95
  else
    tol_dksrc_low = (tol_dksrc - 0.5)/tol_dksrc
  endif
  tol2 = 2.0*CS%kappa_tol_err
  dt_refinements = 5 ! Selected so that 1/2^dt_refinements < 1-tol_dksrc_low

  ! Make contigous copies of 3d inputs
  do k=1,nzc
    hlay_col(k) = hlay(ii,jj,k) ; dz_lay_col(k) = dz_lay(ii,jj,k)
  enddo

  ! Set up Idz as the inverse of layer thicknesses.
  do k=1,nzc ; Idz(k) = 1.0 / dz_lay_col(k) ; enddo

  dist_from_top(1) = 0.0 ; h_from_top(1) = 0.0
  do K=2,nzc
    dist_from_top(K) = dist_from_top(K-1) + dz_lay_col(K-1)
    h_from_top(K) = h_from_top(K-1) + hlay_col(K-1)
  enddo

  ! Find the inverse of the squared distances from the boundaries.
  dist_from_bot = 0.0 ; h_from_bot = 0.0
  do K=nzc,2,-1
    dist_from_bot = dist_from_bot + dz_lay_col(K)
    h_from_bot = h_from_bot + hlay_col(K )
    ! Find the inverse of the squared distances from the boundaries,
    I_L2_bdry(K) = ((dist_from_top(K) + dist_from_bot) * (h_from_top(K) + h_from_bot)) / &
                   ((dist_from_top(K) * dist_from_bot) * (h_from_top(K) * h_from_bot))
    ! reduce the distance by a factor of "lz_rescale"
    I_L2_bdry(K) = I_lz_rescale_sqr*I_L2_bdry(K)
  enddo

  ! This uses half the harmonic mean of thicknesses to provide two estimates
  ! of the boundary between cells, and the inverse of the harmonic mean to
  ! weight the two estimates.  The net effect is that interfaces around thin
  ! layers have thin cells, and the total thickness adds up properly.
  ! The top- and bottom- interfaces have zero thickness, consistent with
  ! adding additional zero thickness layers.
  h_Int(1) = 0.0 ; h_Int(2) = hlay_col(1)
  dz_Int(1) = 0.0 ; dz_Int(2) = dz_lay_col(1)
  do K=2,nzc-1
    Norm = 1.0 / (hlay_col(K)*(hlay_col(K-1)+hlay_col(K+1)) + 2.0*hlay_col(K-1)*hlay_col(K+1))
    wt_a = ((hlay_col(K)+hlay_col(K+1)) * hlay_col(K-1)) * Norm
    wt_b = ((hlay_col(K-1)+hlay_col(K)) * hlay_col(K+1)) * Norm
    h_Int(K) = h_Int(K) + hlay_col(K) * wt_a
    h_Int(K+1) = hlay_col(K) * wt_b
    dz_Int(K) = dz_Int(K) + dz_lay_col(K) * wt_a
    dz_Int(K+1) = dz_lay_col(K) * wt_b
  enddo
  h_Int(nzc) = h_Int(nzc) + hlay_col(nzc) ; h_Int(nzc+1) = 0.0
  dz_Int(nzc) = dz_Int(nzc) + dz_lay_col(nzc) ; dz_Int(nzc+1) = 0.0

  if (GV%Boussinesq) then
    do K=1,nzc+1 ; dz_h_Int(K) = GV%H_to_Z ; enddo
  else
    ! Find an effective average specific volume around an interface.
    dz_h_Int(1:nzc+1) = 0.0
    if (hlay_col(1) > 0.0) dz_h_Int(1) = dz_lay_col(1) / hlay_col(1)
    do K=2,nzc+1
      if (h_Int(K) > 0.0) then
        dz_h_Int(K) = dz_Int(K) / h_Int(K)
      else
        dz_h_Int(K) = dz_h_Int(K-1)
      endif
    enddo
  endif

  ! N2_debug(1) = 0.0 ; N2_debug(nzc+1) = 0.0
  ! do K=2,nzc
  !   N2_debug(K) = max((dbuoy_dT(K) * (T0xdz(k-1)*Idz(k-1) - T0xdz(k)*Idz(k)) + &
  !                      dbuoy_dS(K) * (S0xdz(k-1)*Idz(k-1) - S0xdz(k)*Idz(k))) * &
  !                      I_dz_int(K), 0.0)
  ! enddo

  ! Build contiguous per-column copies of the inputs to calculate_projected_state, to
  ! avoid passing non-contiguous array sections or mismatched-rank whole arrays.
  do K=1,nzc
    u_col(K) = u(ii,jj,K) ; v_col(K) = v(ii,jj,K) ; T_col(K) = T(ii,jj,K) ; Sal_col(K) = Sal(ii,jj,K)
  enddo
  do K=1,nzc+1
    I_dz_int_col(K) = I_dz_int(ii,jj,K)
    dbuoy_dT_col(K) = dbuoy_dT(ii,jj,K) ; dbuoy_dS_col(K) = dbuoy_dS(ii,jj,K)
  enddo
  do K=1,nzc+1 ; kappa_col(K) = kappa(ii,jj,K) ; enddo

  ! This call just calculates N2 and S2.
  call calculate_projected_state(kappa_col, u_col, v_col, T_col, Sal_col, 0.0, nzc, hlay_col, I_dz_int_col, &
                                 dbuoy_dT_col, dbuoy_dS_col, CS%vel_underflow, u_col, v_col, T_col, Sal_col, &
                                 N2, S2, GV, US)
  do K=1,nzc+1
    N2_init(ii,jj,K) = N2(K)
    S2_init(ii,jj,K) = S2(K)
  enddo


! ----------------------------------------------------
! Iterate
! ----------------------------------------------------
  dt_rem = dt
  do K=1,nzc+1
    K_Q(K) = 0.0
    kappa_avg_col(K) = 0.0 ; tke_avg_col(K) = 0.0
    N2_mean_col(K) = 0.0 ; S2_mean_col(K) = 0.0
    local_src_avg(K) = 0.0
    ! Use the grid spacings to scale errors in the source.
    if ( h_Int(K) > 0.0 ) &
      local_src_avg(K) = 0.1 * k0dt * I_dz_int_col(K) / h_Int(K)
  enddo

! call cpu_clock_end(id_clock_setup)

! do itt=1,CS%max_RiNo_it
  do itt=1,CS%max_KS_it

! ----------------------------------------------------
! Calculate new values of u, v, rho, N^2 and S.
! ----------------------------------------------------

  ! call cpu_clock_begin(id_clock_KQ)
    call find_kappa_tke(N2, S2, kappa_col, Idz, h_Int, dz_Int, dz_h_Int, I_L2_bdry, f2, &
                        nzc, CS, GV, US, K_Q, tke_col, kappa_out, kappa_src, local_src)
  ! call cpu_clock_end(id_clock_KQ)

  ! call cpu_clock_begin(id_clock_avg)
    ! Determine the range of non-zero values of kappa_out.
    ks_kappa = GV%ke+1 ; ke_kappa = 0
    do K=2,nzc ; if (kappa_out(K) > 0.0) then
      ks_kappa = K ; exit
    endif ; enddo
    do k=nzc,ks_kappa,-1 ; if (kappa_out(K) > 0.0) then
      ke_kappa = K ; exit
    endif ; enddo
    if (ke_kappa == nzc) kappa_out(nzc+1) = 0.0
  ! call cpu_clock_end(id_clock_avg)

    ! Determine how long to use this value of kappa (dt_now).

  ! call cpu_clock_begin(id_clock_project)
    if ((ke_kappa < ks_kappa) .or. (itt==CS%max_KS_it)) then
      dt_now = dt_rem
    else
      ! Limit dt_now so that |k_src(k)-kappa_src(k)| < tol * local_src(k)
      dt_test = dt_rem
      do K=2,nzc
        tol_max(K) = kappa_src(K) + tol_dksrc * local_src(K)
        tol_min(K) = kappa_src(K) - tol_dksrc_low * local_src(K)
        tol_chg(K) = tol2 * local_src_avg(K)
      enddo

      do itt_dt=1,(CS%max_KS_it+1-itt)/2
        !   The maximum number of times that the time-step is halved in
        ! seeking an acceptable timestep is reduced with each iteration,
        ! so that as the maximum number of iterations is approached, the
        ! whole remaining timestep is used.  Typically, an acceptable
        ! timestep is found long before the minimum is reached, so the
        ! value of max_KS_it may be unimportant, especially if it is large
        ! enough.
        call calculate_projected_state(kappa_out, u_col, v_col, T_col, Sal_col, 0.5*dt_test, nzc, hlay_col, &
                                       I_dz_int_col, dbuoy_dT_col, dbuoy_dS_col, CS%vel_underflow, u_test, v_test, &
                                       T_test, S_test, N2, S2, GV, US, ks_int=ks_kappa, ke_int=ke_kappa)
        valid_dt = .true.
        Idtt = 1.0 / dt_test
        do K=max(ks_kappa-1,2),min(ke_kappa+1,nzc)
          if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
            K_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                       ((Ri_crit*S2(K) - N2(K)) / (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
            if (CS%restrictive_tolerance_check) then
              if ((K_src(K) > min(tol_max(K), kappa_src(K) + Idtt*tol_chg(K))) .or. &
                  (K_src(K) < max(tol_min(K), kappa_src(K) - Idtt*tol_chg(K)))) then
                valid_dt = .false. ; exit
              endif
            else
              if ((K_src(K) > max(tol_max(K), kappa_src(K) + Idtt*tol_chg(K))) .or. &
                  (K_src(K) < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K)))) then
                valid_dt = .false. ; exit
              endif
            endif
          else
            if (0.0 < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K))) then
              valid_dt = .false. ; k_src(K) = 0.0 ; exit
            endif
          endif
        enddo

        if (valid_dt) exit
        dt_test = 0.5*dt_test
      enddo
      if ((dt_test < dt_rem) .and. valid_dt) then
        dt_inc = 0.5*dt_test
        do itt_dt=1,dt_refinements
          call calculate_projected_state(kappa_out, u_col, v_col, T_col, Sal_col, 0.5*(dt_test+dt_inc), nzc, &
                   hlay_col, I_dz_int_col, dbuoy_dT_col, dbuoy_dS_col, CS%vel_underflow, u_test, v_test, &
                   T_test, S_test, N2, S2, GV, US, ks_int=ks_kappa, ke_int=ke_kappa)
          valid_dt = .true.
          Idtt = 1.0 / (dt_test+dt_inc)
          do K=max(ks_kappa-1,2),min(ke_kappa+1,nzc)
            if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
              K_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                         ((Ri_crit*S2(K) - N2(K)) / (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
              if ((K_src(K) > max(tol_max(K), kappa_src(K) + Idtt*tol_chg(K))) .or. &
                  (K_src(K) < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K)))) then
                valid_dt = .false. ; exit
              endif
            else
              if (0.0 < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K))) then
                valid_dt = .false. ; k_src(K) = 0.0 ; exit
              endif
            endif
          enddo

          if (valid_dt) dt_test = dt_test + dt_inc
          dt_inc = 0.5*dt_inc
        enddo
      else
        dt_inc = 0.0
      endif

      dt_now = min(dt_test*(1.0+CS%kappa_tol_err)+dt_inc, dt_rem)
      do K=2,nzc
        local_src_avg(K) = local_src_avg(K) + dt_now * local_src(K)
      enddo
    endif  ! Are all the values of kappa_out 0?
  ! call cpu_clock_end(id_clock_project)

    ! The state has already been projected forward. Now find new values of kappa.

    if (ke_kappa < ks_kappa) then
      ! There is no mixing now, and will not be again.
    ! call cpu_clock_begin(id_clock_avg)
      dt_wt = dt_rem / dt ; dt_rem = 0.0
      do K=1,nzc+1
        kappa_mid(K) = 0.0
        ! This would be here but does nothing.
        ! kappa_avg(K) = kappa_avg(K) + kappa_mid(K)*dt_wt
        tke_avg_col(K) = tke_avg_col(K) + dt_wt*tke_col(K)
      enddo
    ! call cpu_clock_end(id_clock_avg)
    else
    ! call cpu_clock_begin(id_clock_project)
      call calculate_projected_state(kappa_out, u_col, v_col, T_col, Sal_col, dt_now, nzc, hlay_col, I_dz_int_col, &
                                     dbuoy_dT_col, dbuoy_dS_col, CS%vel_underflow, u_test, v_test, &
                                     T_test, S_test, N2, S2, GV, US, ks_int=ks_kappa, ke_int=ke_kappa)
    ! call cpu_clock_end(id_clock_project)

    ! call cpu_clock_begin(id_clock_KQ)
      do K=1,nzc+1 ; K_Q_tmp(K) = K_Q(K) ; enddo
      call find_kappa_tke(N2, S2, kappa_out, Idz, h_Int, dz_Int, dz_h_Int, I_L2_bdry, f2, &
                          nzc, CS, GV, US, K_Q_tmp, tke_pred, kappa_pred)
    ! call cpu_clock_end(id_clock_KQ)

      ks_kappa = GV%ke+1 ; ke_kappa = 0
      do K=1,nzc+1
        kappa_mid(K) = 0.5*(kappa_out(K) + kappa_pred(K))
        if ((kappa_mid(K) > 0.0) .and. (K<ks_kappa)) ks_kappa = K
        if (kappa_mid(K) > 0.0) ke_kappa = K
      enddo

    ! call cpu_clock_begin(id_clock_project)
      call calculate_projected_state(kappa_mid, u_col, v_col, T_col, Sal_col, dt_now, nzc, hlay_col, I_dz_int_col, &
                                     dbuoy_dT_col, dbuoy_dS_col, CS%vel_underflow, u_test, v_test, &
                                     T_test, S_test, N2, S2, GV, US, ks_int=ks_kappa, ke_int=ke_kappa)
    ! call cpu_clock_end(id_clock_project)

    ! call cpu_clock_begin(id_clock_KQ)
      call find_kappa_tke(N2, S2, kappa_out, Idz, h_Int, dz_Int, dz_h_Int, I_L2_bdry, f2, &
                          nzc, CS, GV, US, K_Q, tke_pred, kappa_pred)
    ! call cpu_clock_end(id_clock_KQ)

    ! call cpu_clock_begin(id_clock_avg)
      dt_wt = dt_now / dt ; dt_rem = dt_rem - dt_now
      do K=1,nzc+1
        kappa_mid(K) = 0.5*(kappa_out(K) + kappa_pred(K))
        kappa_avg_col(K) = kappa_avg_col(K) + kappa_mid(K)*dt_wt
        tke_avg_col(K) = tke_avg_col(K) + dt_wt*0.5*(tke_pred(K) + tke_col(K))
        N2_mean_col(K) = N2_mean_col(K) + dt_wt*N2(K)
        S2_mean_col(K) = S2_mean_col(K) + dt_wt*S2(K)
        kappa_col(K) = kappa_pred(K) ! First guess for the next iteration.
      enddo

    ! call cpu_clock_end(id_clock_avg)
    endif

    if (dt_rem > 0.0) then
      ! Update the values of u, v, T, Sal, N2, and S2 for the next iteration.
    ! call cpu_clock_begin(id_clock_project)
      call calculate_projected_state(kappa_mid, u_col, v_col, T_col, Sal_col, dt_now, nzc, hlay_col, I_dz_int_col, &
                                     dbuoy_dT_col, dbuoy_dS_col, CS%vel_underflow, u_col, v_col, T_col, Sal_col, &
                                     N2, S2, GV, US)
    ! call cpu_clock_end(id_clock_project)
    endif

    if (dt_rem <= 0.0) exit

  enddo ! end itt loop

  do K=1,nzc
    u(ii,jj,K) = u_col(K) ; v(ii,jj,K) = v_col(K) ; T(ii,jj,K) = T_col(K) ; Sal(ii,jj,K) = Sal_col(K)
  enddo

  do K=1,nzc+1
    kappa(ii,jj,K) = kappa_col(K) ; tke(ii,jj,K) = tke_col(K)
    kappa_avg(ii,jj,K) = kappa_avg_col(K) ; tke_avg(ii,jj,K) = tke_avg_col(K)
    N2_mean(ii,jj,K) = N2_mean_col(K) ; S2_mean(ii,jj,K) = S2_mean_col(K)
  enddo

end subroutine kappa_shear_column

!>   This subroutine calculates the velocities, temperature and salinity that
!! the water column will have after mixing for dt with diffusivities kappa.  It
!! may also calculate the projected buoyancy frequency and shear.
pure subroutine calculate_projected_state(kappa, u0, v0, T0, S0, dt, nz, dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                     vel_under, u, v, T, Sal, N2, S2, GV, US, ks_int, ke_int)
  !$omp declare target
  integer,               intent(in)    :: nz  !< The number of layers (after eliminating massless
                                              !! layers?).
  real, dimension(nz+1), intent(in)    :: kappa !< The diapycnal diffusivity at interfaces,
                                              !! [H Z T-1 ~> m2 s-1 or Pa s].
  real, dimension(nz),   intent(in)    :: u0  !< The initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(nz),   intent(in)    :: v0  !< The initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(nz),   intent(in)    :: T0  !< The initial temperature [C ~> degC].
  real, dimension(nz),   intent(in)    :: S0  !< The initial salinity [S ~> ppt].
  real,                  intent(in)    :: dt  !< The time step [T ~> s].
  real, dimension(nz),   intent(in)    :: dz  !< The layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz+1), intent(in)    :: I_dz_int !< The inverse of the distance between successive
                                              !! layer centers [Z-1 ~> m-1].
  real, dimension(nz+1), intent(in)    :: dbuoy_dT !< The partial derivative of buoyancy with
                                              !! temperature [Z T-2 C-1 ~> m s-2 degC-1].
  real, dimension(nz+1), intent(in)    :: dbuoy_dS !< The partial derivative of buoyancy with
                                              !! salinity [Z T-2 S-1 ~> m s-2 ppt-1].
  real,                  intent(in)    :: vel_under !< Any velocities that are smaller in magnitude
                                              !! than this value are set to 0 [L T-1 ~> m s-1].
  real, dimension(nz),   intent(inout) :: u   !< The zonal velocity after dt [L T-1 ~> m s-1].
  real, dimension(nz),   intent(inout) :: v   !< The meridional velocity after dt [L T-1 ~> m s-1].
  real, dimension(nz),   intent(inout) :: T   !< The temperature after dt [C ~> degC].
  real, dimension(nz),   intent(inout) :: Sal !< The salinity after dt [S ~> ppt].
  real, dimension(nz+1), intent(inout) :: N2  !< The buoyancy frequency squared at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), intent(inout) :: S2  !< The squared shear at interfaces [T-2 ~> s-2].
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  integer, optional,     intent(in)    :: ks_int !< The topmost k-index with a non-zero diffusivity.
  integer, optional,     intent(in)    :: ke_int !< The bottommost k-index with a non-zero
                                              !! diffusivity.

  ! Local variables
  ! See the note on GPU_nk_max in kappa_shear_column above.
#ifdef __NVCOMPILER_OPENMP_GPU
  real, dimension(GPU_nk_max+1) :: c1 ! A tridiagonal variable [nondim]
#else
  real, dimension(nz+1) :: c1 ! A tridiagonal variable [nondim]
#endif
  real :: a_a, a_b   ! Tridiagonal coupling coefficients [H ~> m or kg m-2]
  real :: b1, b1nz_0 ! Tridiagonal variables [H-1 ~> m-1 or m2 kg-1]
  real :: bd1        ! A term in the denominator of b1 [H ~> m or kg m-2]
  real :: d1         ! A tridiagonal variable [nondim]
  integer :: k, ks, ke

  ks = 1 ; ke = nz
  if (present(ks_int)) ks = max(ks_int-1,1)
  if (present(ke_int)) ke = min(ke_int,nz)

  if (ks > ke) return

  if (dt > 0.0) then
    a_b = dt*(kappa(ks+1)*I_dz_int(ks+1))
    b1 = 1.0 / (dz(ks) + a_b)
    c1(ks+1) = a_b * b1 ; d1 = dz(ks) * b1 ! = 1 - c1

    u(ks) = (b1 * dz(ks))*u0(ks) ; v(ks) = (b1 * dz(ks))*v0(ks)
    T(ks) = (b1 * dz(ks))*T0(ks) ; Sal(ks) = (b1 * dz(ks))*S0(ks)
    do K=ks+1,ke-1
      a_a = a_b
      a_b = dt*(kappa(K+1)*I_dz_int(K+1))
      bd1 = dz(k) + d1*a_a
      b1 = 1.0 / (bd1 + a_b)
      c1(K+1) = a_b * b1 ; d1 = bd1 * b1 ! d1 = 1 - c1

      u(k) = b1 * (dz(k)*u0(k) + a_a*u(k-1))
      v(k) = b1 * (dz(k)*v0(k) + a_a*v(k-1))
      T(k) = b1 * (dz(k)*T0(k) + a_a*T(k-1))
      Sal(k) = b1 * (dz(k)*S0(k) + a_a*Sal(k-1))
    enddo
    !   T and S have insulating boundary conditions, u & v use no-slip
    ! bottom boundary conditions at the solid bottom.

    ! For insulating boundary conditions or mixing simply stopping, use...
    a_a = a_b
    b1 = 1.0 / (dz(ke) + d1*a_a)
    T(ke) = b1 * (dz(ke)*T0(ke) + a_a*T(ke-1))
    Sal(ke) = b1 * (dz(ke)*S0(ke) + a_a*Sal(ke-1))

    !   There is no distinction between the effective boundary conditions for
    ! tracers and velocities if the mixing is separated from the bottom, but if
    ! the mixing goes all the way to the bottom, use no-slip BCs for velocities.
    if (ke == nz) then
      a_b = dt*(kappa(nz+1)*I_dz_int(nz+1))
      b1nz_0 = 1.0 / ((dz(nz) + d1*a_a) + a_b)
    else
      b1nz_0 = b1
    endif
    u(ke) = b1nz_0 * (dz(ke)*u0(ke) + a_a*u(ke-1))
    v(ke) = b1nz_0 * (dz(ke)*v0(ke) + a_a*v(ke-1))
    if (abs(u(ke)) < vel_under) u(ke) = 0.0
    if (abs(v(ke)) < vel_under) v(ke) = 0.0

    do k=ke-1,ks,-1
      u(k) = u(k) + c1(k+1)*u(k+1)
      v(k) = v(k) + c1(k+1)*v(k+1)
      if (abs(u(k)) < vel_under) u(k) = 0.0
      if (abs(v(k)) < vel_under) v(k) = 0.0
      T(k) = T(k) + c1(k+1)*T(k+1)
      Sal(k) = Sal(k) + c1(k+1)*Sal(k+1)
    enddo
  else ! dt <= 0.0
    do k=1,nz
      u(k) = u0(k) ; v(k) = v0(k) ; T(k) = T0(k) ; Sal(k) = S0(k)
      if (abs(u(k)) < vel_under) u(k) = 0.0
      if (abs(v(k)) < vel_under) v(k) = 0.0
    enddo
  endif

  ! Store the squared shear at interfaces
  S2(1) = 0.0 ; S2(nz+1) = 0.0
  if (ks > 1) &
    S2(ks) = (((u(ks)-u0(ks-1))**2) + ((v(ks)-v0(ks-1))**2)) * (US%L_to_Z*I_dz_int(ks))**2
  do K=ks+1,ke
    S2(K) = (((u(k)-u(k-1))**2) + ((v(k)-v(k-1))**2)) * (US%L_to_Z*I_dz_int(K))**2
  enddo
  if (ke<nz) &
    S2(ke+1) = (((u0(ke+1)-u(ke))**2) + ((v0(ke+1)-v(ke))**2)) * (US%L_to_Z*I_dz_int(ke+1))**2

  ! Store the buoyancy frequency at interfaces
  N2(1) = 0.0 ; N2(nz+1) = 0.0
  if (ks > 1) &
    N2(ks) = max(0.0, I_dz_int(ks) * &
      (dbuoy_dT(ks) * (T0(ks-1)-T(ks)) + dbuoy_dS(ks) * (S0(ks-1)-Sal(ks))))
  do K=ks+1,ke
    N2(K) = max(0.0, I_dz_int(K) * &
      (dbuoy_dT(K) * (T(k-1)-T(k)) + dbuoy_dS(K) * (Sal(k-1)-Sal(k))))
  enddo
  if (ke<nz) &
    N2(ke+1) = max(0.0, I_dz_int(ke+1) * &
      (dbuoy_dT(ke+1) * (T(ke)-T0(ke+1)) + dbuoy_dS(ke+1) * (Sal(ke)-S0(ke+1))))

end subroutine calculate_projected_state

!> This subroutine calculates new, consistent estimates of TKE and kappa.
pure subroutine find_kappa_tke(N2, S2, kappa_in, Idz, h_Int, dz_Int, dz_h_Int, I_L2_bdry, f2, &
                          nz, CS, GV, US, K_Q, tke, kappa, kappa_src, local_src)
  !$omp declare target
  integer,               intent(in)    :: nz  !< The number of layers to work on.
  real, dimension(nz+1), intent(in)    :: N2  !< The buoyancy frequency squared at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), intent(in)    :: S2  !< The squared shear at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), intent(in)    :: kappa_in  !< The initial guess at the diffusivity
                                              !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(nz+1), intent(in)    :: h_Int !< The thicknesses associated with interfaces
                                              !! [H ~> m or kg m-2]
  real, dimension(nz+1), intent(in)    :: dz_Int !< The vertical distances around interfaces [Z ~> m]
  real, dimension(nz+1), intent(in)    :: dz_h_Int !< The ratio of the vertical distances to the
                                              !! thickness around an interface [Z H-1 ~> nondim or m3 kg-1].
                                              !! In non-Boussinesq mode this is the specific volume.
  real, dimension(nz+1), intent(in)    :: I_L2_bdry !< The inverse of the squared distance to
                                              !! boundaries [H-1 Z-1 ~> m-2 or m kg-1].
  real, dimension(nz),   intent(in)    :: Idz !< The inverse grid spacing of layers [Z-1 ~> m-1].
  real,                  intent(in)    :: f2  !< The squared Coriolis parameter [T-2 ~> s-2].
  type(Kappa_shear_CS),  intent(in)    :: CS  !< This module's control structure.  Not a pointer:
                                              !! see the note in kappa_shear_column.
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(nz+1), intent(inout) :: K_Q !< The shear-driven diapycnal diffusivity divided by
                                              !! the turbulent kinetic energy per unit mass at
                                              !! interfaces [H T Z-1 ~> s or kg s m-3].
  real, dimension(nz+1), intent(out)   :: tke !< The turbulent kinetic energy per unit mass at
                                              !! interfaces [Z2 T-2 ~> m2 s-2].
  real, dimension(nz+1), intent(out)   :: kappa !< The diapycnal diffusivity at interfaces
                                              !! [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(nz+1), optional, &
                         intent(out)   :: kappa_src !< The source term for kappa [T-1 ~> s-1]
  real, dimension(nz+1), optional, &
                         intent(out)   :: local_src !< The sum of all local sources for kappa
                                              !! [T-1 ~> s-1]
  ! This subroutine calculates new, consistent estimates of TKE and kappa.

  ! Local variables
  ! See the note on GPU_nk_max in kappa_shear_column above.
#ifdef __NVCOMPILER_OPENMP_GPU
  real, dimension(GPU_nk_max) :: &
#else
  real, dimension(nz) :: &
#endif
    aQ, &       ! aQ is the coupling between adjacent interfaces in the TKE equations [H T-1 ~> m s-1 or kg m-2 s-1]
    dQdz        ! Half the partial derivative of TKE with depth [Z T-2 ~> m s-2].
#ifdef __NVCOMPILER_OPENMP_GPU
  real, dimension(GPU_nk_max+1) :: &
#else
  real, dimension(nz+1) :: &
#endif
    dK, &         ! The change in kappa [H Z T-1 ~> m2 s-1 or Pa s].
    dQ, &         ! The change in TKE [Z2 T-2 ~> m2 s-2].
    cQ, cK, &     ! cQ and cK are the upward influences in the tridiagonal and
                  ! hexadiagonal solvers for the TKE and kappa equations [nondim].
    I_Ld2, &      ! 1/Ld^2, where Ld is the effective decay length scale for kappa [H-1 Z-1 ~> m-2 or m kg-1]
    TKE_decay, &  ! The local TKE decay rate [T-1 ~> s-1].
    k_src, &      ! The source term in the kappa equation [T-1 ~> s-1].
    dQmdK, &      ! With Newton's method the change in dQ(k-1) due to dK(k) [Z T H-1 ~> s or m3 s kg-1]
    dKdQ, &       ! With Newton's method the change in dK(k) due to dQ(k) [H Z-1 T-1 ~> s-1 or kg m-3 s-1]
    e1            ! The fractional change in a layer TKE due to a change in the
                  ! TKE of the layer above when all the kappas below are 0 [nondim].
                  ! e1 is nondimensional, and 0 < e1 < 1.
  real :: tke_src       ! The net source of TKE due to mixing against the shear and stratification
                        ! [Z2 T-3 ~> m2 s-3] or [H Z T-3 ~> m2 s-3 or kg m-1 s-3].
                        ! (For convenience, a term involving the non-dissipation of q0 is also included here.)
  real :: bQ            ! The inverse of the pivot in the tridiagonal equations [T H-1 ~> s m-1 or m2 s kg-1]
  real :: bK            ! The inverse of the pivot in the tridiagonal equations [Z-1 ~> m-1].
  real :: bQd1          ! A term in the denominator of bQ [H T-1 ~> m s-1 or kg m-2 s-1]
  real :: bKd1          ! A term in the denominator of bK [Z ~> m].
  real :: cQcomp, cKcomp ! 1 - cQ or 1 - cK in the tridiagonal equations [nondim].
  real :: c_s2          !   The coefficient for the decay of TKE due to
                        ! shear (i.e. proportional to |S|*tke) [nondim].
  real :: c_n2          !   The coefficient for the decay of TKE due to
                        ! stratification (i.e. proportional to N*tke) [nondim].
  real :: Ri_crit       !   The critical shear Richardson number for shear-
                        ! driven mixing [nondim]. The theoretical value is 0.25.
  real :: q0            !   The background level of TKE [Z2 T-2 ~> m2 s-2].
  real :: Ilambda2      ! 1.0 / CS%lambda**2 [nondim]
  real :: TKE_min       !   The minimum value of shear-driven TKE that can be
                        ! solved for [Z2 T-2 ~> m2 s-2].
  real :: kappa0        ! The background diapycnal diffusivity [H Z T-1 ~> m2 s-1 or Pa s]
  real :: kappa_trunc   ! Diffusivities smaller than this are rounded to 0 [H Z T-1 ~> m2 s-1 or Pa s]

  real :: eden1, eden2  ! Variables used in calculating e1 [H Z-2 ~> m-1 or kg m-4]
  real :: I_eden        ! The inverse of the denominator in e1 [Z2 H-1 ~> m or m4 kg-1]
  real :: ome           ! Variables used in calculating e1 [nondim]
  real :: diffusive_src ! The diffusive source in the kappa equation [H T-1 ~> m s-1 or kg m-2 s-1]
  real :: chg_by_k0     ! The value of k_src that leads to an increase of
                        ! kappa_0 if only the diffusive term is a sink [T-1 ~> s-1]
  real :: h_dz_here     ! The ratio of the thicknesses to the vertical distances around an interface
                        ! [H Z-1 ~> nondim or kg m-3].  In non-Boussinesq mode this is the density.

  real :: kappa_mean    ! A mean value of kappa [H Z T-1 ~> m2 s-1 or Pa s]
  real :: Newton_test   ! The value of relative error that will cause the next
                        ! iteration to use Newton's method [nondim].
  ! Temporary variables used in the Newton's method iterations.
  real :: decay_term_k  ! The decay term in the diffusivity equation [Z-1 ~> m-1]
  real :: decay_term_Q  ! The decay term in the TKE equation - proportional to [H Z-1 T-1 ~> s-1 or kg m-3 s-1]
  real :: I_Q           ! The inverse of TKE [T2 Z-2 ~> s2 m-2]
  real :: kap_src       ! A source term in the kappa equation [H T-1 ~> m s-1 or kg m-2 s-1]
  real :: v1            ! A temporary variable proportional to [H Z-1 T-1 ~> s-1 or kg m-3 s-1]
  real :: v2            ! A temporary variable in [Z T-2 ~> m s-2]
  real :: tol_err       ! The tolerance for max_err that determines when to
                        ! stop iterating [nondim].
  real :: Newton_err    ! The tolerance for max_err that determines when to
                        ! start using Newton's method [nondim].  Empirically, an initial
                        ! value of about 0.2 seems to be most efficient.
  real, parameter :: roundoff = 1.0e-16 ! A negligible fractional change in TKE [nondim].
                        ! This could be larger but performance gains are small.

  logical, parameter :: tke_noflux_bottom_BC = .false. ! Specify the boundary conditions
  logical, parameter :: tke_noflux_top_BC = .false.    ! that are applied to the TKE equations.
  logical :: do_Newton    ! If .true., use Newton's method for the next iteration.
  logical :: abort_Newton ! If .true., an Newton's method has encountered a 0
                          ! pivot, and should not have been used.
  logical :: was_Newton   ! The value of do_Newton before checking convergence.
  logical :: within_tolerance ! If .true., all points are within tolerance to
                          ! enable this subroutine to return.
  integer :: ks_src, ke_src ! The range indices that have nonzero k_src.
  integer :: ks_kappa, ke_kappa, ke_tke   ! The ranges of k-indices that are or
  integer :: ks_kappa_prev, ke_kappa_prev ! were being worked on.
  integer :: itt, k, k2

  ! These variables are used only for debugging.
  logical, parameter :: debug_soln = .false.
  real :: K_err_lin ! The imbalance in the K equation [H T-1 ~> m s-1 or kg m-2 s-1]
  real :: Q_err_lin ! The imbalance in the Q equation [H Z T-3 ~> m2 s-3 or kg m-1 s-3]
  real, dimension(nz+1) :: &
    I_Ld2_debug, & ! A separate version of I_Ld2 for debugging [H-1 Z-1 ~> m-2 or m kg-1].
    kappa_prev, & ! The value of kappa at the start of the current iteration [H Z T-1 ~> m2 s-1 or Pa s]
    TKE_prev   ! The value of TKE at the start of the current iteration [Z2 T-2 ~> m2 s-2].

  c_N2 = CS%C_N**2 ; c_S2 = CS%C_S**2
  q0 = CS%TKE_bg ; kappa0 = CS%kappa_0
  TKE_min = max(CS%TKE_bg, 1.0E-20*US%m_to_Z**2*US%T_to_s**2)
  Ri_crit = CS%Rino_crit
  Ilambda2 = 1.0 / CS%lambda**2
  kappa_trunc = CS%kappa_trunc
  do_Newton = .false. ; abort_Newton = .false.
  tol_err = CS%kappa_tol_err
  Newton_err = 0.2     ! This initial value may be automatically reduced later.

  ks_kappa = 2 ; ke_kappa = nz ; ks_kappa_prev = 2 ; ke_kappa_prev = nz

  ke_src = 0 ; ks_src = nz+1
  do K=2,nz
    if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
!       Ri = N2(K) / S2(K)
!       k_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
!                  ((Ri_crit - Ri) / (Ri_crit + CS%FRi_curvature*Ri))
      K_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                 ((Ri_crit*S2(K) - N2(K)) / (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
      ke_src = K
      if (ks_src > k) ks_src = K
    else
      k_src(K) = 0.0
    endif
  enddo

  ! If there is no source anywhere, return kappa(K) = 0.
  if (ks_src > ke_src) then
    do K=1,nz+1
      kappa(K) = 0.0 ; K_Q(K) = 0.0 ; tke(K) = TKE_min
    enddo
    if (present(kappa_src)) then ; do K=1,nz+1 ; kappa_src(K) = 0.0 ; enddo ; endif
    if (present(local_src)) then ; do K=1,nz+1 ; local_src(K) = 0.0 ; enddo ; endif
    return
  endif

  do K=1,nz+1
    kappa(K) = kappa_in(K)
!     TKE_decay(K) = c_n*sqrt(N2(K)) + c_s*sqrt(S2(K)) ! The expression in JHL.
    TKE_decay(K) = sqrt(c_n2*N2(K) + c_s2*S2(K))
    if ((kappa(K) > 0.0) .and. (K_Q(K) > 0.0)) then
      TKE(K) = kappa(K) / K_Q(K) ! Perhaps take the max with TKE_min
    else
      TKE(K) = TKE_min
    endif
  enddo
  ! Apply boundary conditions to kappa.
  kappa(1) = 0.0 ; kappa(nz+1) = 0.0

  ! Calculate the term (e1) that allows changes in TKE to be calculated quickly
  ! below the deepest nonzero value of kappa.  If kappa = 0, below interface
  ! k-1, the final changes in TKE are related by dQ(K+1) = e1(K+1)*dQ(K).
  eden2 = kappa0 * Idz(nz)
  if (tke_noflux_bottom_BC) then
    eden1 = h_Int(nz+1)*TKE_decay(nz+1)
    I_eden = 1.0 / (eden2 + eden1)
    e1(nz+1) = eden2 * I_eden ; ome = eden1 * I_eden
  else
    e1(nz+1) = 0.0 ; ome = 1.0
  endif
  do k=nz,2,-1
    eden1 = h_Int(K)*TKE_decay(K) + ome * eden2
    eden2 = kappa0 * Idz(k-1)
    I_eden = 1.0 / (eden2 + eden1)
    e1(K) = eden2 * I_eden ; ome = eden1 * I_eden ! = 1-e1
  enddo
  e1(1) = 0.0


  ! Iterate here to convergence to within some tolerance of order tol_err.
  do itt=1,CS%max_RiNo_it

  ! ----------------------------------------------------
  ! Calculate TKE
  ! ----------------------------------------------------

    if (debug_soln) then ; do K=1,nz+1 ; kappa_prev(K) = kappa(K) ; TKE_prev(K) = TKE(K) ; enddo ; endif

    if (.not.do_Newton) then
      !   Use separate steps of the TKE and kappa equations, that are
      ! explicit in the nonlinear source terms, implicit in a linearized
      ! version of the nonlinear sink terms, and implicit in the linear
      ! terms.

      ke_tke = max(ke_kappa,ke_kappa_prev)+1
      ! aQ is the coupling between adjacent interfaces [Z T-1 ~> m s-1].
      do k=1,min(ke_tke,nz)
        aQ(k) = (0.5*(kappa(K)+kappa(K+1)) + kappa0) * Idz(k)
      enddo
      dQ(1) = -TKE(1)
      if (tke_noflux_top_BC) then
        tke_src = dz_h_Int(1)*kappa0*S2(1) + q0 * TKE_decay(1) ! Uses that kappa(1) = 0
        bQd1 = h_Int(1) * TKE_decay(1)
        bQ = 1.0 / (bQd1 +  aQ(1))
        tke(1) = bQ * (h_Int(1)*tke_src)
        cQ(2) = aQ(1) * bQ ; cQcomp = bQd1 * bQ ! = 1 - cQ
      else
        tke(1) = q0 ; cQ(2) = 0.0 ; cQcomp = 1.0
      endif
      do K=2,ke_tke-1
        dQ(K) = -TKE(K)
        tke_src = dz_h_Int(K)*(kappa(K) + kappa0)*S2(K) + q0*TKE_decay(K)
        bQd1 = h_Int(K)*(TKE_decay(K) + dz_h_Int(K)*N2(K)*K_Q(K)) + cQcomp*aQ(k-1)
        bQ = 1.0 / (bQd1 + aQ(k))
        tke(K) = bQ * (h_Int(K)*tke_src + aQ(k-1)*tke(K-1))
        cQ(K+1) = aQ(k) * bQ ; cQcomp = bQd1 * bQ ! = 1 - cQ
      enddo
      if ((ke_tke == nz+1) .and. .not.(tke_noflux_bottom_BC)) then
        tke(nz+1) = TKE_min
        dQ(nz+1) = 0.0
      else
        k = ke_tke
        tke_src = dz_h_Int(K)*kappa0*S2(K) + q0*TKE_decay(K) ! Uses that kappa(ke_tke) = 0
        if (K == nz+1) then
          dQ(K) = -TKE(K)
          bQ = 1.0 / (h_Int(K)*TKE_decay(K) + cQcomp*aQ(k-1))
          tke(K) = max(TKE_min, bQ * (h_Int(K)*tke_src + aQ(k-1)*tke(K-1)))
          dQ(K) = tke(K) + dQ(K)
        else
          bQ = 1.0 / ((h_Int(K)*TKE_decay(K) + cQcomp*aQ(k-1)) + aQ(k))
          cQ(K+1) = aQ(k) * bQ
          ! Account for all changes deeper in the water column.
          dQ(K) = -TKE(K)
          tke(K) = max((bQ * (h_Int(K)*tke_src + aQ(k-1)*tke(K-1)) + &
                        cQ(K+1)*(tke(K+1) - e1(K+1)*tke(K))) / (1.0 - cQ(K+1)*e1(K+1)), TKE_min)
          dQ(K) = tke(K) + dQ(K)

          ! Adjust TKE deeper in the water column in case ke_tke increases.
          ! This might not be strictly necessary?
          do K=ke_tke+1,nz+1
            dQ(K) = e1(K)*dQ(K-1)
            tke(K) = max(tke(K) + dQ(K), TKE_min)
            if (abs(dQ(K)) < roundoff*tke(K)) exit
          enddo
          do K2=K+1,nz
            if (dQ(K2) == 0.0) exit
            dQ(K2) = 0.0
          enddo
        endif
      endif
      do K=ke_tke-1,1,-1
        tke(K) = max(tke(K) + cQ(K+1)*tke(K+1), TKE_min)
        dQ(K) = tke(K) + dQ(K)
      enddo

  ! ----------------------------------------------------
  ! Calculate kappa, here defined at interfaces.
  ! ----------------------------------------------------

      ke_kappa_prev = ke_kappa ; ks_kappa_prev = ks_kappa

      dK(1) = 0.0 ! kappa takes boundary values of 0.
      cK(2) = 0.0 ; cKcomp = 1.0
      if (itt == 1) then ; do K=2,nz
        I_Ld2(K) = dz_h_Int(K)*(N2(K)*Ilambda2 + f2) / tke(K) + I_L2_bdry(K)
      enddo ; endif
      do K=2,nz
        dK(K) = -kappa(K)
        if (itt>1) &
          I_Ld2(K) = dz_h_Int(K)*(N2(K)*Ilambda2 + f2) / tke(K) + I_L2_bdry(K)
        bKd1 = h_Int(K)*I_Ld2(K) + cKcomp*Idz(k-1)
        bK = 1.0 / (bKd1 + Idz(k))

        kappa(K) = bK * (Idz(k-1)*kappa(K-1) + h_Int(K) * K_src(K))
        cK(K+1) = Idz(k) * bK ; cKcomp = bKd1 * bK ! = 1 - cK(K+1)

        ! Neglect values that are smaller than kappa_trunc.
        if (kappa(K) < cKcomp*kappa_trunc) then
          kappa(K) = 0.0
          if (K > ke_src) then ; ke_kappa = k-1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (kappa(K) < 2.0*cKcomp*kappa_trunc) then
          kappa(K) = 2.0 * (kappa(K) - cKcomp*kappa_trunc)
        endif
      enddo
      K_Q(ke_kappa) = kappa(ke_kappa) / tke(ke_kappa)
      dK(ke_kappa) = dK(ke_kappa) + kappa(ke_kappa)
      do K=ke_kappa+2,ke_kappa_prev
        dK(K) = -kappa(K) ; kappa(K) = 0.0 ; K_Q(K) = 0.0
      enddo
      do K=ke_kappa-1,2,-1
        kappa(K) = kappa(K) + cK(K+1)*kappa(K+1)
        ! Neglect values that are smaller than kappa_trunc.
        if (kappa(K) <= kappa_trunc) then
          kappa(K) = 0.0
          if (K < ks_src) then ; ks_kappa = k+1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (kappa(K) < 2.0*kappa_trunc) then
          kappa(K) = 2.0 * (kappa(K) - kappa_trunc)
        endif

        dK(K) = dK(K) + kappa(K)
        K_Q(K) = kappa(K) / tke(K)
      enddo
      do K=ks_kappa_prev,ks_kappa-2 ; kappa(K) = 0.0 ; K_Q(K) = 0.0 ; enddo

    else ! do_Newton is .true.
!   Once the solutions are close enough, use a Newton's method solver of the
!  whole system to accelerate convergence.
      ks_kappa_prev = ks_kappa ; ke_kappa_prev = ke_kappa ; ke_kappa = nz
      ks_kappa = 2
      dK(1) = 0.0 ; cK(2) = 0.0 ; cKcomp = 1.0 ; dKdQ(1) = 0.0
      aQ(1) = (0.5*(kappa(1)+kappa(2))+kappa0) * Idz(1)
      dQdz(1) = 0.5*(TKE(1) - TKE(2))*Idz(1)
      if (tke_noflux_top_BC) then
        tke_src = h_Int(1) * (kappa0*dz_h_Int(1)*S2(1) - (TKE(1) - q0)*TKE_decay(1)) - &
                  aQ(1) * (TKE(1) - TKE(2))

        bQ = 1.0 / (aQ(1) + h_Int(1)*TKE_decay(1))
        cQ(2) = aQ(1) * bQ
        cQcomp = (h_Int(1)*TKE_decay(1)) * bQ ! = 1 - cQ(2)
        dQmdK(2) = -dQdz(1) * bQ
        dQ(1) = bQ * tke_src
      else
        dQ(1) = 0.0 ; cQ(2) = 0.0 ; cQcomp = 1.0 ; dQmdK(2) = 0.0
      endif
      do K=2,nz
        I_Q = 1.0 / TKE(K)
        I_Ld2(K) = (N2(K)*Ilambda2 + f2) * dz_h_Int(K)*I_Q + I_L2_bdry(K)

        kap_src = h_Int(K) * (K_src(K) - I_Ld2(K)*kappa(K)) + &
                            Idz(k-1)*(kappa(K-1)-kappa(K)) - Idz(k)*(kappa(K)-kappa(K+1))

        ! Ensure that the pivot is always positive, and that 0 <= cK <= 1.
        ! Otherwise do not use Newton's method.
        decay_term_k = -Idz(k-1)*dQmdK(K)*dKdQ(K-1) + h_Int(K)*I_Ld2(K)
        if (decay_term_k < 0.0) then ; abort_Newton = .true. ; exit ; endif
        bK = 1.0 / (Idz(k) + Idz(k-1)*cKcomp + decay_term_k)

        cK(K+1) = bK * Idz(k)
        cKcomp = bK * (Idz(k-1)*cKcomp + decay_term_k) ! = 1-cK(K+1)
        if (CS%dKdQ_iteration_bug) then
          dKdQ(K) = bK * (Idz(k-1)*dKdQ(K-1)*cQ(K) + &
                      US%m_to_Z*(N2(K)*Ilambda2 + f2) * I_Q**2 * kappa(K) )
        else
          dKdQ(K) = bK * (Idz(k-1)*dKdQ(K-1)*cQ(K) + &
                      dz_Int(K)*(N2(K)*Ilambda2 + f2) * I_Q**2 * kappa(K) )
        endif
        dK(K) = bK * (kap_src + Idz(k-1)*dK(K-1) + Idz(k-1)*dKdQ(K-1)*dQ(K-1))

        ! Truncate away negligibly small values of kappa.
        if (dK(K) <= cKcomp*(kappa_trunc - kappa(K))) then
          dK(K) = -cKcomp*kappa(K)
!         if (K > ke_src) then ; ke_kappa = k-1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (dK(K) < cKcomp*(2.0*kappa_trunc - kappa(K))) then
          dK(K) = 2.0 * dK(K) - cKcomp*(2.0*kappa_trunc - kappa(K))
        endif

        ! Solve for dQ(K)...
        aQ(k) = (0.5*(kappa(K)+kappa(K+1))+kappa0) * Idz(k)
        dQdz(k) = 0.5*(TKE(K) - TKE(K+1))*Idz(k)
        tke_src = h_Int(K) * ((dz_h_Int(K) * ((kappa(K) + kappa0)*S2(K) - kappa(k)*N2(K))) - &
                               (TKE(k) - q0)*TKE_decay(k)) - &
                  (aQ(k) * (TKE(K) - TKE(K+1)) - aQ(k-1) * (TKE(K-1) - TKE(K)))
        v1 = aQ(k-1) + dQdz(k-1)*dKdQ(K-1)
        v2 = (v1*dQmdK(K) + dQdz(k-1)*cK(K)) + &
             ((dQdz(k-1) - dQdz(k)) + dz_Int(K)*(S2(K) - N2(K)))

        ! Ensure that the pivot is always positive, and that 0 <= cQ <= 1.
        ! Otherwise do not use Newton's method.
        decay_term_Q = h_Int(K)*TKE_decay(K) - dQdz(k-1)*dKdQ(K-1)*cQ(K) - v2*dKdQ(K)
        if (decay_term_Q < 0.0) then ; abort_Newton = .true. ; exit ; endif
        bQ = 1.0 / (aQ(k) + (cQcomp*aQ(k-1) + decay_term_Q))

        cQ(K+1) = aQ(k) * bQ
        cQcomp = (cQcomp*aQ(k-1) + decay_term_Q) * bQ
        dQmdK(K+1) = (v2 * cK(K+1) - dQdz(k)) * bQ

        ! Ensure that TKE+dQ will not drop below 0.5*TKE.
        dQ(K) = max(bQ * ((v1 * dQ(K-1) + dQdz(k-1)*dK(k-1)) + &
                          (v2 * dK(K) + tke_src)), cQcomp*(-0.5*TKE(K)))

        ! Check whether the next layer will be affected by any nonzero kappas.
        if ((itt > 1) .and. (K > ke_src) .and. (dK(K) == 0.0) .and. &
            ((kappa(K) + kappa(K+1)) == 0.0)) then
        ! Could also do  .and. (bQ*abs(tke_src) < roundoff*TKE(K)) then
          ke_kappa = k-1 ; exit
        endif
      enddo
      if ((ke_kappa == nz) .and. (.not. abort_Newton)) then
        dK(nz+1) = 0.0 ; dKdQ(nz+1) = 0.0
        if (tke_noflux_bottom_BC) then
          K = nz+1
          tke_src = h_Int(K) * (kappa0*dz_h_Int(K)*S2(K) - (TKE(K) - q0)*TKE_decay(K)) + &
                    aQ(k-1) * (TKE(K-1) - TKE(K))

          v1 = aQ(k-1) + dQdz(k-1)*dKdQ(K-1)
          decay_term_Q = max(0.0, h_Int(K)*TKE_decay(K) - dQdz(k-1)*dKdQ(K-1)*cQ(K))
          if (decay_term_Q < 0.0) then
            abort_Newton = .true.
          else
            bQ = 1.0 / (aQ(k) + (cQcomp*aQ(k-1) + decay_term_Q))
          ! Ensure that TKE+dQ will not drop below 0.5*TKE.
            dQ(K) = max(bQ * ((v1 * dQ(K-1) + dQdz(k-1)*dK(K-1)) + tke_src), -0.5*TKE(K))
            TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          endif
        else
          dQ(nz+1) = 0.0
        endif
      elseif (.not. abort_Newton) then
        ! Alter the first-guess determination of dQ(K).
        dQ(ke_kappa+1) = dQ(ke_kappa+1) / (1.0 - cQ(ke_kappa+2)*e1(ke_kappa+2))
        TKE(ke_kappa+1) = max(TKE(ke_kappa+1) + dQ(ke_kappa+1), TKE_min)
        do k=ke_kappa+2,nz+1
          if (debug_soln .and. (K < nz+1)) then
          ! Ignore this source?
            aQ(k) = (0.5*(kappa(K)+kappa(K+1))+kappa0) * Idz(k)
        !    tke_src_norm = ((kappa0*dz_Int(K)*S2(K) - h_Int(K)*(TKE(K)-q0)*TKE_decay(K)) - &
        !                   (aQ(k) * (TKE(K) - TKE(K+1)) - aQ(k-1) * (TKE(K-1) - TKE(K))) ) / &
        !                   (aQ(k) + (aQ(k-1) + h_Int(K)*TKE_decay(K)))
          endif
          dK(K) = 0.0
        ! Ensure that TKE+dQ will not drop below 0.5*TKE.
          dQ(K) = max(e1(K)*dQ(K-1),-0.5*TKE(K))
          TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          if (abs(dQ(K)) < roundoff*TKE(K)) exit
        enddo
        if (debug_soln) then ; do K2=K+1,nz+1 ; dQ(K2) = 0.0 ; dK(K2) = 0.0 ; enddo ; endif
      endif
      if (.not. abort_Newton) then
        do K=ke_kappa,2,-1
          ! Ensure that TKE+dQ will not drop below 0.5*TKE.
          dQ(K) = max(dQ(K) + (cQ(K+1)*dQ(K+1) + dQmdK(K+1) * dK(K+1)), -0.5*TKE(K))
          TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          dK(K) = dK(K) + (cK(K+1)*dK(K+1) + dKdQ(K) * dQ(K))
          ! Truncate away negligibly small values of kappa.
          if (dK(K) <= kappa_trunc - kappa(K)) then
            dK(K) = -kappa(K)
            kappa(K) = 0.0
            if ((K < ks_src) .and. (K+1 > ks_kappa)) ks_kappa = K+1
          elseif (dK(K) < 2.0*kappa_trunc - kappa(K)) then
            dK(K) =  2.0*dK(K) - (2.0*kappa_trunc - kappa(K))
            kappa(K) = max(kappa(K) + dK(K), 0.0) ! The max is for paranoia.
            if (K<=ks_kappa) ks_kappa = 2
          else
            kappa(K) = kappa(K) + dK(K)
            if (K<=ks_kappa) ks_kappa = 2
          endif
        enddo
        dQ(1) = max(dQ(1) + cQ(2)*dQ(2) + dQmdK(2) * dK(2), TKE_min - TKE(1))
        TKE(1) = max(TKE(1) + dQ(1), TKE_min)
        dK(1) = 0.0
      endif

      ! Check these solutions for consistency.
      !  The unit conversions here have not been carefully tested.
      if (debug_soln) then ; do K=2,nz
        ! In these equations, K_err_lin and Q_err_lin should be at round-off levels
        ! compared with the dominant terms, perhaps, h_Int*I_Ld2*kappa and
        ! h_Int*TKE_decay*TKE.  The exception is where, either 1) the decay term has been
        ! been increased to ensure a positive pivot, or 2) negative TKEs have been
        ! truncated, or 3) small or negative kappas have been rounded toward 0.
        I_Q = 1.0 / TKE(K)
        I_Ld2_debug(K) = (N2(K)*Ilambda2 + f2) * dz_h_Int(K)*I_Q + I_L2_bdry(K)

        kap_src = h_Int(K) * (k_src(K) - I_Ld2(K)*kappa_prev(K)) + &
                            (Idz(k-1)*(kappa_prev(k-1)-kappa_prev(k)) - &
                             Idz(k)*(kappa_prev(k)-kappa_prev(k+1)))
        K_err_lin = -Idz(k-1)*(dK(K-1)-dK(K)) + Idz(k)*(dK(K)-dK(K+1)) + &
                     h_Int(K)*I_Ld2_debug(K)*dK(K) - kap_src - &
                     dz_Int(K)*(N2(K)*Ilambda2 + f2)*I_Q**2*kappa_prev(K) * dQ(K)

        h_dz_here = 0.0 ; if (abs(dz_h_Int(K)) > 0.0) h_dz_here = 1.0 / dz_h_Int(K)
        tke_src = h_Int(K) * ((kappa_prev(K) + kappa0)*S2(K) - &
                     kappa_prev(K)*N2(K) - (TKE_prev(K) - q0)*h_dz_here*TKE_decay(K)) - &
                  (aQ(k) * (TKE_prev(K) - TKE_prev(K+1)) - aQ(k-1) * (TKE_prev(K-1) - TKE_prev(K)))
        Q_err_lin = tke_src + (aQ(k-1) * (dQ(K-1)-dQ(K)) - aQ(k) * (dQ(k)-dQ(k+1))) - &
                    0.5*(TKE_prev(K)-TKE_prev(K+1))*Idz(k)  * (dK(K) + dK(K+1)) - &
                    0.5*(TKE_prev(K)-TKE_prev(K-1))*Idz(k-1)* (dK(K-1) + dK(K)) + &
                    dz_Int(K) * (dK(K) * (S2(K) - N2(K)) - dQ(K)*TKE_decay(K))
      enddo ; endif

    endif  ! End of the Newton's method solver.

    ! Test kappa for convergence...
    if ((tol_err < Newton_err) .and. (.not.abort_Newton)) then
      !  A lower tolerance is used to switch to Newton's method than to switch back.
      Newton_test = Newton_err ; if (do_Newton) Newton_test = 2.0*Newton_err
      was_Newton = do_Newton
      within_tolerance = .true. ; do_Newton = .true.
      do K=min(ks_kappa,ks_kappa_prev),max(ke_kappa,ke_kappa_prev)
        kappa_mean = kappa0 + (kappa(K) - 0.5*dK(K))
        if (abs(dK(K)) > Newton_test * kappa_mean) then
          if (do_Newton) abort_Newton = .true.
          within_tolerance = .false. ; do_Newton = .false. ; exit
        elseif (abs(dK(K)) > tol_err * kappa_mean) then
          within_tolerance = .false. ; if (.not.do_Newton) exit
        endif
        if (abs(dQ(K)) > Newton_test*(tke(K) - 0.5*dQ(K))) then
          if (do_Newton) abort_Newton = .true.
          do_Newton = .false. ; if (.not.within_tolerance) exit
        endif
      enddo

    else  ! Newton's method will not be used again, so no need to check.
      within_tolerance = .true.
      do K=min(ks_kappa,ks_kappa_prev),max(ke_kappa,ke_kappa_prev)
        if (abs(dK(K)) > tol_err * (kappa0 + (kappa(K) - 0.5*dK(K)))) then
          within_tolerance = .false. ; exit
        endif
      enddo
    endif

    if (abort_Newton) then
      do_Newton = .false. ; abort_Newton = .false.
      ! We went to Newton too quickly last time, so restrict the tolerance.
      Newton_err = 0.5*Newton_err
      ke_kappa_prev = nz
      do K=2,nz ; K_Q(K) = kappa(K) / max(TKE(K), TKE_min) ; enddo
    endif

    if (within_tolerance) exit

  enddo

  if (do_Newton) then  ! K_Q needs to be calculated.
    do K=1,ks_kappa-1 ;  K_Q(K) = 0.0 ; enddo
    do K=ks_kappa,ke_kappa ; K_Q(K) = kappa(K) / TKE(K) ; enddo
    do K=ke_kappa+1,nz+1 ; K_Q(K) = 0.0 ; enddo
  endif

  if (present(local_src)) then
    local_src(1) = 0.0 ; local_src(nz+1) = 0.0
    do K=2,nz
      diffusive_src = Idz(k-1)*(kappa(K-1)-kappa(K)) + Idz(k)*(kappa(K+1)-kappa(K))
      chg_by_k0 = kappa0 * ((Idz(k-1)+Idz(k)) / h_Int(K) + I_Ld2(K))
      if (diffusive_src <= 0.0) then
        local_src(K) = K_src(K) + chg_by_k0
      else
        local_src(K) = (K_src(K) + chg_by_k0) + diffusive_src / h_Int(K)
      endif
    enddo
  endif
  if (present(kappa_src)) then
    kappa_src(1) = 0.0 ; kappa_src(nz+1) = 0.0
    do K=2,nz
      kappa_src(K) = K_src(K)
    enddo
  endif

end subroutine find_kappa_tke

!> This subroutine initializes the parameters that regulate shear-driven mixing
function kappa_shear_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(Kappa_shear_CS),    pointer       :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  logical :: kappa_shear_init !< True if module is to be used, False otherwise

  ! Local variables
  real :: KD_normal ! The KD of the main model, read here only as a parameter
                    ! for setting the default of KD_SMOOTH [Z2 T-1 ~> m2 s-1]
  real :: kappa_0_default ! The default value for KD_KAPPA_SHEAR_0 [Z2 T-1 ~> m2 s-1]
  logical :: merge_mixedlayer
  integer :: number_of_OBC_segments
  logical :: debug_shear
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: just_read ! If true, this module is not used, so only read the parameters.
#ifdef __NVCOMPILER_OPENMP_GPU
  !   On the device the whole domain is taken as a single block, so that every column is an
  ! independent task in one kernel launch.
  integer, parameter :: default_niblock = 0 !< Default i block size, 0 being the full domain [nondim].
  integer, parameter :: default_njblock = 0 !< Default j block size, 0 being the full domain [nondim].
#else
  !   A block size of 1 in both horizontal directions reduces the per-column working arrays to a
  ! single column, which is how this scheme is posed on the CPU.
  integer, parameter :: default_niblock = 0 !< Default i block size [nondim].
  integer, parameter :: default_njblock = 1 !< Default j block size [nondim].
#endif
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_kappa_shear"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "kappa_shear_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  !   The Jackson-Hallberg-Legg shear mixing parameterization uses the following
  ! 6 nondimensional coefficients.  That paper gives 3 best fit parameter sets.
  !    Ri_Crit  Rate    FRi_Curv  K_buoy  TKE_N  TKE_Shear
  ! p1: 0.25    0.089    -0.97     0.82    0.24    0.14
  ! p2: 0.30    0.085    -0.94     0.86    0.26    0.13
  ! p3: 0.35    0.088    -0.89     0.81    0.28    0.12
  !   Future research will reveal how these should be modified to take
  ! subgridscale inhomogeneity into account.

! Set default, read and log parameters
  call get_param(param_file, mdl, "USE_JACKSON_PARAM", kappa_shear_init, default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, &
    "Parameterization of shear-driven turbulence following Jackson, Hallberg and Legg, JPO 2008", &
    log_to_all=.true., debugging=kappa_shear_init, all_default=.not.kappa_shear_init)
  call get_param(param_file, mdl, "USE_JACKSON_PARAM", kappa_shear_init, &
                 "If true, use the Jackson-Hallberg-Legg (JPO 2008) "//&
                 "shear mixing parameterization.", default=.false.)
  just_read = .not.kappa_shear_init
  call get_param(param_file, mdl, "VERTEX_SHEAR", CS%KS_at_vertex, &
                 "If true, do the calculations of the shear-driven mixing "//&
                 "at the cell vertices (i.e., the vorticity points).", &
                 default=.false., do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_SHEAR_NIBLOCK", CS%niblock, &
                 "The i-direction block size used to hold the per-column working arrays in the "//&
                 "shear-driven mixing calculations.  The default 0 setting dynamically uses the "//&
                 "full computational domain width.", &
                 default=default_niblock, layoutParam=.true.)
  call get_param(param_file, mdl, "KAPPA_SHEAR_NJBLOCK", CS%njblock, &
                 "The j-direction block size used to hold the per-column working arrays in the "//&
                 "shear-driven mixing calculations.  The default 0 setting dynamically uses the "//&
                 "full computational domain width.", &
                 default=default_njblock, layoutParam=.true.)
  if (CS%niblock < 0) call MOM_error(FATAL, "KAPPA_SHEAR_NIBLOCK must be nonnegative; "//&
                                            "use 0 to select the full computational domain.")
  if (CS%njblock < 0) call MOM_error(FATAL, "KAPPA_SHEAR_NJBLOCK must be nonnegative; "//&
                                            "use 0 to select the full computational domain.")
  !   The vertex scheme interpolates thicknesses to the u- and v-points and then reads the j+1 and
  ! i+1 neighbours of those arrays, so consecutive blocks have to overlap by one row and one
  ! column.  That leaves niblock-1 columns and njblock-1 rows actually written per block, which is
  ! impossible with a block size of 1.  This is the one case where the requested block size cannot
  ! be honoured, so say so rather than failing or silently reading past the end of a block.  The
  ! tracer-point scheme has no such stencil and is left alone.
  if (CS%KS_at_vertex .and. .not.just_read) then
    if (CS%niblock == 1) then
      call MOM_error(WARNING, "KAPPA_SHEAR_NIBLOCK must be >= 2, or 0, when VERTEX_SHEAR is "//&
                     "true, because consecutive blocks overlap by one column.  Changing the i "//&
                     "block size from 1 to 2 for this run.")
      CS%niblock = 2
    endif
    if (CS%njblock == 1) then
      call MOM_error(WARNING, "KAPPA_SHEAR_NJBLOCK must be >= 2, or 0, when VERTEX_SHEAR is "//&
                     "true, because consecutive blocks overlap by one row.  Changing the j "//&
                     "block size from 1 to 2 for this run.")
      CS%njblock = 2
    endif
  endif
  call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
  call get_param(param_file, mdl, "VERTEX_SHEAR_VISCOSITY_BUG", CS%VS_viscosity_bug, &
                 "If true, use a bug in vertex shear that zeros out viscosities at "//&
                 "vertices on coastlines.", &
                 default=enable_bugs, do_not_log=just_read.or.(.not.CS%KS_at_vertex))
  call get_param(param_file, mdl, "OBC_NUMBER_OF_SEGMENTS", number_of_OBC_segments, &
                 default=0, do_not_log=.true.)
  call get_param(param_file, mdl, "VERTEX_SHEAR_OBC_BUG", CS%vertex_shear_OBC_bug, &
                 "If false, use extra masking when interpolating thicknesses to velocity "//&
                 "points for setting up the shear velocities at vertices to avoid using "//&
                 "external thicknesses at open boundaries.  When OBCs are not in use, "//&
                 "this parameter does not change answers, but true is more efficient.", &
                 default=enable_bugs, &
                 do_not_log=just_read.or.(.not.CS%KS_at_vertex).or.(number_of_OBC_segments<=0))
                 ! Use OBC settings to set the default for VERTEX_SHEAR_OBC_BUG?
  call get_param(param_file, mdl, "VERTEX_SHEAR_GEOMETRIC_MEAN", CS%VS_GeometricMean, &
                 "If true, use a geometric mean for moving diffusivity from "//&
                 "vertices to tracer points.  False uses algebraic mean.", &
                 default=.false., do_not_log=just_read.or.(.not.CS%KS_at_vertex))
  call get_param(param_file, mdl, "VERTEX_SHEAR_THICKNESS_MEAN", CS%VS_ThicknessMean, &
                 "If true, apply thickness weighting to horizontal averagings of diffusivity "//&
                 "to tracer points in the kappa shear solver.", &
                 default=.false.)
  if (CS%VS_GeometricMean) then
    call get_param(param_file, mdl, "VERTEX_SHEAR_GEOMETRIC_MEAN_KDMIN", &
                   CS%VS_GeoMean_Kdmin, "If using the geometric mean in vertex shear, "//&
                   "use this minimum value for Kd. This is an ad-hoc parameter, the "//&
                   "diffusivities on the edge of shear regions are sensitive to the choice.",&
                   units="m2 s-1", default=0.0, scale=GV%m2_s_to_HZ_T, do_not_log=just_read)
  endif
  call get_param(param_file, mdl, "RINO_CRIT", CS%RiNo_crit, &
                 "The critical Richardson number for shear mixing.", &
                 units="nondim", default=0.25, do_not_log=just_read)
  call get_param(param_file, mdl, "SHEARMIX_RATE", CS%Shearmix_rate, &
                 "A nondimensional rate scale for shear-driven entrainment. "//&
                 "Jackson et al find values in the range of 0.085-0.089.", &
                 units="nondim", default=0.089, do_not_log=just_read)
  call get_param(param_file, mdl, "MAX_RINO_IT", CS%max_RiNo_it, &
                 "The maximum number of iterations that may be used to "//&
                 "estimate the Richardson number driven mixing.", &
                 units="nondim", default=50, do_not_log=just_read)
  call get_param(param_file, mdl, "KD", KD_normal, &
                 units="m2 s-1", scale=US%m2_s_to_Z2_T, default=0.0, do_not_log=.true.)
  kappa_0_default = max(Kd_normal, 1.0e-7*US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "KD_KAPPA_SHEAR_0", CS%kappa_0, &
                 "The background diffusivity that is used to smooth the "//&
                 "density and shear profiles before solving for the "//&
                 "diffusivities.  The default is the greater of KD and 1e-7 m2 s-1.", &
                 units="m2 s-1", default=kappa_0_default*US%Z2_T_to_m2_s, scale=GV%m2_s_to_HZ_T, &
                 do_not_log=just_read)
  call get_param(param_file, mdl, "KD_SEED_KAPPA_SHEAR", CS%kappa_seed, &
                 "A moderately large seed value of diapycnal diffusivity that is used as a "//&
                 "starting turbulent diffusivity in the iterations to find an energetically "//&
                 "constrained solution for the shear-driven diffusivity.", &
                 units="m2 s-1", default=1.0, scale=GV%m2_s_to_HZ_T)
  call get_param(param_file, mdl, "KD_TRUNC_KAPPA_SHEAR", CS%kappa_trunc, &
                 "The value of shear-driven diffusivity that is considered negligible "//&
                 "and is rounded down to 0. The default is 1% of KD_KAPPA_SHEAR_0.", &
                 units="m2 s-1", default=0.01*CS%kappa_0*GV%HZ_T_to_m2_s, scale=GV%m2_s_to_HZ_T, &
                 do_not_log=just_read)
  call get_param(param_file, mdl, "FRI_CURVATURE", CS%FRi_curvature, &
                 "The nondimensional curvature of the function of the "//&
                 "Richardson number in the kappa source term in the "//&
                 "Jackson et al. scheme.", units="nondim", default=-0.97, do_not_log=just_read)
  call get_param(param_file, mdl, "TKE_N_DECAY_CONST", CS%C_N, &
                 "The coefficient for the decay of TKE due to "//&
                 "stratification (i.e. proportional to N*tke). "//&
                 "The values found by Jackson et al. are 0.24-0.28.", &
                 units="nondim", default=0.24, do_not_log=just_read)
!  call get_param(param_file, mdl, "LAYER_KAPPA_STAGGER", CS%layer_stagger, &
!                 default=.false., do_not_log=just_read)
  call get_param(param_file, mdl, "TKE_SHEAR_DECAY_CONST", CS%C_S, &
                 "The coefficient for the decay of TKE due to shear (i.e. "//&
                 "proportional to |S|*tke). The values found by Jackson "//&
                 "et al. are 0.14-0.12.", units="nondim", default=0.14, do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_BUOY_SCALE_COEF", CS%lambda, &
                 "The coefficient for the buoyancy length scale in the "//&
                 "kappa equation.  The values found by Jackson et al. are "//&
                 "in the range of 0.81-0.86.", units="nondim", default=0.82, do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_N_OVER_S_SCALE_COEF2", CS%lambda2_N_S, &
                 "The square of the ratio of the coefficients of the "//&
                 "buoyancy and shear scales in the diffusivity equation, "//&
                 "Set this to 0 (the default) to eliminate the shear scale. "//&
                 "This is only used if USE_JACKSON_PARAM is true.", &
                 units="nondim", default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl, "LZ_RESCALE", CS%lz_rescale, &
                 "A coefficient to rescale the distance to the nearest solid boundary. "//&
                 "This adjustment is to account for regions where 3 dimensional turbulence "//&
                 "prevents the growth of shear instabilities [nondim].", &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "KAPPA_SHEAR_TOL_ERR", CS%kappa_tol_err, &
                 "The fractional error in kappa that is tolerated. "//&
                 "Iteration stops when changes between subsequent "//&
                 "iterations are smaller than this everywhere in a "//&
                 "column.  The peak diffusivities usually converge most "//&
                 "rapidly, and have much smaller errors than this.", &
                 units="nondim", default=0.1, do_not_log=just_read)
  call get_param(param_file, mdl, "TKE_BACKGROUND", CS%TKE_bg, &
                 "A background level of TKE used in the first iteration "//&
                 "of the kappa equation.  TKE_BACKGROUND could be 0.", &
                 units="m2 s-2", default=0.0, scale=US%m_to_Z**2*US%T_to_s**2)
  call get_param(param_file, mdl, "KAPPA_SHEAR_ELIM_MASSLESS", CS%eliminate_massless, &
                 "If true, massless layers are merged with neighboring "//&
                 "massive layers in this calculation.  The default is "//&
                 "true and I can think of no good reason why it should "//&
                 "be false. This is only used if USE_JACKSON_PARAM is true.", &
                 default=.true., do_not_log=just_read)
  call get_param(param_file, mdl, "MAX_KAPPA_SHEAR_IT", CS%max_KS_it, &
                 "The maximum number of iterations that may be used to "//&
                 "estimate the time-averaged diffusivity.", &
                 default=13, do_not_log=just_read)
  call get_param(param_file, mdl, "PRANDTL_TURB", CS%Prandtl_turb, &
                 "The turbulent Prandtl number applied to shear instability.", &
                 units="nondim", default=1.0, do_not_log=just_read)
  call get_param(param_file, mdl, "VEL_UNDERFLOW", CS%vel_underflow, &
                 "A negligibly small velocity magnitude below which velocity components are set "//&
                 "to 0.  A reasonable value might be 1e-30 m/s, which is less than an "//&
                 "Angstrom divided by the age of the universe.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T, do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_SHEAR_MAX_KAP_SRC_CHG", CS%kappa_src_max_chg, &
                 "The maximum permitted increase in the kappa source within an iteration relative "//&
                 "to the local source; this must be greater than 1.  The lower limit for the "//&
                 "permitted fractional decrease is (1 - 0.5/kappa_src_max_chg).  These limits "//&
                 "could perhaps be made dynamic with an improved iterative solver.", &
                 default=10.0, units="nondim", do_not_log=just_read)

  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true., do_not_log=just_read)
  call get_param(param_file, mdl, "DEBUG_KAPPA_SHEAR", debug_shear, &
                 "If true, write debugging data for the kappa-shear code.", &
                 default=.false., debuggingParam=.true., do_not_log=.true.)
  if (debug_shear) CS%debug = .true.
  call get_param(param_file, mdl, "KAPPA_SHEAR_VERTEX_PSURF_BUG", CS%psurf_bug, &
                 "If true, do a simple average of the cell surface pressures to get a pressure "//&
                 "at the corner if VERTEX_SHEAR=True.  Otherwise mask out any land points in "//&
                 "the average.", default=.false., do_not_log=(just_read .or. (.not.CS%KS_at_vertex)))

  call get_param(param_file, mdl, "KAPPA_SHEAR_ITER_BUG", CS%dKdQ_iteration_bug, &
                 "If true, use an older, dimensionally inconsistent estimate of the "//&
                 "derivative of diffusivity with energy in the Newton's method iteration.  "//&
                 "The bug causes undercorrections when dz > 1 m.", default=.false., do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_SHEAR_ALL_LAYER_TKE_BUG", CS%all_layer_TKE_bug, &
                 "If true, report back the latest estimate of TKE instead of the time average "//&
                 "TKE when there is mass in all layers.  Otherwise always report the time "//&
                 "averaged TKE, as is currently done when there are some massless layers.", &
                 default=.false., do_not_log=just_read)
  call get_param(param_file, mdl, "USE_RESTRICTIVE_TOLERANCE_CHECK", CS%restrictive_tolerance_check, &
                 "If true, uses the more restrictive tolerance check to determine if a timestep "//&
                 "is acceptable for the KS_it outer iteration loop.  False uses the original less "//&
                 "restrictive check.", default=.false., do_not_log=just_read)
!    id_clock_KQ = cpu_clock_id('Ocean KS kappa_shear', grain=CLOCK_ROUTINE)
!    id_clock_avg = cpu_clock_id('Ocean KS avg', grain=CLOCK_ROUTINE)
!    id_clock_project = cpu_clock_id('Ocean KS project', grain=CLOCK_ROUTINE)
!    id_clock_setup = cpu_clock_id('Ocean KS setup', grain=CLOCK_ROUTINE)

  CS%nkml = 1
  if (GV%nkml>0) then
    call get_param(param_file, mdl, "KAPPA_SHEAR_MERGE_ML",merge_mixedlayer, &
                 "If true, combine the mixed layers together before solving the "//&
                 "kappa-shear equations.", default=.true., do_not_log=just_read)
    if (merge_mixedlayer) CS%nkml = GV%nkml
  endif

! Forego remainder of initialization if not using this scheme
  if (.not. kappa_shear_init) return

  CS%diag => diag

  CS%id_Kd_shear = register_diag_field('ocean_model','Kd_shear', diag%axesTi, Time, &
      'Shear-driven Diapycnal Diffusivity at horizontal tracer points', 'm2 s-1', conversion=GV%HZ_T_to_m2_s)
  if (CS%KS_at_vertex) then
    CS%id_TKE = register_diag_field('ocean_model','TKE_shear', diag%axesBi, Time, &
       'Shear-driven Turbulent Kinetic Energy at horizontal vertices', 'm2 s-2', conversion=US%Z_to_m**2*US%s_to_T**2)
    CS%id_Kd_vertex = register_diag_field('ocean_model','Kd_shear_vertex', diag%axesBi, Time, &
         'Shear-driven Diapycnal Diffusivity at horizontal vertices', 'm2 s-1', conversion=GV%HZ_T_to_m2_s)
    CS%id_S2_init = register_diag_field('ocean_model','S2_shear_in', diag%axesBi, Time, &
         'Interface shear squared at horizontal vertices, as input to kappa-shear', 's-2', conversion=US%s_to_T**2)
    CS%id_N2_init = register_diag_field('ocean_model','N2_shear_in', diag%axesBi, Time, &
         'Interface stratification at horizontal vertices, as input to kappa-shear', 's-2', conversion=US%s_to_T**2)
    CS%id_S2_mean = register_diag_field('ocean_model','S2_shear_mean', diag%axesBi, Time, &
         'Interface shear squared at horizontal vertices, averaged over timestep in kappa-shear', &
         's-2', conversion=US%s_to_T**2)
    CS%id_N2_mean = register_diag_field('ocean_model','N2_shear_mean', diag%axesBi, Time, &
         'Interface stratification at horizontal vertices, averaged over timestep in kappa-shear', &
         's-2', conversion=US%s_to_T**2)
  else
    CS%id_TKE = register_diag_field('ocean_model','TKE_shear', diag%axesTi, Time, &
         'Shear-driven Turbulent Kinetic Energy at horizontal tracer points', &
         'm2 s-2', conversion=US%Z_to_m**2*US%s_to_T**2)
    CS%id_S2_init = register_diag_field('ocean_model','S2_shear_in', diag%axesTi, Time, &
         'Interface shear squared at horizontal tracer points, as input to kappa-shear', 's-2', conversion=US%s_to_T**2)
    CS%id_N2_init = register_diag_field('ocean_model','N2_shear_in', diag%axesTi, Time, &
         'Interface stratification at horizontal tracer points, as input to kappa-shear', &
         's-2', conversion=US%s_to_T**2)
    CS%id_S2_mean = register_diag_field('ocean_model','S2_shear_mean', diag%axesTi, Time, &
         'Interface shear squared at horizontal tracer points, averaged over timestep in kappa-shear', &
         's-2', conversion=US%s_to_T**2)
    CS%id_N2_mean = register_diag_field('ocean_model','N2_shear_mean', diag%axesTi, Time, &
         'Interface stratification at horizontal tracer points, averaged ove timestep in kappa-shear', &
         's-2', conversion=US%s_to_T**2)
  endif

  ! !$omp target enter data map(to: CS)

#ifdef __NVCOMPILER_OPENMP_GPU
  if (kappa_shear_init .and. (GV%ke > GPU_nk_max)) call MOM_error(FATAL, &
    "kappa_shear_init: GPU builds of kappa_shear require GV%ke <= GPU_nk_max because the "//&
    "column routines use fixed-size local arrays on the device (this applies to the "//&
    "vertex scheme too, which shares those routines); increase GPU_nk_max in "//&
    "MOM_kappa_shear.F90 or use a CPU build.")
#endif

end function kappa_shear_init

!> This function indicates to other modules whether the Jackson et al shear mixing
!! parameterization will be used without needing to duplicate the log entry.
logical function kappa_shear_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters

  ! Local variables
  character(len=40)  :: mdl = "MOM_kappa_shear"  ! This module's name.
  ! This function reads the parameter "USE_JACKSON_PARAM" and returns its value.

  call get_param(param_file, mdl, "USE_JACKSON_PARAM", kappa_shear_is_used, &
                 default=.false., do_not_log=.true.)
end function kappa_shear_is_used

!> This function indicates to other modules whether the Jackson et al shear mixing parameterization
!! will be used at the vertices without needing to duplicate the log entry.  It returns false if
!! the Jackson et al scheme is not used or if it is used via calculations at the tracer points.
logical function kappa_shear_at_vertex(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters

  ! Local variables
  character(len=40)  :: mdl = "MOM_kappa_shear"  ! This module's name.
  logical :: do_kappa_shear
  ! This function returns true only if the parameters "USE_JACKSON_PARAM" and "VERTEX_SHEAR" are both true.

  kappa_shear_at_vertex = .false.

  call get_param(param_file, mdl, "USE_JACKSON_PARAM", do_kappa_shear, &
                 default=.false., do_not_log=.true.)
  if (do_Kappa_Shear) &
    call get_param(param_file, mdl, "VERTEX_SHEAR", kappa_shear_at_vertex, &
                 "If true, do the calculations of the shear-driven mixing "//&
                 "at the cell vertices (i.e., the vorticity points).", &
                 default=.false., do_not_log=.true.)

end function kappa_shear_at_vertex

!> \namespace mom_kappa_shear
!!
!! By Laura Jackson and Robert Hallberg, 2006-2008
!!
!!   This file contains the subroutines that determine the diapycnal
!! diffusivity driven by resolved shears, as specified by the
!! parameterizations described in Jackson and Hallberg (JPO, 2008).
!!
!!   The technique by which the 6 equations (for kappa, TKE, u, v, T,
!! and S) are solved simultaneously has been dramatically revised
!! from the previous version. The previous version was not converging
!! in some cases, especially near the surface mixed layer, while the
!! revised version does.  The revised version solves for kappa and
!! TKE with shear and stratification fixed, then marches the density
!! and velocities forward with an adaptive (and aggressive) time step
!! in a predictor-corrector-corrector emulation of a trapezoidal
!! scheme.  Run-time-settable parameters determine the tolerance to
!! which the kappa and TKE equations are solved and the minimum time
!! step that can be taken.

end module MOM_kappa_shear
