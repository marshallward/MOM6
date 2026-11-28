! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

#include "do_concurrent_compat.h"

!> Implements vertical viscosity (vertvisc)
module MOM_vert_friction

use MOM_domains,       only : pass_var, To_All, Omit_corners
use MOM_domains,       only : pass_vector, Scalar_Pair
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : post_product_u, post_product_sum_u
use MOM_diag_mediator, only : post_product_v, post_product_sum_v
use MOM_diag_mediator, only : diag_ctrl, query_averaging_enabled
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,       only : To_North, To_East
use MOM_debugging,     only : uvchksum, hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : mech_forcing, find_ustar
use MOM_get_input,     only : directories
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : MOM_read_data, slasher
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary, only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_PointAccel,    only : write_u_accel, write_v_accel, PointAccel_init
use MOM_PointAccel,    only : PointAccel_CS
use MOM_time_manager,  only : time_type, time_minus_signed
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, vertvisc_type
use MOM_variables,     only : cont_diag_ptrs, accel_diag_ptrs
use MOM_variables,     only : ocean_internal_state
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS
use MOM_set_visc,      only : set_v_at_u, set_u_at_v
use MOM_lateral_mixing_coeffs, only : VarMix_CS

use CVMix_kpp,         only : cvmix_kpp_composite_Gshape

implicit none ; private

#include <MOM_memory.h>

public vertvisc, vertvisc_remnant, vertvisc_coef
public vertvisc_limit_vel, vertvisc_init, vertvisc_end
public updateCFLtruncationValue
public vertFPmix

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure with parameters and memory for the MOM_vert_friction module
type, public :: vertvisc_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: Hmix            !< The mixed layer thickness [Z ~> m].
  real    :: Hmix_stress     !< The mixed layer thickness over which the wind
                             !! stress is applied with direct_stress [H ~> m or kg m-2].
  real    :: Kvml_invZ2      !< The extra vertical viscosity scale in [H Z T-1 ~> m2 s-1 or Pa s] in a
                             !! surface mixed layer with a characteristic thickness given by Hmix,
                             !! and scaling proportional to (Hmix/z)^2, where z is the distance
                             !! from the surface; this can get very large with thin layers.
  real    :: Kv              !< The interior vertical viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real    :: Hbbl            !< The static bottom boundary layer thickness [Z ~> m].
  real    :: Hbbl_gl90       !< The static bottom boundary layer thickness used for GL90 [Z ~> m].
  real    :: Kv_extra_bbl    !< An extra vertical viscosity in the bottom boundary layer of thickness
                             !! Hbbl when there is not a bottom drag law in use [H Z T-1 ~> m2 s-1 or Pa s].
  real    :: vonKar          !< The von Karman constant as used for mixed layer viscosity [nondim]

  logical :: use_GL90_in_SSW !< If true, use the GL90 parameterization in stacked shallow water mode (SSW).
                             !! The calculation of the GL90 viscosity coefficient uses the fact that in SSW
                             !! we simply have 1/N^2 = h/g^prime, where g^prime is the reduced gravity.
                             !! This identity does not generalize to non-SSW setups.
  logical :: use_GL90_N2     !< If true, use GL90 vertical viscosity coefficient that is depth-independent;
                             !! this corresponds to a kappa_GM that scales as N^2 with depth.
  real    :: kappa_gl90      !< The scalar diffusivity used in the GL90 vertical viscosity scheme
                             !! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
  logical :: read_kappa_gl90 !< If true, read a file containing the spatially varying kappa_gl90
  real    :: alpha_gl90      !< Coefficient used to compute a depth-independent GL90 vertical
                             !! viscosity via Kv_gl90 = alpha_gl90 * f^2. Note that the implied
                             !! Kv_gl90 corresponds to a kappa_gl90 that scales as N^2 with depth.
                             !! [H Z T ~> m2 s or kg s m-1]
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  real    :: CFL_trunc       !< Velocity components will be truncated when they
                             !! are large enough that the corresponding CFL number
                             !! exceeds this value [nondim].
  real    :: CFL_report      !< The value of the CFL number that will cause the
                             !! accelerations to be reported [nondim].  CFL_report
                             !! will often equal CFL_trunc.
  real    :: truncRampTime   !< The time-scale over which to ramp up the value of
                             !! CFL_trunc from CFL_truncS to CFL_truncE [T ~> s]
  real    :: CFL_truncS      !< The start value of CFL_trunc [nondim]
  real    :: CFL_truncE      !< The end/target value of CFL_trunc [nondim]
  logical :: CFLrampingIsActivated = .false. !< True if the ramping has been initialized
  type(time_type) :: rampStartTime !< The time at which the ramping of CFL_trunc starts

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NK_INTERFACE_) :: &
    a_u                !< The u-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NK_INTERFACE_) :: &
    a_u_gl90           !< The u-drag coefficient associated with GL90 across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    h_u                !< The effective layer thickness at u-points [H ~> m or kg m-2].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NK_INTERFACE_) :: &
    a_v                !< The v-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NK_INTERFACE_) :: &
    a_v_gl90           !< The v-drag coefficient associated with GL90 across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    h_v                !< The effective layer thickness at v-points [H ~> m or kg m-2].
  real, pointer, dimension(:,:) :: a1_shelf_u => NULL() !< The u-momentum coupling coefficient under
                           !! ice shelves [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress under shelves.
  real, pointer, dimension(:,:) :: a1_shelf_v => NULL() !< The v-momentum coupling coefficient under
                           !! ice shelves [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress under shelves.

  logical :: split          !< If true, use the split time stepping scheme.
  logical :: bottomdraglaw  !< If true, the  bottom stress is calculated with a
                            !! drag law c_drag*|u|*u. The velocity magnitude
                            !! may be an assumed value or it may be based on the
                            !! actual velocity in the bottommost HBBL, depending
                            !! on whether linear_drag is true.
  logical :: harmonic_visc  !< If true, the harmonic mean thicknesses are used
                            !! to calculate the viscous coupling between layers
                            !! except near the bottom.  Otherwise the arithmetic
                            !! mean thickness is used except near the bottom.
  real    :: harm_BL_val    !< A scale to determine when water is in the boundary
                            !! layers based solely on harmonic mean thicknesses
                            !! for the purpose of determining the extent to which
                            !! the thicknesses used in the viscosities are upwinded [nondim].
  logical :: direct_stress  !< If true, the wind stress is distributed over the topmost Hmix_stress
                            !! of fluid, and an added mixed layer viscosity or a physically based
                            !! boundary layer turbulence parameterization is not needed for stability.
  logical :: dynamic_viscous_ML  !< If true, use the results from a dynamic
                            !! calculation, perhaps based on a bulk Richardson
                            !! number criterion, to determine the mixed layer
                            !! thickness for viscosity.
  logical :: fixed_LOTW_ML  !< If true, use a Law-of-the-wall prescription for the mixed layer
                            !! viscosity within a boundary layer that is the lesser of Hmix and the
                            !! total depth of the ocean in a column.
  logical :: apply_LOTW_floor !< If true, use a Law-of-the-wall prescription to set a lower bound
                            !! on the viscous coupling between layers within the surface boundary
                            !! layer, based the distance of interfaces from the surface.  This only
                            !! acts when there are large changes in the thicknesses of successive
                            !! layers or when the viscosity is set externally and the wind stress
                            !! has subsequently increased.
  integer :: answer_date    !< The vintage of the order of arithmetic and expressions in the viscous
                            !! calculations.  Values below 20190101 recover the answers from the end
                            !! of 2018, while higher values use expressions that do not use an
                            !! arbitrary and hard-coded maximum viscous coupling coefficient between
                            !! layers.  In non-Boussinesq cases, values below 20230601 recover a
                            !! form of the viscosity within  the mixed layer that breaks up the
                            !! magnitude of the wind stress with BULKMIXEDLAYER, DYNAMIC_VISCOUS_ML
                            !! or FIXED_DEPTH_LOTW_ML, but not LOTW_VISCOUS_ML_FLOOR.
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  integer :: niblock        !< Horizontal i block size; 0 uses the full computational extent.
  integer :: njblock        !< Horizontal j block size; 0 uses the full computational extent.
  integer :: nkml           !< The number of layers in the mixed layer.
  integer, pointer :: ntrunc !< The number of times the velocity has been
                            !! truncated since the last call to write_energy.
  character(len=200) :: u_trunc_file  !< The complete path to a file in which a column of
                            !! u-accelerations are written if velocity truncations occur.
  character(len=200) :: v_trunc_file !< The complete path to a file in which a column of
                            !! v-accelerations are written if velocity truncations occur.
  logical :: StokesMixing   !< If true, do Stokes drift mixing via the Lagrangian current
                            !! (Eulerian plus Stokes drift).  False by default and set
                            !! via STOKES_MIXING_COMBINED.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  real, allocatable, dimension(:,:) :: kappa_gl90_2d !< 2D kappa_gl90 at h-points [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]

  !>@{ Diagnostic identifiers
  integer :: id_du_dt_visc = -1, id_dv_dt_visc = -1, id_du_dt_visc_gl90 = -1, id_dv_dt_visc_gl90 = -1
  integer :: id_GLwork = -1
  integer :: id_au_vv = -1, id_av_vv = -1, id_au_gl90_vv = -1, id_av_gl90_vv = -1
  integer :: id_du_dt_str = -1, id_dv_dt_str = -1
  integer :: id_h_u = -1, id_h_v = -1, id_hML_u = -1 , id_hML_v = -1
  integer :: id_Omega_w2x = -1, id_FPtau2s  = -1 , id_FPtau2w = -1
  integer :: id_uE_h  = -1, id_vE_h  = -1
  integer :: id_uStk  = -1, id_vStk  = -1
  integer :: id_uStk0 = -1, id_vStk0 = -1
  integer :: id_uInc_h= -1, id_vInc_h= -1
  integer :: id_taux_bot = -1, id_tauy_bot = -1
  integer :: id_Kv_slow = -1, id_Kv_u = -1, id_Kv_v = -1
  integer :: id_Kv_gl90_u = -1, id_Kv_gl90_v = -1
  ! integer :: id_hf_du_dt_visc    = -1, id_hf_dv_dt_visc    = -1
  integer :: id_h_du_dt_visc    = -1, id_h_dv_dt_visc    = -1
  integer :: id_hf_du_dt_visc_2d = -1, id_hf_dv_dt_visc_2d = -1
  integer :: id_h_du_dt_str    = -1, id_h_dv_dt_str    = -1
  integer :: id_du_dt_str_visc_rem = -1, id_dv_dt_str_visc_rem = -1
  !>@}

  type(PointAccel_CS), pointer :: PointAccel_CSp => NULL() !< A pointer to the control structure
                              !! for recording accelerations leading to velocity truncations

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes
end type vertvisc_CS

contains

!> Add nonlocal stress increments to ui^n and vi^n.
subroutine vertFPmix(ui, vi, uold, vold, hbl_h, h, forces, dt, lpost, Cemp_NL, G, GV, US, CS, OBC, Waves)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: ui     !< Zonal velocity after vertvisc [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vi     !< Meridional velocity after vertvisc [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uold   !< Old Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vold   !< Old Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: hbl_h !<  boundary layer depth [H ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h       !< Layer thicknesses [H ~> m or kg m-2]
  type(mech_forcing),      intent(in) :: forces  !< A structure with the driving mechanical forces
  real,                    intent(in) :: dt      !< Time increment [T ~> s]
  real,                    intent(in) :: Cemp_NL !< empirical coefficient of non-local momentum mixing [nondim]
  logical,                 intent(in) :: lpost   !< Compute and make available FPMix diagnostics
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(vertvisc_CS),       intent(inout) :: CS   !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer    :: OBC     !< Open boundary condition structure
  type(wave_parameters_CS), &
                   optional, pointer  :: Waves   !< Container for wave/Stokes information

  ! local variables
  real, dimension(SZIB_(G),SZJ_(G))  :: hbl_u   !< boundary layer depth (u-pts) [H ~> m]
  real, dimension(SZI_(G),SZJB_(G))  :: hbl_v   !< boundary layer depth (v-pts) [H ~> m]
  real, dimension(SZIB_(G),SZJ_(G))  :: taux_u  !< kinematic zonal wind stress (u-pts) [L Z T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJB_(G))  :: tauy_v  !< kinematic merid wind stress (v-pts) [L Z T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G))   :: uS0     !< surface zonal Stokes drift h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G))   :: vS0     !< surface zonal Stokes drift h-pts [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uE_u    !< zonal Eulerian u-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uE_h    !< zonal Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vE_v    !< merid Eulerian v-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vE_h    !< merid Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uInc_u  !< zonal Eulerian u-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uInc_h  !< zonal Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vInc_v  !< merid Eulerian v-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vInc_h  !< merid Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uStk    !< zonal Stokes Drift (h-pts) [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vStk    !< merid Stokes Drift (h-pts) [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)+1) :: omega_tau2s !< angle stress to shear (h-pts) [rad]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)+1) :: omega_tau2w !< angle stress to wind  (h-pts) [rad]
  real :: omega_tmp, omega_s2x, omega_tau2x                    !< temporary angle wrt the x axis [rad]
  real :: Irho0        !< Inverse of the mean density [R-1 ~> m3 kg-1]
  real :: pi           !< ! The ratio of the circumference of a circle to its diameter [nondim]
  real :: tmp_u, tmp_v !< temporary ocean mask weights on u and v points [nondim]
  real :: fexp         !< temporary exponential function [nondim]
  real :: sigma        !< temporary normalize boundary layer coordinate [nondim]
  real :: Gat1, Gsig, dGdsig !< Shape parameters [nondim]
  real :: du, dv       !< Intermediate velocity differences [L T-1 ~> m s-1]
  real :: depth        !< Cumulative of thicknesses [H ~> m]
  integer :: b, kp1, k, nz !< band and vertical indices
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq !< horizontal indices

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke

  pi = 4. * atan2(1.,1.)
  Irho0 = 1.0 / GV%Rho0

  ! initialize arrays
  uE_h(:,:,:)   = 0.0
  vE_h(:,:,:)   = 0.0
  uE_u(:,:,:)   = 0.0
  vE_v(:,:,:)   = 0.0
  vInc_v(:,:,:) = 0.0
  uInc_u(:,:,:) = 0.0
  call pass_var(hbl_h , G%Domain, halo=1)

  ! u-points
  do j = js,je
    do I = Isq,Ieq
      taux_u(I,j)  = forces%taux(I,j) * Irho0
      if ( (G%mask2dCu(I,j) > 0.5) ) then
        ! h to u-pts
        tmp_u  = MAX (1.0 ,(G%mask2dT(i,j) + G%mask2dT(i+1,j) ) )
        hbl_u(I,j) = ((G%mask2dT(i,j) * hbl_h(i,j)) + (G%mask2dT(i+1,j) * hbl_h(i+1,j))) / tmp_u
        depth   = 0.
        Gat1  = 0.
        do k=1, nz
          ! cell center
          depth = depth + 0.5*CS%h_u(I,j,k)
          uE_u(I,j,k) = ui(I,j,k) - waves%Us_x(I,j,k)
          if ( depth < hbl_u(I,j) )     then
            sigma = depth / hbl_u(i,j)
            ! cell bottom
            depth = depth + 0.5*CS%h_u(I,j,k)
            call cvmix_kpp_composite_Gshape(sigma,Gat1,Gsig,dGdsig)
            ! nonlocal boundary-layer increment
            uInc_u(I,j,k)  = dt * Cemp_NL * taux_u(I,j) * dGdsig / (hbl_u(I,j) * GV%H_to_Z)
            ui(I,j,k) = ui(I,j,k) + uInc_u(I,j,k)
          else
            uInc_u(I,j,k) = 0.0
          endif
        enddo
      else
        do k=1, nz
          uInc_u(I,j,k) = 0.0
        enddo
      endif
    enddo
  enddo

  ! v-points
  do J = Jsq,Jeq
    do i = is,ie
      tauy_v(i,J)  = forces%tauy(i,J) * Irho0
      if ( (G%mask2dCv(i,J) > 0.5) ) then
        ! h to v-pts
        tmp_v  = max( 1.0 ,(G%mask2dT(i,j) + G%mask2dT(i,j+1)))
        hbl_v(i,J) = (G%mask2dT(i,j) * hbl_h(i,J) + G%mask2dT(i,j+1) * hbl_h(i,j+1)) / tmp_v
        depth = 0.
        Gat1  = 0.
        do k=1, nz
          ! cell center
          depth = depth + 0.5* CS%h_v(i,J,k)
          vE_v(i,J,k) = vi(i,J,k) - waves%Us_y(i,J,k)
          if ( depth < hbl_v(i,J) )    then
            sigma = depth / hbl_v(i,J)
            ! cell bottom
            depth = depth + 0.5* CS%h_v(i,J,k)
            call cvmix_kpp_composite_Gshape(sigma,Gat1,Gsig,dGdsig)
            ! nonlocal boundary-layer increment
            vInc_v(i,J,k) = dt * Cemp_NL * tauy_v(i,J) * dGdsig / (hbl_v(i,J) * GV%H_to_Z)
            vi(i,J,k) = vi(i,J,k) + vInc_v(i,J,k)
          else
            vInc_v(i,J,k)  = 0.0
          endif
        enddo
      else
        do k=1, nz
          vInc_v(i,J,k)  = 0.0
        enddo
      endif
    enddo
  enddo

  ! Compute and store diagnostics, only during the corrector step.
  if (lpost)  then
    call pass_vector(uE_u  ,  vE_v  , G%Domain, To_All)
    call pass_vector(uInc_u, vInc_v , G%Domain, To_All)
    uStk = 0.0
    vStk = 0.0
    uS0  = 0.0
    vS0  = 0.0

    do j = js,je
      do i = is,ie
        if (G%mask2dT(i,j) > 0.5)  then
          ! u to h-pts
          tmp_u  = max( 1.0 ,(G%mask2dCu(i,j) + G%mask2dCu(i-1,j)))
          ! v to h-pts
          tmp_v  = max( 1.0 ,(G%mask2dCv(i,j) + G%mask2dCv(i,j-1)))
          do k = 1,nz
            uE_h(i,j,k)   = (G%mask2dCu(i,j) *   uE_u(i,j,k) + G%mask2dCu(i-1,j) *   uE_u(i-1,j,k)) / tmp_u
            uInc_h(i,j,k) = (G%mask2dCu(i,j) * uInc_u(i,j,k) + G%mask2dCu(i-1,j) * uInc_u(i-1,j,k)) / tmp_u
            vE_h(i,j,k)   = (G%mask2dCv(i,j) *   vE_v(i,j,k) + G%mask2dCv(i,j-1) *   vE_v(i,j-1,k)) / tmp_v
            vInc_h(i,j,k) = (G%mask2dCv(i,j) * vInc_v(i,j,k) + G%mask2dCv(i,j-1) * vInc_v(i,j-1,k)) / tmp_v
          enddo
          ! Wind, Stress and Shear align at surface
          Omega_tau2w(i,j,:) = 0.0
          Omega_tau2s(i,j,:) = 0.0
          do k = 1,nz
            kp1 = min( nz , k+1)
            du = uE_h(i,j,k) - uE_h(i,j,kp1)
            dv = vE_h(i,j,k) - vE_h(i,j,kp1)
            omega_s2x = atan2(dv, du)

            du = du + uInc_h(i,j,k) - uInc_h(i,j,kp1)
            dv = dv + vInc_h(i,j,k) - vInc_h(i,j,kp1)
            omega_tau2x = atan2(dv, du)
            omega_tmp = omega_tau2x - forces%omega_w2x(i,j)

            if ( (omega_tmp  >   pi   ) ) omega_tmp = omega_tmp - 2.*pi
            if ( (omega_tmp  < (0.-pi)) )  omega_tmp = omega_tmp + 2.*pi
            Omega_tau2w(i,j,kp1) = omega_tmp

            omega_tmp = omega_tau2x - omega_s2x
            if ( (omega_tmp  >   pi   ) )  omega_tmp = omega_tmp - 2.*pi
            if ( (omega_tmp  < (0.-pi)) )  omega_tmp = omega_tmp + 2.*pi
            Omega_tau2s(i,j,kp1) = omega_tmp

          enddo
        endif

        ! Stokes drift
        do b=1,waves%NumBands
          uS0(i,j)  = uS0(i,j) + waves%UStk_Hb(i,j,b)    ! or forces%UStkb(i,j,b)
          vS0(i,j)  = vS0(i,j) + waves%VStk_Hb(i,j,b)    ! or forces%VStkb(i,j,b)
        enddo
        depth = 0.0
        do k = 1,nz
          do b  = 1, waves%NumBands
            ! cell center
            fexp = exp(-2. * waves%WaveNum_Cen(b) * (depth+0.5*h(i,j,k)) )
            uStk(i,j,k) = uStk(i,j,k) + waves%UStk_Hb(i,j,b) * fexp
            vStk(i,j,k) = vStk(i,j,k) + waves%VStk_Hb(i,j,b) * fexp
          enddo
          ! cell bottom
          depth = depth + h(i,j,k)
        enddo
      enddo
    enddo

    ! post FPmix diagnostics
    if (CS%id_uE_h    > 0) call post_data(CS%id_uE_h     , uE_h   , CS%diag)
    if (CS%id_vE_h    > 0) call post_data(CS%id_vE_h   , vE_h   , CS%diag)
    if (CS%id_uInc_h  > 0) call post_data(CS%id_uInc_h , uInc_h , CS%diag)
    if (CS%id_vInc_h  > 0) call post_data(CS%id_vInc_h , vInc_h , CS%diag)
    if (CS%id_FPtau2s > 0) call post_data(CS%id_FPtau2s, Omega_tau2s, CS%diag)
    if (CS%id_FPtau2w > 0) call post_data(CS%id_FPtau2w, Omega_tau2w, CS%diag)
    if (CS%id_uStk0   > 0) call post_data(CS%id_uStk0  , uS0 , CS%diag)
    if (CS%id_vStk0   > 0) call post_data(CS%id_vStk0  , vS0    , CS%diag)
    if (CS%id_uStk    > 0) call post_data(CS%id_uStk   , uStk   , CS%diag)
    if (CS%id_vStk    > 0) call post_data(CS%id_vStk   , vStk   , CS%diag)
    if (CS%id_Omega_w2x > 0) call post_data(CS%id_Omega_w2x, forces%omega_w2x, CS%diag)

  endif

end subroutine vertFPmix


!> Compute coupling coefficient associated with vertical viscosity parameterization as in Greatbatch and Lamb
!! (1990), hereafter referred to as the GL90 vertical viscosity parameterization. This vertical viscosity scheme
!! redistributes momentum in the vertical, and is the equivalent of the Gent & McWilliams (1990) parameterization,
!! but in a TWA (thickness-weighted averaged) set of equations. The vertical viscosity coefficient nu is computed
!! from kappa_GM via thermal wind balance, and the following relation:
!! nu = kappa_GM * f^2 / N^2.
!! In the following subroutine kappa_GM is assumed either (a) constant or (b) horizontally varying. In both cases,
!! (a) and  (b), one can additionally impose an EBT structure in the vertical for kappa_GM.
!! A third possible formulation of nu is depth-independent:
!! nu = f^2 * alpha
!! The latter formulation would be equivalent to a kappa_GM that varies as N^2 with depth.
!! The vertical viscosity del_z ( nu del_z u) is applied to the momentum equation with stress-free boundary
!! conditions at the top and bottom.
!!
!! In SSW mode, we have 1/N^2 = h/g'. The coupling coefficient is therefore equal to
!! a_cpl_gl90 = nu / h = kappa_GM * f^2 / g'
!! or
!! a_cpl_gl90 = nu / h = f^2 * alpha / h

subroutine find_coupling_coef_gl90(a_cpl_gl90, hvel, do_ij, z_i, G, GV, CS, VarMix, work_on_u, &
                                   isb, ieb, jsb, jeb, nii, njj)
  type(ocean_grid_type), intent(in) :: G        !< Grid structure.
  type(verticalGrid_type), intent(in) :: GV     !< Vertical grid structure.
  integer, intent(in) :: nii                   !< Zonal block storage extent
  integer, intent(in) :: njj                   !< Meridional block storage extent
  real, dimension(nii,njj,SZK_(GV)+1),intent(inout) :: a_cpl_gl90 !< Coupling coefficient associated
                                                !! with GL90 across interfaces; is not
                                                !! included in a_cpl [H T-1 ~> m s-1 or Pa s m-1].
  real, dimension(nii,njj,SZK_(GV)), intent(in) :: hvel !< Distance between interfaces
                                                !! at velocity points [Z ~> m]
  logical, intent(in) :: do_ij(nii,njj)         !< If true, calculate this column
  real, dimension(nii,njj,SZK_(GV)+1), intent(in) :: z_i !< Estimate of interface heights above the
                                                !! bottom, normalized by the GL90 bottom
                                                !! boundary layer thickness [nondim]
  type(vertvisc_cs), intent(in) :: CS           !< Vertical viscosity control structure
  type(VarMix_CS), intent(in) :: VarMix         !< Variable mixing coefficients
  logical, intent(in) :: work_on_u              !< If true, u-points are being calculated,
                                                !! otherwise they are v-points.
  integer, intent(in) :: isb                    !< Start of the block's i index range
  integer, intent(in) :: ieb                    !< End of the block's i index range
  integer, intent(in) :: jsb                    !< Start of the block's j index range
  integer, intent(in) :: jeb                    !< End of the block's j index range

  ! local variables
  logical :: kdgl90_use_vert_struct  ! use vertical structure for GL90 coefficient
  integer :: i, j, k, nz, ii, jj, iie, jje
  real    :: f2         !< Squared Coriolis parameter at a velocity grid point [T-2 ~> s-2].
  real    :: h_neglect  ! A vertical distance that is so small it is usually lost in roundoff error
                        ! and can be neglected [Z ~> m].
  real    :: botfn      ! A function that is 1 at the bottom and small far from it [nondim]
  real    :: z2         ! The distance from the bottom, normalized by Hbbl_gl90 [nondim]

  nz = GV%ke
  iie = ieb - isb + 1
  jje = jeb - jsb + 1
  h_neglect = GV%dZ_subroundoff
  kdgl90_use_vert_struct = .false.

  if (VarMix%use_variable_mixing) then
    kdgl90_use_vert_struct = allocated(VarMix%kdgl90_struct)
  endif

  do concurrent (K=1:nz+1, jj=1:jje, ii=1:iie, do_ij(ii,jj))
    a_cpl_gl90(ii,jj,K) = 0.
  enddo

  do concurrent (K=2:nz, jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(i,j,f2,z2,botfn))
    i = isb + ii - 1
    j = jsb + jj - 1

    if (work_on_u) then
      ! compute coupling coefficient at u-points
      f2 = 0.25 * (G%CoriolisBu(I,J-1) + G%CoriolisBu(I,J))**2
      if (CS%use_GL90_N2) then
        a_cpl_gl90(ii,jj,K) = 2. * f2 * CS%alpha_gl90 / (hvel(ii,jj,k) + hvel(ii,jj,k-1) + h_neglect)
      else
        if (CS%read_kappa_gl90) then
          a_cpl_gl90(ii,jj,K) = f2 * 0.5 * (CS%kappa_gl90_2d(i,j) + CS%kappa_gl90_2d(i+1,j)) / GV%g_prime(K)
        else
          a_cpl_gl90(ii,jj,K) = f2 * CS%kappa_gl90 / GV%g_prime(K)
        endif
        if (kdgl90_use_vert_struct) then
          a_cpl_gl90(ii,jj,K) = a_cpl_gl90(ii,jj,K) * 0.5 &
              * (VarMix%kdgl90_struct(i,j,k-1) + VarMix%kdgl90_struct(i+1,j,k-1))
        endif
      endif
      ! botfn determines when a point is within the influence of the GL90 bottom boundary layer,
      ! going from 1 at the bottom to 0 in the interior.
      z2 = z_i(ii,jj,k)
      botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

      a_cpl_gl90(ii,jj,K) = a_cpl_gl90(ii,jj,K) * (1. - botfn)
    else
      ! compute viscosities at v-points
      f2 = 0.25 * (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J))**2

      if (CS%use_GL90_N2) then
        a_cpl_gl90(ii,jj,K) = 2. * f2 * CS%alpha_gl90 / (hvel(ii,jj,k) + hvel(ii,jj,k-1) + h_neglect)
      else
        if (CS%read_kappa_gl90) then
          a_cpl_gl90(ii,jj,K) = f2 * 0.5 * (CS%kappa_gl90_2d(i,j) + CS%kappa_gl90_2d(i,j+1)) / GV%g_prime(K)
        else
          a_cpl_gl90(ii,jj,K) = f2 * CS%kappa_gl90 / GV%g_prime(K)
        endif
        if (kdgl90_use_vert_struct) then
          a_cpl_gl90(ii,jj,K) = a_cpl_gl90(ii,jj,K) * 0.5 &
              * (VarMix%kdgl90_struct(i,j,k-1) + VarMix%kdgl90_struct(i,j+1,k-1))
        endif
      endif
      ! botfn determines when a point is within the influence of the GL90 bottom boundary layer,
      ! going from 1 at the bottom to 0 in the interior.
      z2 = z_i(ii,jj,k)
      botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

      a_cpl_gl90(ii,jj,K) = a_cpl_gl90(ii,jj,K) * (1. - botfn)
    endif
  enddo
end subroutine find_coupling_coef_gl90


!> Perform a fully implicit vertical diffusion
!! of momentum.  Stress top and bottom boundary conditions are used.
!!
!! This is solving the tridiagonal system
!! \f[ \left(h_k + a_{k + 1/2} + a_{k - 1/2} + r_k\right) u_k^{n+1}
!!     = h_k u_k^n + a_{k + 1/2} u_{k+1}^{n+1} + a_{k - 1/2} u_{k-1}^{n+1} \f]
!! where \f$a_{k + 1/2} = \Delta t \nu_{k + 1/2} / h_{k + 1/2}\f$
!! is the <em>interfacial coupling thickness per time step</em>,
!! encompassing background viscosity as well as contributions from
!! enhanced mixed and bottom layer viscosities.
!! $r_k$ is a Rayleigh drag term due to channel drag.
!! There is an additional stress term on the right-hand side
!! if DIRECT_STRESS is true, applied to the surface layer.
subroutine vertvisc(u, v, h, forces, visc, dt, OBC, ADp, CDp, G, GV, US, CS, &
                    taux_bot, tauy_bot, fpmix, Waves, return_bottom_stress)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(mech_forcing),    intent(in)      :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),   intent(inout)   :: visc   !< Viscosities and bottom drag
  real,                  intent(in)      :: dt     !< Time increment [T ~> s]
  type(ocean_OBC_type),  pointer         :: OBC    !< Open boundary condition structure
  type(accel_diag_ptrs), intent(inout)   :: ADp    !< Accelerations in the momentum
                                                   !! equations for diagnostics
  type(cont_diag_ptrs),  intent(inout)   :: CDp    !< Continuity equation terms
  type(vertvisc_CS),     intent(inout)   :: CS     !< Vertical viscosity control structure
  real, dimension(SZIB_(G),SZJ_(G)), &
                   optional, intent(out) :: taux_bot !< Zonal bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                   optional, intent(out) :: tauy_bot !< Meridional bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  logical,         optional, intent(in)  :: fpmix !< fpmix along Eulerian shear
  type(wave_parameters_CS), &
                   optional, pointer     :: Waves !< Container for wave/Stokes information
  logical, optional, intent(in) :: return_bottom_stress !< If false, leave bottom stresses on device
                                                      !! except for requested diagnostics. Default true.

  ! Fields from forces used in this subroutine:
  !   taux: Zonal wind stress [R L Z T-2 ~> Pa].
  !   tauy: Meridional wind stress [R L Z T-2 ~> Pa].

  ! Local variables

  real :: Hmix             ! The mixed layer thickness over which stress
                           ! is applied with direct_stress [H ~> m or kg m-2].
  real :: I_Hmix           ! The inverse of Hmix [H-1 ~> m-1 or m2 kg-1].
  real :: Idt              ! The inverse of the time step [T-1 ~> s-1].
  real :: dt_Rho0          ! The time step divided by the mean density [T H Z-1 R-1 ~> s m3 kg-1 or s].
  real :: h_neglect        ! A thickness that is so small it is usually lost
                           ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: stress           !   The surface stress times the time step, divided
                           ! by the density [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: accel_underflow  ! An acceleration magnitude that is so small that values that are less
                           ! than this are diagnosed as 0 [L T-2 ~> m s-2].
  real :: zDS              ! Accumulated thickness used with direct_stress [H ~> m or kg m-2]
  real :: h_a              ! Layer thickness used with direct_stress [H ~> m or kg m-2]
  real :: hfr              ! Temporary ratio of thicknesses used with direct_stress [nondim]
  real :: surface_stress(SZIB_(G), SZJB_(G))
    ! The same as stress, unless the wind stress is applied as a body force
    ! [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, allocatable, dimension(:,:,:) :: KE_term ! A term in the kinetic energy budget
                                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:) :: KE_u ! The area integral of a KE term in a layer at u-points
                                              ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real, allocatable, dimension(:,:,:) :: KE_v ! The area integral of a KE term in a layer at v-points
                                              ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]

  logical :: DoStokesMixing
  logical :: lfpmix
  logical :: host_bottom_stress ! Whether the caller needs host bottom stresses.

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  integer :: nIIB, njj, nii, nJJB ! Resolved horizontal block sizes at u- and v-points

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  if (CS%id_GLwork > 0) then
    allocate(KE_u(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke))
    allocate(KE_v(G%isd:G%ied,G%JsdB:G%JedB,GV%ke))
    allocate(KE_term(G%isd:G%ied,G%jsd:G%jed,GV%ke))
    !$omp target enter data map(alloc: KE_u, KE_v, KE_term)
    do concurrent (k=1:nz, j=G%jsd:G%jed, I=G%IsdB:G%IedB)
      KE_u(I,j,k) = 0.0
    enddo
    do concurrent (k=1:nz, J=G%JsdB:G%JedB, i=G%isd:G%ied)
      KE_v(i,J,k) = 0.0
    enddo
    do concurrent (k=1:nz, j=G%jsd:G%jed, i=G%isd:G%ied)
      KE_term(i,j,k) = 0.0
    enddo
    if (.not.G%symmetric) &
      call create_group_pass(CS%pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
  endif

  if (CS%direct_stress) then
    Hmix = CS%Hmix_stress
    I_Hmix = 1.0 / Hmix
  endif
  dt_Rho0 = dt / GV%H_to_RZ
  h_neglect = GV%H_subroundoff
  Idt = 1.0 / dt

  accel_underflow = CS%vel_underflow * Idt

  !Check if Stokes mixing allowed if requested (present and associated)
  DoStokesMixing=.false.
  if (CS%StokesMixing) then
    if (present(Waves)) DoStokesMixing = associated(Waves)
    if (.not. DoStokesMixing) &
      call MOM_error(FATAL, "Stokes Mixing called without associated Waves Control Structure")
  endif
  lfpmix = .false.
  if ( present(fpmix) ) lfpmix = fpmix
  host_bottom_stress = .true.
  if (present(return_bottom_stress)) host_bottom_stress = return_bottom_stress

  ! Keep forcing, tendencies and solver work resident until the host limiter/reporting boundary.
  !$omp target data map(to: u, v, h, forces)
  !$omp target enter data map(alloc: surface_stress)
  !$omp target enter data map(to: forces%taux, forces%tauy)
  if (DoStokesMixing .or. lfpmix) then
    !$omp target enter data map(to: Waves)
    !$omp target enter data map(to: Waves%Us_x, Waves%Us_y)
  endif
  !   Update the zonal velocity component using a modification of a standard
  ! tridiagonal solver.

  ! WGL: Brandon Reichl says the following is obsolete. u(I,j,k) already
  ! includes Stokes.
  ! When mixing down Eulerian current + Stokes drift add before calling solver
  if (DoStokesMixing) then
    do concurrent (k=1:nz, j=G%jsc:G%jec, I=Isq:Ieq, G%mask2dCu(I,j) > 0.)
      u(I,j,k) = u(I,j,k) + Waves%Us_x(I,j,k)
    enddo
  endif

  if (lfpmix) then
    do concurrent (k=1:nz, j=G%jsc:G%jec, I=Isq:Ieq, G%mask2dCu(I,j) > 0.)
      u(I,j,k) = u(I,j,k) - Waves%Us_x(I,j,k)
    enddo
  endif

  if (associated(ADp%du_dt_visc)) then
    do concurrent (k=1:nz, j=G%jsc:G%jec, I=Isq:Ieq)
      ADp%du_dt_visc(I,j,k) = u(I,j,k)
    enddo
  endif

  if (associated(ADp%du_dt_visc_gl90)) then
    do concurrent (k=1:nz, j=G%jsc:G%jec, I=Isq:Ieq)
      ADp%du_dt_visc_gl90(I,j,k) = u(I,j,k)
    enddo
  endif

  if (associated(ADp%du_dt_str)) then
    do concurrent (k=1:nz, j=G%jsc:G%jec, I=Isq:Ieq)
      ADp%du_dt_str(I,j,k) = 0.0
    enddo
  endif

  !   One option is to have the wind stress applied as a body force
  ! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
  ! the wind stress is applied as a stress boundary condition.
  if (CS%direct_stress) then
    do concurrent (j=G%jsc:G%jec, I=Isq:Ieq, G%mask2dCu(I,j) > 0.) &
        DO_LOCALITY(local(k,zDS,stress,h_a,hfr))
      surface_stress(I,j) = 0.0
      zDS = 0.0
      stress = dt_Rho0 * forces%taux(I,j)
      do k=1,nz
        h_a = 0.5 * (h(i,j,k) + h(i+1,j,k)) + h_neglect
        hfr = 1.0 ; if ((zDS+h_a) > Hmix) hfr = (Hmix - zDS) / h_a
        u(I,j,k) = u(I,j,k) + I_Hmix * hfr * stress
        if (associated(ADp%du_dt_str)) ADp%du_dt_str(I,j,k) = (I_Hmix * hfr * stress) * Idt
        zDS = zDS + h_a ; if (zDS >= Hmix) exit
      enddo
    enddo
  else
    do concurrent (j=G%jsc:G%jec, I=Isq:Ieq)
      surface_stress(I,j) = dt_Rho0 * (G%mask2dCu(I,j)*forces%taux(I,j))
    enddo
  endif

  ! perform forward elimination on the tridiagonal system
  !
  ! denote the diagonal of the system as b_k, the subdiagonal as a_k
  ! and the superdiagonal as c_k. The right-hand side terms are d_k.
  !
  ! ignoring the Rayleigh drag contribution,
  ! we have a_k = -dt * a_u(k)
  !         b_k = h_u(k) + dt * (a_u(k) + a_u(k+1))
  !         c_k = -dt * a_u(k+1)
  !
  ! for forward elimination, we want to:
  ! calculate c'_k = - c_k                / (b_k + a_k c'_(k-1))
  ! and       d'_k = (d_k - a_k d'_(k-1)) / (b_k + a_k c'_(k-1))
  ! where c'_1 = c_1/b_1 and d'_1 = d_1/b_1
  !
  ! This form is mathematically equivalent to Thomas' tridiagonal matrix algorithm, but it
  ! does not suffer from the acute sensitivity to truncation errors of the Thomas algorithm
  ! because it involves no subtraction, as discussed by Schopf & Loughe, MWR, 1995.
  !
  ! b1 is the denominator term 1 / (b_k + a_k c'_(k-1))
  ! b_denom_1 is (b_k + a_k + c_k) - a_k(1 - c'_(k-1))
  !            = (b_k + c_k + c'_(k-1))
  ! this is done so that d1 = b1 * b_denom_1 = 1 - c'_(k-1)
  ! c1(k) is -c'_(k - 1)
  ! and the right-hand-side is destructively updated to be d'_k

  nIIB = CS%niblock
  if (nIIB == 0) nIIB = Ieq - Isq + 1
  njj = CS%njblock
  if (njj == 0) njj = je - js + 1
  call vertvisc_solve_u(u, surface_stress, visc, ADp, dt, Idt, accel_underflow, G, GV, CS, nIIB, njj)

  ! compute vertical velocity tendency that arises from GL90 viscosity;
  ! follow tridiagonal solve method as above; to avoid corrupting u,
  ! use ADp%du_dt_visc_gl90 as a placeholder for updated u (due to GL90) until last do loop
  if ((CS%id_du_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    if (associated(ADp%du_dt_visc_gl90)) then
      call vertvisc_gl90_u(ADp, KE_u, dt, Idt, accel_underflow, (CS%id_GLwork > 0), G, GV, CS, nIIB, njj)
    endif
  endif

  if (associated(ADp%du_dt_visc)) then
    do concurrent (k=1:nz, j=G%jsc:G%jec, I=Isq:Ieq)
      ADp%du_dt_visc(I,j,k) = (u(I,j,k) - ADp%du_dt_visc(I,j,k)) * Idt

      if (abs(ADp%du_dt_visc(I,j,k)) < accel_underflow) &
        ADp%du_dt_visc(I,j,k) = 0.0
    enddo
  endif

  if (allocated(visc%taux_shelf)) then
    do concurrent (j=G%jsc:G%jec, I=Isq:Ieq)
      visc%taux_shelf(I,j) = -GV%H_to_RZ * CS%a1_shelf_u(I,j) * u(I,j,1) ! - u_shelf?
    enddo
  endif

  if (present(taux_bot)) then
    do concurrent (j=G%jsc:G%jec, I=Isq:Ieq)
      taux_bot(I,j) = GV%H_to_RZ * (u(I,j,nz) * CS%a_u(I,j,nz+1))
    enddo

    if (allocated(visc%Ray_u)) then
      do concurrent (j=G%jsc:G%jec) DO_LOCALITY(local(k))
        do k=1,nz
          do concurrent (I=Isq:Ieq)
            taux_bot(I,j) = taux_bot(I,j) + GV%H_to_RZ * (visc%Ray_u(I,j,k) * u(I,j,k))
          enddo
        enddo
      enddo
    endif
  endif

  ! When mixing down Eulerian current + Stokes drift subtract after calling solver
  if (DoStokesMixing) then
    do concurrent (k=1:nz, j=G%jsc:G%jec, I=Isq:Ieq, G%mask2dCu(I,j) > 0.)
      u(I,j,k) = u(I,j,k) - Waves%Us_x(I,j,k)
    enddo
  endif

  if (lfpmix) then
    do concurrent (k=1:nz, j=G%jsc:G%jec, I=Isq:Ieq, G%mask2dCu(I,j) > 0.)
      u(I,j,k) = u(I,j,k) + Waves%Us_x(I,j,k)
    enddo
  endif

  ! == Now work on the meridional velocity component.

  ! When mixing down Eulerian current + Stokes drift add before calling solver
  if (DoStokesMixing) then
    do concurrent (k=1:nz, J=Jsq:Jeq, i=is:ie, G%mask2dCv(i,J) > 0.)
      v(i,j,k) = v(i,j,k) + Waves%Us_y(i,j,k)
    enddo
  endif

  if (lfpmix) then
    do concurrent (k=1:nz, J=Jsq:Jeq, i=is:ie, G%mask2dCv(i,J) > 0.)
      v(i,j,k) = v(i,j,k) - Waves%Us_y(i,j,k)
    enddo
  endif

  if (associated(ADp%dv_dt_visc)) then
    do concurrent (k=1:nz, J=Jsq:Jeq, i=is:ie)
      ADp%dv_dt_visc(i,J,k) = v(i,J,k)
    enddo
  endif

  if (associated(ADp%dv_dt_visc_gl90)) then
    do concurrent (k=1:nz, J=Jsq:Jeq, i=is:ie)
      ADp%dv_dt_visc_gl90(i,J,k) = v(i,J,k)
    enddo
  endif

  if (associated(ADp%dv_dt_str)) then
    do concurrent (k=1:nz, J=Jsq:Jeq, i=is:ie)
      ADp%dv_dt_str(i,J,k) = 0.0
    enddo
  endif

  !   One option is to have the wind stress applied as a body force
  ! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
  ! the wind stress is applied as a stress boundary condition.
  if (CS%direct_stress) then
    do concurrent (J=Jsq:Jeq, i=is:ie, G%mask2dCv(i,J) > 0.) &
        DO_LOCALITY(local(k,zDS,stress,h_a,hfr))
      surface_stress(i,J) = 0.0
      zDS = 0.0
      stress = dt_Rho0 * forces%tauy(i,J)
      do k=1,nz
        h_a = 0.5 * (h(i,J,k) + h(i,J+1,k)) + h_neglect
        hfr = 1.0 ; if ((zDS+h_a) > Hmix) hfr = (Hmix - zDS) / h_a
        v(i,J,k) = v(i,J,k) + I_Hmix * hfr * stress
        if (associated(ADp%dv_dt_str)) ADp%dv_dt_str(i,J,k) = (I_Hmix * hfr * stress) * Idt
        zDS = zDS + h_a ; if (zDS >= Hmix) exit
      enddo
    enddo
  else
    do concurrent (J=Jsq:Jeq, i=is:ie)
      surface_stress(i,J) = dt_Rho0 * (G%mask2dCv(i,J) * forces%tauy(i,J))
    enddo
  endif

  nii = CS%niblock
  if (nii == 0) nii = ie - is + 1
  nJJB = CS%njblock
  if (nJJB == 0) nJJB = Jeq - Jsq + 1
  call vertvisc_solve_v(v, surface_stress, visc, ADp, dt, Idt, accel_underflow, G, GV, CS, nii, nJJB)

  ! compute vertical velocity tendency that arises from GL90 viscosity;
  ! follow tridiagonal solve method as above; to avoid corrupting v,
  ! use ADp%dv_dt_visc_gl90 as a placeholder for updated v (due to GL90) until last do loop
  if ((CS%id_dv_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    if (associated(ADp%dv_dt_visc_gl90)) then
      call vertvisc_gl90_v(ADp, KE_v, dt, Idt, accel_underflow, (CS%id_GLwork > 0), G, GV, CS, nii, nJJB)
    endif
  endif

  if (associated(ADp%dv_dt_visc)) then
    do concurrent (k=1:nz, J=Jsq:Jeq, i=is:ie)
      ADp%dv_dt_visc(i,J,k) = (v(i,J,k) - ADp%dv_dt_visc(i,J,k))*Idt
      if (abs(ADp%dv_dt_visc(i,J,k)) < accel_underflow) ADp%dv_dt_visc(i,J,k) = 0.0
    enddo
  endif

  if (allocated(visc%tauy_shelf)) then
    do concurrent (J=Jsq:Jeq, i=is:ie)
      visc%tauy_shelf(i,J) = -GV%H_to_RZ * CS%a1_shelf_v(i,J) * v(i,J,1) ! - v_shelf?
    enddo
  endif

  if (present(tauy_bot)) then
    do concurrent (J=Jsq:Jeq, i=is:ie)
      tauy_bot(i,J) = GV%H_to_RZ * (v(i,J,nz) * CS%a_v(i,J,nz+1))
    enddo

    if (allocated(visc%Ray_v)) then
      do concurrent (J=Jsq:Jeq) DO_LOCALITY(local(k))
        do k=1,nz
          do concurrent (i=is:ie)
            tauy_bot(i,J) = tauy_bot(i,J) + GV%H_to_RZ * (visc%Ray_v(i,J,k)*v(i,J,k))
          enddo
        enddo
      enddo
    endif
  endif

  ! When mixing down Eulerian current + Stokes drift subtract after calling solver
  if (DoStokesMixing) then
    do concurrent (k=1:nz, J=Jsq:Jeq, i=is:ie, G%mask2dCv(i,J) > 0.)
      v(i,J,k) = v(i,J,k) - Waves%Us_y(i,J,k)
    enddo
  endif

  if (lfpmix) then
    do concurrent (k=1:nz, J=Jsq:Jeq, i=is:ie, G%mask2dCv(i,J) > 0.)
      v(i,J,k) = v(i,J,k) + Waves%Us_y(i,J,k)
    enddo
  endif

  ! Calculate the KE source from GL90 vertical viscosity [H L2 T-3 ~> m3 s-3].
  ! We do the KE-rate calculation here (rather than in MOM_diagnostics) to ensure
  ! a sign-definite term. MOM_diagnostics does not have access to the velocities
  ! and thicknesses used in the vertical solver, but rather uses a time-mean
  ! barotropic transport [uv]h.
  if (CS%id_GLwork > 0) then
    if (.not.G%symmetric) then
      !$omp target update from(KE_u, KE_v)
      call do_group_pass(CS%pass_KE_uv, G%domain)
      !$omp target update to(KE_u, KE_v)
    endif
    do concurrent (k=1:nz, j=js:je, i=is:ie)
      KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
          * (KE_u(I,j,k) + KE_u(I-1,j,k) + KE_v(i,J,k) + KE_v(i,J-1,k))
    enddo
    !$omp target update from(KE_term)
    call post_data(CS%id_GLwork, KE_term, CS%diag)
    !$omp target exit data map(release: KE_u, KE_v, KE_term)
  endif

  ! Acceleration reports and diagnostic posting consume host tendencies.
  if (associated(ADp%du_dt_visc)) then
    !$omp target update from(ADp%du_dt_visc)
  endif
  if (associated(ADp%dv_dt_visc)) then
    !$omp target update from(ADp%dv_dt_visc)
  endif
  if (associated(ADp%du_dt_visc_gl90)) then
    !$omp target update from(ADp%du_dt_visc_gl90)
  endif
  if (associated(ADp%dv_dt_visc_gl90)) then
    !$omp target update from(ADp%dv_dt_visc_gl90)
  endif
  if (associated(ADp%du_dt_str)) then
    !$omp target update from(ADp%du_dt_str)
  endif
  if (associated(ADp%dv_dt_str)) then
    !$omp target update from(ADp%dv_dt_str)
  endif
  if (present(taux_bot)) then
    if (host_bottom_stress .or. &
        (query_averaging_enabled(CS%diag) .and. (CS%id_taux_bot > 0))) then
      !$omp target update from(taux_bot)
    endif
  endif
  if (present(tauy_bot)) then
    if (host_bottom_stress .or. &
        (query_averaging_enabled(CS%diag) .and. (CS%id_tauy_bot > 0))) then
      !$omp target update from(tauy_bot)
    endif
  endif
  if (allocated(visc%taux_shelf)) then
    !$omp target update from(visc%taux_shelf)
  endif
  if (allocated(visc%tauy_shelf)) then
    !$omp target update from(visc%tauy_shelf)
  endif
  if (DoStokesMixing .or. lfpmix) then
    !$omp target exit data map(release: Waves%Us_x, Waves%Us_y)
    !$omp target exit data map(release: Waves)
  endif
  !$omp target exit data map(release: forces%taux, forces%tauy)
  !$omp target exit data map(release: surface_stress)

  call vertvisc_limit_vel(u, v, h, ADp, CDp, forces, visc, dt, G, GV, US, CS)

  ! Fence the device limiter before host OBC handling and subsequent host consumers.
  !$omp target update from(u, v)

  ! Here the velocities associated with open boundary conditions are applied.
  if (associated(OBC)) then
    do n=1,OBC%number_of_segments
      if (OBC%segment(n)%specified) then
        if (OBC%segment(n)%is_N_or_S) then
          J = OBC%segment(n)%HI%JsdB
          do k=1,nz ; do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
            v(i,J,k) = OBC%segment(n)%normal_vel(i,J,k)
          enddo ; enddo
        elseif (OBC%segment(n)%is_E_or_W) then
          I = OBC%segment(n)%HI%IsdB
          do k=1,nz ; do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
            u(I,j,k) = OBC%segment(n)%normal_vel(I,j,k)
          enddo ; enddo
        endif
      endif
    enddo
    !$omp target update to(u, v)
  endif
  !$omp end target data

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_du_dt_visc > 0) &
      call post_data(CS%id_du_dt_visc, ADp%du_dt_visc, CS%diag)
    if (CS%id_du_dt_visc_gl90 > 0) &
      call post_data(CS%id_du_dt_visc_gl90, ADp%du_dt_visc_gl90, CS%diag)
    if (CS%id_dv_dt_visc > 0) &
      call post_data(CS%id_dv_dt_visc, ADp%dv_dt_visc, CS%diag)
    if (CS%id_dv_dt_visc_gl90 > 0) &
      call post_data(CS%id_dv_dt_visc_gl90, ADp%dv_dt_visc_gl90, CS%diag)
    if (present(taux_bot) .and. (CS%id_taux_bot > 0)) &
      call post_data(CS%id_taux_bot, taux_bot, CS%diag)
    if (present(tauy_bot) .and. (CS%id_tauy_bot > 0)) &
      call post_data(CS%id_tauy_bot, tauy_bot, CS%diag)
    if (CS%id_du_dt_str > 0) &
      call post_data(CS%id_du_dt_str, ADp%du_dt_str, CS%diag)
    if (CS%id_dv_dt_str > 0) &
      call post_data(CS%id_dv_dt_str, ADp%dv_dt_str, CS%diag)

    if (associated(ADp%du_dt_visc) .and. associated(ADp%dv_dt_visc)) then
      ! Diagnostics of the fractional thicknesses times momentum budget terms
      ! 3D diagnostics of hf_du(dv)_dt_visc are commented because there is no clarity on proper remapping grid option.
      ! The code is retained for debugging purposes in the future.
      !if (CS%id_hf_du_dt_visc > 0) &
      !  call post_product_u(CS%id_hf_du_dt_visc, ADp%du_dt_visc, ADp%diag_hfrac_u, G, nz, CS%diag)
      !if (CS%id_hf_dv_dt_visc > 0) &
      !  call post_product_v(CS%id_hf_dv_dt_visc, ADp%dv_dt_visc, ADp%diag_hfrac_v, G, nz, CS%diag)

      ! Diagnostics for thickness-weighted vertically averaged viscous accelerations
      if (CS%id_hf_du_dt_visc_2d > 0) &
        call post_product_sum_u(CS%id_hf_du_dt_visc_2d, ADp%du_dt_visc, ADp%diag_hfrac_u, G, nz, CS%diag)
      if (CS%id_hf_dv_dt_visc_2d > 0) &
        call post_product_sum_v(CS%id_hf_dv_dt_visc_2d, ADp%dv_dt_visc, ADp%diag_hfrac_v, G, nz, CS%diag)

      ! Diagnostics for thickness x viscous accelerations
      if (CS%id_h_du_dt_visc > 0) call post_product_u(CS%id_h_du_dt_visc, ADp%du_dt_visc, ADp%diag_hu, G, nz, CS%diag)
      if (CS%id_h_dv_dt_visc > 0) call post_product_v(CS%id_h_dv_dt_visc, ADp%dv_dt_visc, ADp%diag_hv, G, nz, CS%diag)
    endif

    if (associated(ADp%du_dt_str) .and.  associated(ADp%dv_dt_str)) then
      ! Diagnostics for thickness x wind stress accelerations
      if (CS%id_h_du_dt_str > 0) call post_product_u(CS%id_h_du_dt_str, ADp%du_dt_str, ADp%diag_hu, G, nz, CS%diag)
      if (CS%id_h_dv_dt_str > 0) call post_product_v(CS%id_h_dv_dt_str, ADp%dv_dt_str, ADp%diag_hv, G, nz, CS%diag)

      ! Diagnostics for wind stress accelerations multiplied by visc_rem_[uv],
      if (CS%id_du_dt_str_visc_rem > 0) &
        call post_product_u(CS%id_du_dt_str_visc_rem, ADp%du_dt_str, ADp%visc_rem_u, G, nz, CS%diag)
      if (CS%id_dv_dt_str_visc_rem > 0) &
        call post_product_v(CS%id_dv_dt_str_visc_rem, ADp%dv_dt_str, ADp%visc_rem_v, G, nz, CS%diag)
    endif
  endif
end subroutine vertvisc


!> Solve vertical viscosity at zonal velocity points in horizontal blocks.
!! Forward elimination and back substitution retain their serial vertical order.
subroutine vertvisc_solve_u(u, surface_stress, visc, ADp, dt, Idt, accel_underflow, G, GV, CS, nIIB, njj)
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  real, intent(inout) :: u(SZIB_(G),SZJ_(G),SZK_(GV))
    !< Zonal velocity [L T-1 ~> m s-1]
  real, intent(in) :: surface_stress(SZIB_(G),SZJB_(G))
    !< Applied surface stress times the time step divided by density
    !! [H L T-1 ~> m2 s-1 or kg m-1 s-1]
  type(vertvisc_type), intent(in) :: visc
    !< Viscosities and bottom drag
  type(accel_diag_ptrs), intent(in) :: ADp
    !< Acceleration diagnostic pointers
  real, intent(in) :: dt
    !< Time step [T ~> s]
  real, intent(in) :: Idt
    !< Inverse time step [T-1 ~> s-1]
  real, intent(in) :: accel_underflow
    !< Acceleration cutoff [L T-2 ~> m s-2]
  type(vertvisc_CS), intent(in) :: CS
    !< Vertical viscosity control structure
  integer, intent(in) :: nIIB
    !< Zonal block size at zonal velocity points
  integer, intent(in) :: njj
    !< Meridional block size at zonal velocity points

  real :: c1(nIIB,njj,SZK_(GV))
    ! Tridiagonal back substitution coefficient [nondim]
  real :: b1(nIIB,njj)
    ! Inverse denominator [H-1 ~> m-1 or m2 kg-1]
  real :: d1(nIIB,njj)
    ! Elimination factor [nondim]
  real :: Ray
    ! Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
  real :: b_denom_1
    ! First term in the denominator of b1 [H ~> m or kg m-2]
  logical :: do_ij(nIIB,njj)
    ! Wet-point mask
  logical :: have_Ray, have_str
    ! Availability of Rayleigh drag and stress diagnostics
  integer :: I, j, k, IsbB, IebB, jsb, jeb, II, jj, IIe, jje, nz
    ! Global and block-local indices and bounds

  nz = GV%ke
  have_Ray = allocated(visc%Ray_u)
  have_str = associated(ADp%du_dt_str)

  !$omp target enter data map(alloc: c1, b1, d1, do_ij)

  do jsb=G%jsc,G%jec,njj ; do IsbB=G%IscB,G%IecB,nIIB
    IebB = min(IsbB + nIIB - 1, G%IecB)
    jeb = min(jsb + njj - 1, G%jec)
    IIe = IebB - IsbB + 1
    jje = jeb - jsb + 1

    do concurrent (jj=1:jje, II=1:IIe) DO_LOCALITY(local(I,j,Ray,b_denom_1))
      I = IsbB + II - 1
      j = jsb + jj - 1

      do_ij(II,jj) = (G%mask2dCu(I,j) > 0.)
      if (do_ij(II,jj)) then
        Ray = 0.
        if (have_Ray) Ray = visc%Ray_u(I,j,1)
        b_denom_1 = CS%h_u(I,j,1) + dt * (Ray + CS%a_u(I,j,1))
        b1(II,jj) = 1. / (b_denom_1 + dt * CS%a_u(I,j,2))
        d1(II,jj) = b_denom_1 * b1(II,jj)
        u(I,j,1) = b1(II,jj) * (CS%h_u(I,j,1) * u(I,j,1) + surface_stress(I,j))
        if (have_str) &
          ADp%du_dt_str(I,j,1) = b1(II,jj) * (CS%h_u(I,j,1) * ADp%du_dt_str(I,j,1) + surface_stress(I,j) * Idt)
      endif
    enddo

    do concurrent (jj=1:jje) DO_LOCALITY(local(k))
      do k=2,nz
        do concurrent (II=1:IIe) DO_LOCALITY(local(I,j,Ray,b_denom_1))
          I = IsbB + II - 1
          j = jsb + jj - 1

          if (do_ij(II,jj)) then
            Ray = 0.
            if (have_Ray) Ray = visc%Ray_u(I,j,k)
            c1(II,jj,k) = dt * CS%a_u(I,j,K) * b1(II,jj)
            b_denom_1 = CS%h_u(I,j,k) + dt * (Ray + CS%a_u(I,j,K) * d1(II,jj))
            b1(II,jj) = 1. / (b_denom_1 + dt * CS%a_u(I,j,K+1))
            d1(II,jj) = b_denom_1 * b1(II,jj)
            u(I,j,k) = (CS%h_u(I,j,k) * u(I,j,k) + dt * CS%a_u(I,j,K) * u(I,j,k-1)) * b1(II,jj)
            if (have_str) &
              ADp%du_dt_str(I,j,k) = (CS%h_u(I,j,k) * ADp%du_dt_str(I,j,k) &
                  + dt * CS%a_u(I,j,K) * ADp%du_dt_str(I,j,k-1)) * b1(II,jj)
          endif
        enddo
        !### Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
        if (dt < 0) exit
      enddo

      if (have_str) then
        do concurrent (II=1:IIe) DO_LOCALITY(local(I,j))
          I = IsbB + II - 1
          j = jsb + jj - 1

          if (do_ij(II,jj)) then
            if (abs(ADp%du_dt_str(I,j,nz)) < accel_underflow) ADp%du_dt_str(I,j,nz) = 0.
          endif
        enddo
      endif

      do k=nz-1,1,-1
        do concurrent (II=1:IIe) DO_LOCALITY(local(I,j))
          I = IsbB + II - 1
          j = jsb + jj - 1

          if (do_ij(II,jj)) then
            u(I,j,k) = u(I,j,k) + c1(II,jj,k+1) * u(I,j,k+1)
            if (have_str) then
              ADp%du_dt_str(I,j,k) = ADp%du_dt_str(I,j,k) + c1(II,jj,k+1) * ADp%du_dt_str(I,j,k+1)
              if (abs(ADp%du_dt_str(I,j,k)) < accel_underflow) ADp%du_dt_str(I,j,k) = 0.0
            endif
          endif
        enddo
      enddo
    enddo
  enddo ; enddo
  !$omp target exit data map(release: c1, b1, d1, do_ij)
end subroutine vertvisc_solve_u


!> Solve vertical viscosity at meridional velocity points in horizontal blocks.
!! Forward elimination and back substitution retain their serial vertical order.
subroutine vertvisc_solve_v(v, surface_stress, visc, ADp, dt, Idt, accel_underflow, G, GV, CS, nii, nJJB)
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  real, intent(inout) :: v(SZI_(G),SZJB_(G),SZK_(GV))
    !< Meridional velocity [L T-1 ~> m s-1]
  real, intent(in) :: surface_stress(SZIB_(G),SZJB_(G))
    !< Applied surface stress times the time step divided by density
    !! [H L T-1 ~> m2 s-1 or kg m-1 s-1]
  type(vertvisc_type), intent(in) :: visc
    !< Viscosities and bottom drag
  type(accel_diag_ptrs), intent(in) :: ADp
    !< Acceleration diagnostic pointers
  real, intent(in) :: dt
    !< Time step [T ~> s]
  real, intent(in) :: Idt
    !< Inverse time step [T-1 ~> s-1]
  real, intent(in) :: accel_underflow
    !< Acceleration cutoff [L T-2 ~> m s-2]
  type(vertvisc_CS), intent(in) :: CS
    !< Vertical viscosity control structure
  integer, intent(in) :: nii
    !< Zonal block size at meridional velocity points
  integer, intent(in) :: nJJB
    !< Meridional block size at meridional velocity points

  real :: c1(nii,nJJB,SZK_(GV))
    ! Tridiagonal back substitution coefficient [nondim]
  real :: b1(nii,nJJB)
    ! Inverse denominator [H-1 ~> m-1 or m2 kg-1]
  real :: d1(nii,nJJB)
    ! Elimination factor [nondim]
  real :: Ray
    ! Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
  real :: b_denom_1
    ! First term in the denominator of b1 [H ~> m or kg m-2]
  logical :: do_ij(nii,nJJB)
    ! Wet-point mask
  logical :: have_Ray, have_str
    ! Availability of Rayleigh drag and stress diagnostics
  integer :: i, J, k, isb, ieb, JsbB, JebB, ii, JJ, iie, JJe, nz
    ! Global and block-local indices and bounds

  nz = GV%ke
  have_Ray = allocated(visc%Ray_v)
  have_str = associated(ADp%dv_dt_str)
  !$omp target enter data map(alloc: c1, b1, d1, do_ij)

  do JsbB=G%JscB,G%JecB,nJJB ; do isb=G%isc,G%iec,nii
    ieb = min(isb + nii - 1, G%iec)
    JebB = min(JsbB + nJJB - 1, G%JecB)
    iie = ieb - isb + 1
    JJe = JebB - JsbB + 1

    do concurrent (JJ=1:JJe, ii=1:iie) DO_LOCALITY(local(i,J,Ray,b_denom_1))
      i = isb + ii - 1
      J = JsbB + JJ - 1

      do_ij(ii,JJ) = (G%mask2dCv(i,J) > 0.)
      if (do_ij(ii,JJ)) then
        Ray = 0.
        if (have_Ray) Ray = visc%Ray_v(i,J,1)
        b_denom_1 = CS%h_v(i,J,1) + dt * (Ray + CS%a_v(i,J,1))
        b1(ii,JJ) = 1.0 / (b_denom_1 + dt*CS%a_v(i,J,2))
        d1(ii,JJ) = b_denom_1 * b1(ii,JJ)
        v(i,J,1) = b1(ii,JJ) * (CS%h_v(i,J,1) * v(i,J,1) + surface_stress(i,J))

        if (have_str) then
          ADp%dv_dt_str(i,J,1) = b1(ii,JJ) * (CS%h_v(i,J,1) * ADp%dv_dt_str(i,J,1) + surface_stress(i,J) * Idt)
        endif
      endif
    enddo

    do concurrent (JJ=1:JJe) DO_LOCALITY(local(k))
      do k=2,nz
        do concurrent (ii=1:iie) DO_LOCALITY(local(i,J,Ray,b_denom_1))
          i = isb + ii - 1
          J = JsbB + JJ - 1

          if (do_ij(ii,JJ)) then
            Ray = 0.
            if (have_Ray) Ray = visc%Ray_v(i,J,k)
            c1(ii,JJ,k) = dt * CS%a_v(i,J,K) * b1(ii,JJ)
            b_denom_1 = CS%h_v(i,J,k) + dt * (Ray + CS%a_v(i,J,K) * d1(ii,JJ))
            b1(ii,JJ) = 1. / (b_denom_1 + dt * CS%a_v(i,J,K+1))
            d1(ii,JJ) = b_denom_1 * b1(ii,JJ)
            v(i,J,k) = (CS%h_v(i,J,k) * v(i,J,k) + dt * CS%a_v(i,J,K) * v(i,J,k-1)) * b1(ii,JJ)
            if (have_str) &
              ADp%dv_dt_str(i,J,k) = (CS%h_v(i,J,k) * ADp%dv_dt_str(i,J,k) &
                  + dt * CS%a_v(i,J,K) * ADp%dv_dt_str(i,J,k-1)) * b1(ii,JJ)
          endif
        enddo

        !### Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
        if (dt < 0) exit
      enddo
      if (have_str) then
        do concurrent (ii=1:iie, do_ij(ii,JJ)) DO_LOCALITY(local(i,J))
          i = isb + ii - 1
          J = JsbB + JJ - 1

          if (abs(ADp%dv_dt_str(i,J,nz)) < accel_underflow) ADp%dv_dt_str(i,J,nz) = 0.0
        enddo
      endif
      do k=nz-1,1,-1
        do concurrent (ii=1:iie, do_ij(ii,JJ)) DO_LOCALITY(local(i,J))
          i = isb + ii - 1
          J = JsbB + JJ - 1

          v(i,J,k) = v(i,J,k) + c1(ii,JJ,k+1) * v(i,J,k+1)
          if (have_str) then
            ADp%dv_dt_str(i,J,k) = ADp%dv_dt_str(i,J,k) + c1(ii,JJ,k+1) * ADp%dv_dt_str(i,J,k+1)
            if (abs(ADp%dv_dt_str(i,J,k)) < accel_underflow) ADp%dv_dt_str(i,J,k) = 0.0
          endif
        enddo
      enddo
    enddo
  enddo ; enddo
  !$omp target exit data map(release: c1, b1, d1, do_ij)
end subroutine vertvisc_solve_v


!> Diagnose the GL90 acceleration and kinetic energy source at zonal points in horizontal blocks.
subroutine vertvisc_gl90_u(ADp, KE_u, dt, Idt, accel_underflow, do_KE, G, GV, CS, nIIB, njj)
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  type(accel_diag_ptrs), intent(in) :: ADp
    !< Acceleration diagnostic pointers; GL90 tendency initially holds the old velocity
  real, allocatable, intent(inout) :: KE_u(:,:,:)
    !< Area-integrated kinetic energy source [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real, intent(in) :: dt
    !< Time step [T ~> s]
  real, intent(in) :: Idt
    !< Inverse time step [T-1 ~> s-1]
  real, intent(in) :: accel_underflow
    !< Acceleration cutoff [L T-2 ~> m s-2]
  logical, intent(in) :: do_KE
    !< If true, calculate the kinetic energy source in KE_u
  type(vertvisc_CS), intent(in) :: CS
    !< Vertical viscosity control structure
  integer, intent(in) :: nIIB
    !< Zonal block size at zonal velocity points
  integer, intent(in) :: njj
    !< Meridional block size at zonal velocity points

  real :: c1(nIIB,njj,SZK_(GV))
    ! Back substitution coefficient [nondim]
  real :: b1(nIIB,njj)
    ! Inverse denominator [H-1 ~> m-1 or m2 kg-1]
  real :: d1(nIIB,njj)
    ! Elimination factor [nondim]
  real :: b_denom_1
    ! First term in the denominator of b1 [H ~> m or kg m-2]
  logical :: do_ij(nIIB,njj)
    ! Wet-point mask
  integer :: I, j, k, IsbB, IebB, jsb, jeb, II, jj, IIe, jje, nz
    ! Global and block-local indices and bounds

  nz = GV%ke
  !$omp target enter data map(alloc: c1, b1, d1, do_ij)

  do jsb=G%jsc,G%jec,njj ; do IsbB=G%IscB,G%IecB,nIIB
    IebB = min(IsbB + nIIB - 1, G%IecB)
    jeb = min(jsb + njj - 1, G%jec)
    IIe = IebB - IsbB + 1
    jje = jeb - jsb + 1
    do concurrent (jj=1:jje, II=1:IIe) DO_LOCALITY(local(I,j,b_denom_1))
      I = IsbB + II - 1
      j = jsb + jj - 1

      do_ij(II,jj) = (G%mask2dCu(I,j) > 0.)
      if (do_ij(II,jj)) then
        b_denom_1 = CS%h_u(I,j,1) ! CS%a_u_gl90(I,j,1) is zero.
        b1(II,jj) = 1.0 / (b_denom_1 + dt * CS%a_u_gl90(I,j,2))
        d1(II,jj) = b_denom_1 * b1(II,jj)
        ADp%du_dt_visc_gl90(I,j,1) = b1(II,jj) * (CS%h_u(I,j,1) * ADp%du_dt_visc_gl90(I,j,1))
      endif
    enddo
    do concurrent (jj=1:jje) DO_LOCALITY(local(k))
      do k=2,nz
        do concurrent (II=1:IIe) DO_LOCALITY(local(I,j,b_denom_1))
          I = IsbB + II - 1
          j = jsb + jj - 1

          if (do_ij(II,jj)) then
            c1(II,jj,k) = dt * CS%a_u_gl90(I,j,K) * b1(II,jj)
            b_denom_1 = CS%h_u(I,j,k) + dt * (CS%a_u_gl90(I,j,K)*d1(II,jj))
            b1(II,jj) = 1.0 / (b_denom_1 + dt * CS%a_u_gl90(I,j,K+1))
            d1(II,jj) = b_denom_1 * b1(II,jj)
            ADp%du_dt_visc_gl90(I,j,k) = (CS%h_u(I,j,k) * ADp%du_dt_visc_gl90(I,j,k) &
                + dt * CS%a_u_gl90(I,j,K) * ADp%du_dt_visc_gl90(I,j,k-1)) * b1(II,jj)
          endif
        enddo
      enddo
      ! Back substitute to solve for new velocities, held by ADp%du_dt_visc_gl90.
      do k=nz-1,1,-1
        do concurrent (II=1:IIe) DO_LOCALITY(local(I,j))
          I = IsbB + II - 1
          j = jsb + jj - 1

          if (do_ij(II,jj)) ADp%du_dt_visc_gl90(I,j,k) = &
              ADp%du_dt_visc_gl90(I,j,k) + c1(II,jj,k+1) * ADp%du_dt_visc_gl90(I,j,k+1)
        enddo
      enddo
    enddo ! jj-loop over rows within a block

    ! Now fill ADp%du_dt_visc_gl90(I,j,k) with the actual velocity tendency due to GL90.
    ! On the RHS, ADp%du_dt_visc(I,j,k) holds the original velocity u(I,j,k),
    ! and ADp%du_dt_visc_gl90(I,j,k) holds the updated velocity due to GL90.
    do concurrent (k=1:nz, jj=1:jje, II=1:IIe) DO_LOCALITY(local(I,j))
      I = IsbB + II - 1
      j = jsb + jj - 1

      if (do_ij(II,jj)) then
        ADp%du_dt_visc_gl90(I,j,k) = (ADp%du_dt_visc_gl90(I,j,k) - ADp%du_dt_visc(I,j,k)) * Idt
        if (abs(ADp%du_dt_visc_gl90(I,j,k)) < accel_underflow) ADp%du_dt_visc_gl90(I,j,k) = 0.0
      endif
    enddo
    if (do_KE) then
      ! To compute energetics, multiply by u*h, where u is the velocity before the update.
      ! ADp%du_dt_visc(I,j,k) still holds that original velocity u(I,j,k).
      do concurrent (k=1:nz, jj=1:jje, II=1:IIe) DO_LOCALITY(local(I,j))
        I = IsbB + II - 1
        j = jsb + jj - 1

        if (do_ij(II,jj)) &
          KE_u(I,j,k) = ADp%du_dt_visc(I,j,k) * CS%h_u(I,j,k) * G%areaCu(I,j) * ADp%du_dt_visc_gl90(I,j,k)
      enddo
    endif
  enddo ; enddo ! Horizontal blocks
  !$omp target exit data map(release: c1, b1, d1, do_ij)
end subroutine vertvisc_gl90_u


!> Diagnose the GL90 acceleration and kinetic energy source at meridional points in horizontal blocks.
subroutine vertvisc_gl90_v(ADp, KE_v, dt, Idt, accel_underflow, do_KE, G, GV, CS, nii, nJJB)
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  type(accel_diag_ptrs), intent(in) :: ADp
    !< Acceleration diagnostic pointers; GL90 tendency initially holds the old velocity
  real, allocatable, intent(inout) :: KE_v(:,:,:)
    !< Area-integrated kinetic energy source [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real, intent(in) :: dt
    !< Time step [T ~> s]
  real, intent(in) :: Idt
    !< Inverse time step [T-1 ~> s-1]
  real, intent(in) :: accel_underflow
    !< Acceleration cutoff [L T-2 ~> m s-2]
  logical, intent(in) :: do_KE
    !< If true, calculate the kinetic energy source in KE_v
  type(vertvisc_CS), intent(in) :: CS
    !< Vertical viscosity control structure
  integer, intent(in) :: nii
    !< Zonal block size at meridional velocity points
  integer, intent(in) :: nJJB
    !< Meridional block size at meridional velocity points

  real :: c1(nii,nJJB,SZK_(GV))
    ! Back substitution coefficient [nondim]
  real :: b1(nii,nJJB)
    ! Inverse denominator [H-1 ~> m-1 or m2 kg-1]
  real :: d1(nii,nJJB)
    ! Elimination factor [nondim]
  real :: b_denom_1
    ! First term in the denominator of b1 [H ~> m or kg m-2]
  logical :: do_ij(nii,nJJB)
    ! Wet-point mask
  integer :: i, J, k, isb, ieb, JsbB, JebB, ii, JJ, iie, JJe, nz
    ! Global and block-local indices and bounds

  nz = GV%ke
  !$omp target enter data map(alloc: c1, b1, d1, do_ij)

  do JsbB=G%JscB,G%JecB,nJJB ; do isb=G%isc,G%iec,nii
    ieb = min(isb + nii - 1, G%iec)
    JebB = min(JsbB + nJJB - 1, G%JecB)
    iie = ieb - isb + 1
    JJe = JebB - JsbB + 1
    do concurrent (JJ=1:JJe, ii=1:iie) DO_LOCALITY(local(i,J,b_denom_1))
      i = isb + ii - 1
      J = JsbB + JJ - 1

      do_ij(ii,JJ) = (G%mask2dCv(i,J) > 0.)
      if (do_ij(ii,JJ)) then
        b_denom_1 = CS%h_v(i,J,1) ! CS%a_v_gl90(i,J,1) is zero.
        b1(ii,JJ) = 1.0 / (b_denom_1 + dt*CS%a_v_gl90(i,J,2))
        d1(ii,JJ) = b_denom_1 * b1(ii,JJ)
        ADp%dv_dt_visc_gl90(i,J,1) = b1(ii,JJ) * (CS%h_v(i,J,1) * ADp%dv_dt_visc_gl90(i,J,1))
      endif
    enddo
    do concurrent (JJ=1:JJe) DO_LOCALITY(local(k))
      do k=2,nz
        do concurrent (ii=1:iie) DO_LOCALITY(local(i,J,b_denom_1))
          i = isb + ii - 1
          J = JsbB + JJ - 1

          if (do_ij(ii,JJ)) then
            c1(ii,JJ,k) = dt * CS%a_v_gl90(i,J,K) * b1(ii,JJ)
            b_denom_1 = CS%h_v(i,J,k) + dt * (CS%a_v_gl90(i,J,K) * d1(ii,JJ))
            b1(ii,JJ) = 1.0 / (b_denom_1 + dt * CS%a_v_gl90(i,J,K+1))
            d1(ii,JJ) = b_denom_1 * b1(ii,JJ)
            ADp%dv_dt_visc_gl90(i,J,k) = (CS%h_v(i,J,k) * ADp%dv_dt_visc_gl90(i,J,k) &
                + dt * CS%a_v_gl90(i,J,K) * ADp%dv_dt_visc_gl90(i,J,k-1)) * b1(ii,JJ)
          endif
        enddo
      enddo
      ! Back substitute to solve for new velocities, held by ADp%dv_dt_visc_gl90.
      do k=nz-1,1,-1
        do concurrent (ii=1:iie) DO_LOCALITY(local(i,J))
          i = isb + ii - 1
          J = JsbB + JJ - 1

          if (do_ij(ii,JJ)) ADp%dv_dt_visc_gl90(i,J,k) = &
              ADp%dv_dt_visc_gl90(i,J,k) + c1(ii,JJ,k+1) * ADp%dv_dt_visc_gl90(i,J,k+1)
        enddo
      enddo
    enddo ! JJ-loop over rows within a block

    ! Now fill ADp%dv_dt_visc_gl90(i,J,k) with the actual velocity tendency due to GL90.
    ! On the RHS, ADp%dv_dt_visc(i,J,k) holds the original velocity v(i,J,k),
    ! and ADp%dv_dt_visc_gl90(i,J,k) holds the updated velocity due to GL90.
    do concurrent (k=1:nz, JJ=1:JJe, ii=1:iie) DO_LOCALITY(local(i,J))
      i = isb + ii - 1
      J = JsbB + JJ - 1

      if (do_ij(ii,JJ)) then
        ADp%dv_dt_visc_gl90(i,J,k) = (ADp%dv_dt_visc_gl90(i,J,k) - ADp%dv_dt_visc(i,J,k)) * Idt
        if (abs(ADp%dv_dt_visc_gl90(i,J,k)) < accel_underflow) ADp%dv_dt_visc_gl90(i,J,k) = 0.0
      endif
    enddo
    if (do_KE) then
      ! To compute energetics, multiply by v*h, where v is the velocity before the update.
      ! ADp%dv_dt_visc(i,J,k) still holds that original velocity v(i,J,k).
      do concurrent (k=1:nz, JJ=1:JJe, ii=1:iie) DO_LOCALITY(local(i,J))
        i = isb + ii - 1
        J = JsbB + JJ - 1

        if (do_ij(ii,JJ)) &
          KE_v(i,J,k) = ADp%dv_dt_visc(i,J,k) * CS%h_v(i,J,k) * G%areaCv(i,J) * ADp%dv_dt_visc_gl90(i,J,k)
      enddo
    endif
  enddo ; enddo ! Horizontal blocks
  !$omp target exit data map(release: c1, b1, d1, do_ij)
end subroutine vertvisc_gl90_v


!> Calculate the fraction of momentum originally in a layer that remains in the water column
!! after a time-step of viscosity, equivalently the fraction of a time-step's worth of
!! barotropic acceleration that a layer experiences after viscosity is applied.
subroutine vertvisc_remnant(visc, visc_rem_u, visc_rem_v, dt, G, GV, US, CS)
  type(ocean_grid_type), intent(in)   :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV   !< Ocean vertical grid structure
  type(vertvisc_type),   intent(in)   :: visc !< Viscosities and bottom drag
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_u !< Fraction of a time-step's worth of a
                                              !! barotropic acceleration that a layer experiences after
                                              !! viscosity is applied in the zonal direction [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_v !< Fraction of a time-step's worth of a
                                              !! barotropic acceleration that a layer experiences after
                                              !! viscosity is applied in the meridional direction [nondim]
  real,                  intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  type(vertvisc_CS),     intent(in)    :: CS  !< Vertical viscosity control structure

  ! Local variables

  integer :: nIIB, njj, nii, nJJB ! Resolved horizontal block sizes at u- and v-points.

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(remnant): "// &
         "Module must be initialized before it is used.")

  ! Find the zonal viscous remnant using a modification of a standard tridagonal solver.

  nIIB = CS%niblock
  if (nIIB == 0) nIIB = G%IecB - G%IscB + 1
  njj = CS%njblock
  if (njj == 0) njj = G%jec - G%jsc + 1
  call vertvisc_remnant_u(visc_rem_u, visc, dt, G, GV, CS, nIIB, njj)

  ! Now find the meridional viscous remnant using the robust tridiagonal solver.

  nii = CS%niblock
  if (nii == 0) nii = G%iec - G%isc + 1
  nJJB = CS%njblock
  if (nJJB == 0) nJJB = G%JecB - G%JscB + 1
  call vertvisc_remnant_v(visc_rem_v, visc, dt, G, GV, CS, nii, nJJB)

  if (CS%debug) then
    !$omp target update from(visc_rem_u, visc_rem_v)
    call uvchksum("visc_rem_[uv]", visc_rem_u, visc_rem_v, G%HI, haloshift=0, &
                  scalar_pair=.true.)
  endif
end subroutine vertvisc_remnant


!> Find the zonal viscous remnant in horizontal blocks, preserving vertical solve order.
subroutine vertvisc_remnant_u(visc_rem_u, visc, dt, G, GV, CS, nIIB, njj)
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  real, intent(inout) :: visc_rem_u(SZIB_(G),SZJ_(G),SZK_(GV))
    !< Fraction of barotropic acceleration retained after viscosity [nondim]
  type(vertvisc_type), intent(in) :: visc
    !< Viscosities and bottom drag
  real, intent(in) :: dt
    !< Time step [T ~> s]
  type(vertvisc_CS), intent(in) :: CS
    !< Vertical viscosity control structure
  integer, intent(in) :: nIIB
    !< Zonal block size at zonal velocity points
  integer, intent(in) :: njj
    !< Meridional block size at zonal velocity points

  real :: c1(nIIB,njj,SZK_(GV))
    ! Back substitution coefficient [nondim]
  real :: b1(nIIB,njj)
    ! Inverse denominator [H-1 ~> m-1 or m2 kg-1]
  real :: d1(nIIB,njj)
    ! Elimination factor [nondim]
  real :: Ray
    ! Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
  real :: b_denom_1
    ! Denominator term [H ~> m or kg m-2]
  logical :: do_ij(nIIB,njj)
    ! Wet-point mask
  logical :: have_Ray
    ! Availability of Rayleigh drag
  integer :: I, j, k, IsbB, IebB, jsb, jeb, II, jj, IIe, jje, nz
    ! Global and block-local indices and bounds

  nz = GV%ke
  have_Ray = allocated(visc%Ray_u)

  !$omp target enter data map(alloc: c1, b1, d1, do_ij)

  do jsb=G%jsc,G%jec,njj ; do IsbB=G%IscB,G%IecB,nIIB
    IebB = min(IsbB + nIIB - 1, G%IecB)
    jeb = min(jsb + njj - 1, G%jec)
    IIe = IebB - IsbB + 1
    jje = jeb - jsb + 1
    do concurrent (jj=1:jje, II=1:IIe) DO_LOCALITY(local(I,j,Ray,b_denom_1))
      I = IsbB + II - 1
      j = jsb + jj - 1

      do_ij(II,jj) = (G%mask2dCu(I,j) > 0.)
      if (do_ij(II,jj)) then
        Ray = 0.
        if (have_Ray) Ray = visc%Ray_u(I,j,1)
        b_denom_1 = CS%h_u(I,j,1) + dt * (Ray + CS%a_u(I,j,1))
        b1(II,jj) = 1.0 / (b_denom_1 + dt * CS%a_u(I,j,2))
        d1(II,jj) = b_denom_1 * b1(II,jj)
        visc_rem_u(I,j,1) = b1(II,jj) * CS%h_u(I,j,1)
      endif
    enddo
    do concurrent (jj=1:jje) DO_LOCALITY(local(k))
      do k=2,nz
        do concurrent (II=1:IIe) DO_LOCALITY(local(I,j,Ray,b_denom_1))
          I = IsbB + II - 1
          j = jsb + jj - 1

          if (do_ij(II,jj)) then
            Ray = 0.
            if (have_Ray) Ray = visc%Ray_u(I,j,k)
            c1(II,jj,k) = dt * CS%a_u(I,j,K) * b1(II,jj)
            b_denom_1 = CS%h_u(I,j,k) + dt * (Ray + CS%a_u(I,j,K) * d1(II,jj))
            b1(II,jj) = 1.0 / (b_denom_1 + dt * CS%a_u(I,j,K+1))
            d1(II,jj) = b_denom_1 * b1(II,jj)
            visc_rem_u(I,j,k) = (CS%h_u(I,j,k) + dt * CS%a_u(I,j,K) * visc_rem_u(I,j,k-1)) * b1(II,jj)
          endif
        enddo
        !### Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
        if (dt < 0) exit
      enddo
      do k=nz-1,1,-1
        do concurrent (II=1:IIe) DO_LOCALITY(local(I,j))
          I = IsbB + II - 1
          j = jsb + jj - 1

          if (do_ij(II,jj)) visc_rem_u(I,j,k) = visc_rem_u(I,j,k) + c1(II,jj,k+1) * visc_rem_u(I,j,k+1)
        enddo
      enddo
    enddo ! jj-loop over rows within a block
  enddo ; enddo ! Horizontal blocks
  !$omp target exit data map(release: c1, b1, d1, do_ij)
end subroutine vertvisc_remnant_u


!> Find the meridional viscous remnant in horizontal blocks, preserving vertical solve order.
subroutine vertvisc_remnant_v(visc_rem_v, visc, dt, G, GV, CS, nii, nJJB)
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  real, intent(inout) :: visc_rem_v(SZI_(G),SZJB_(G),SZK_(GV))
    !< Fraction of barotropic acceleration retained after viscosity [nondim]
  type(vertvisc_type), intent(in) :: visc
    !< Viscosities and bottom drag
  real, intent(in) :: dt
    !< Time step [T ~> s]
  type(vertvisc_CS), intent(in) :: CS
    !< Vertical viscosity control structure
  integer, intent(in) :: nii
    !< Zonal block size at meridional velocity points
  integer, intent(in) :: nJJB
    !< Meridional block size at meridional velocity points

  real :: c1(nii,nJJB,SZK_(GV))
    ! Back substitution coefficient [nondim]
  real :: b1(nii,nJJB)
    ! Inverse denominator [H-1 ~> m-1 or m2 kg-1]
  real :: d1(nii,nJJB)
    ! Elimination factor [nondim]
  real :: Ray
    ! Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
  real :: b_denom_1
    ! Denominator term [H ~> m or kg m-2]
  logical :: do_ij(nii,nJJB)
    ! Wet-point mask
  logical :: have_Ray
    ! Availability of Rayleigh drag
  integer :: i, J, k, isb, ieb, JsbB, JebB, ii, JJ, iie, JJe, nz
    ! Global and block-local indices and bounds

  nz = GV%ke
  have_Ray = allocated(visc%Ray_v)
  !$omp target enter data map(alloc: c1, b1, d1, do_ij)

  do JsbB=G%JscB,G%JecB,nJJB ; do isb=G%isc,G%iec,nii
    ieb = min(isb + nii - 1, G%iec)
    JebB = min(JsbB + nJJB - 1, G%JecB)
    iie = ieb - isb + 1
    JJe = JebB - JsbB + 1
    do concurrent (JJ=1:JJe, ii=1:iie) DO_LOCALITY(local(i,J,Ray,b_denom_1))
      i = isb + ii - 1
      J = JsbB + JJ - 1

      do_ij(ii,JJ) = (G%mask2dCv(i,J) > 0.)
      if (do_ij(ii,JJ)) then
        Ray = 0.
        if (have_Ray) Ray = visc%Ray_v(i,J,1)
        b_denom_1 = CS%h_v(i,J,1) + dt * (Ray + CS%a_v(i,J,1))
        b1(ii,JJ) = 1.0 / (b_denom_1 + dt*CS%a_v(i,J,2))
        d1(ii,JJ) = b_denom_1 * b1(ii,JJ)
        visc_rem_v(i,J,1) = b1(ii,JJ) * CS%h_v(i,J,1)
      endif
    enddo
    do concurrent (JJ=1:JJe) DO_LOCALITY(local(k))
      do k=2,nz
        do concurrent (ii=1:iie) DO_LOCALITY(local(i,J,Ray,b_denom_1))
          i = isb + ii - 1
          J = JsbB + JJ - 1

          if (do_ij(ii,JJ)) then
            Ray = 0.
            if (have_Ray) Ray = visc%Ray_v(i,J,k)
            c1(ii,JJ,k) = dt * CS%a_v(i,J,K) * b1(ii,JJ)
            b_denom_1 = CS%h_v(i,J,k) + dt * (Ray + CS%a_v(i,J,K) * d1(ii,JJ))
            b1(ii,JJ) = 1.0 / (b_denom_1 + dt * CS%a_v(i,J,K+1))
            d1(ii,JJ) = b_denom_1 * b1(ii,JJ)
            visc_rem_v(i,J,k) = (CS%h_v(i,J,k) + dt * CS%a_v(i,J,K) * visc_rem_v(i,J,k-1)) * b1(ii,JJ)
          endif
        enddo
        !### Force FMA evaluation of b1 by blocking lookahead with an impossible branch.
        if (dt < 0) exit
      enddo
      do k=nz-1,1,-1
        do concurrent (ii=1:iie) DO_LOCALITY(local(i,J))
          i = isb + ii - 1
          J = JsbB + JJ - 1

          if (do_ij(ii,JJ)) visc_rem_v(i,J,k) = visc_rem_v(i,J,k) + c1(ii,JJ,k+1) * visc_rem_v(i,J,k+1)
        enddo
      enddo
    enddo ! JJ-loop over rows within a block
  enddo ; enddo ! Horizontal blocks
  !$omp target exit data map(release: c1, b1, d1, do_ij)
end subroutine vertvisc_remnant_v


!> Calculate the coupling coefficients (CS%a_u, CS%a_v, CS%a_u_gl90, CS%a_v_gl90)
!! and effective layer thicknesses (CS%h_u and CS%h_v) for later use in the
!! applying the implicit vertical viscosity via vertvisc().
subroutine vertvisc_coef(u, v, h, dz, forces, visc, tv, dt, G, GV, US, CS, OBC, VarMix)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz     !< Vertical distance across layers [Z ~> m]
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any available
                                                   !! thermodynamic fields.
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       intent(inout) :: CS     !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary condition structure
  type(VarMix_CS),         intent(in) :: VarMix !< Variable mixing coefficients
  ! Field from forces used in this subroutine:
  !   ustar: the friction velocity [Z T-1 ~> m s-1], used here as the mixing
  !     velocity in the mixed layer if NKML > 1 in a bulk mixed layer.

  ! Local variables

  real, dimension(SZI_(G),SZJ_(G)) :: &
    Ustar_2d    ! The wind friction velocity, calculated using the Boussinesq reference density or
                ! the time-evolving surface density in non-Boussinesq mode [Z T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: hML_u ! Diagnostic of the mixed layer depth at u points [Z ~> m].
  real, allocatable, dimension(:,:) :: hML_v ! Diagnostic of the mixed layer depth at v points [Z ~> m].
  real, allocatable, dimension(:,:,:) :: Kv_u ! Total vertical viscosity at u-points in
                                              ! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, allocatable, dimension(:,:,:) :: Kv_v ! Total vertical viscosity at v-points in
                                              ! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, allocatable, dimension(:,:,:) :: Kv_gl90_u ! GL90 vertical viscosity at u-points in
                                              ! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, allocatable, dimension(:,:,:) :: Kv_gl90_v ! GL90 vertical viscosity at v-points in
                                              ! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  integer :: nIIB, njj, nii, nJJB ! Resolved horizontal block sizes at u- and v-points.
  integer :: nz ! Number of vertical layers.

  nz = GV%ke

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(coef): "// &
         "Module must be initialized before it is used.")

  if (CS%id_Kv_u > 0) allocate(Kv_u(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke), source=0.0)

  if (CS%id_Kv_v > 0) allocate(Kv_v(G%isd:G%ied,G%JsdB:G%JedB,GV%ke), source=0.0)

  if (CS%id_Kv_gl90_u > 0) allocate(Kv_gl90_u(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke), source=0.0)

  if (CS%id_Kv_gl90_v > 0) allocate(Kv_gl90_v(G%isd:G%ied,G%JsdB:G%JedB,GV%ke), source=0.0)

  if (CS%debug .or. (CS%id_hML_u > 0)) allocate(hML_u(G%IsdB:G%IedB,G%jsd:G%jed), source=0.0)
  if (CS%debug .or. (CS%id_hML_v > 0)) allocate(hML_v(G%isd:G%ied,G%JsdB:G%JedB), source=0.0)

  if ((allocated(visc%taux_shelf) .or. associated(forces%frac_shelf_u)) .and. &
      .not.associated(CS%a1_shelf_u)) then
    allocate(CS%a1_shelf_u(G%IsdB:G%IedB,G%jsd:G%jed), source=0.0)
    !$omp target enter data map(to: CS%a1_shelf_u)
  endif
  if ((allocated(visc%tauy_shelf) .or. associated(forces%frac_shelf_v)) .and. &
      .not.associated(CS%a1_shelf_v)) then
    allocate(CS%a1_shelf_v(G%isd:G%ied,G%JsdB:G%JedB), source=0.0)
    !$omp target enter data map(to: CS%a1_shelf_v)
  endif

  call find_ustar(forces, tv, Ustar_2d, G, GV, US, halo=1)

  ! find_ustar and the callers remain on the host; upload their current inputs once per call.
  !$omp target data map(to: u, v, h, dz, Ustar_2d, forces, tv, VarMix)
  if (associated(forces%frac_shelf_u)) then
    !$omp target enter data map(to: forces%frac_shelf_u)
  endif
  if (associated(forces%frac_shelf_v)) then
    !$omp target enter data map(to: forces%frac_shelf_v)
  endif
  if (allocated(tv%SpV_avg)) then
    !$omp target enter data map(to: tv%SpV_avg)
  endif
  ! Shear fields still have host producers in diabatic mixing, ALE and halo passes.
  ! Their storage is persistent; refresh values without creating a new mapping.
  if (associated(visc%Kv_shear)) then
    !$omp target update to(visc%Kv_shear)
  endif
  if (associated(visc%Kv_shear_Bu)) then
    !$omp target update to(visc%Kv_shear_Bu)
  endif
  if (allocated(VarMix%kdgl90_struct)) then
    !$omp target enter data map(to: VarMix%kdgl90_struct)
  endif
  if (associated(OBC)) then
    !$omp target enter data map(to: OBC)
    !$omp target enter data map(to: OBC%segnum_u, OBC%segnum_v)
  endif
  if (allocated(hML_u)) then
    !$omp target enter data map(to: hML_u)
  endif
  if (allocated(hML_v)) then
    !$omp target enter data map(to: hML_v)
  endif
  if (allocated(Kv_u)) then
    !$omp target enter data map(to: Kv_u)
  endif
  if (allocated(Kv_v)) then
    !$omp target enter data map(to: Kv_v)
  endif
  if (allocated(Kv_gl90_u)) then
    !$omp target enter data map(to: Kv_gl90_u)
  endif
  if (allocated(Kv_gl90_v)) then
    !$omp target enter data map(to: Kv_gl90_v)
  endif

  nIIB = CS%niblock
  if (nIIB == 0) nIIB = G%IecB - G%IscB + 1
  njj = CS%njblock
  if (njj == 0) njj = G%jec - G%jsc + 1
  call vertvisc_coef_u(u, h, dz, forces, visc, tv, dt, G, GV, US, CS, OBC, VarMix, Ustar_2d, &
                       hML_u, Kv_u, Kv_gl90_u, nIIB, njj)
  nii = CS%niblock
  if (nii == 0) nii = G%iec - G%isc + 1
  nJJB = CS%njblock
  if (nJJB == 0) nJJB = G%JecB - G%JscB + 1
  call vertvisc_coef_v(v, h, dz, forces, visc, tv, dt, G, GV, US, CS, OBC, VarMix, Ustar_2d, &
                       hML_v, Kv_v, Kv_gl90_v, nii, nJJB)

  if (allocated(hML_u)) then
    !$omp target exit data map(from: hML_u)
  endif
  if (allocated(hML_v)) then
    !$omp target exit data map(from: hML_v)
  endif
  if (allocated(Kv_u)) then
    !$omp target exit data map(from: Kv_u)
  endif
  if (allocated(Kv_v)) then
    !$omp target exit data map(from: Kv_v)
  endif
  if (allocated(Kv_gl90_u)) then
    !$omp target exit data map(from: Kv_gl90_u)
  endif
  if (allocated(Kv_gl90_v)) then
    !$omp target exit data map(from: Kv_gl90_v)
  endif
  if (associated(OBC)) then
    !$omp target exit data map(release: OBC%segnum_u, OBC%segnum_v)
    !$omp target exit data map(release: OBC)
  endif
  if (allocated(VarMix%kdgl90_struct)) then
    !$omp target exit data map(release: VarMix%kdgl90_struct)
  endif
  if (allocated(tv%SpV_avg)) then
    !$omp target exit data map(release: tv%SpV_avg)
  endif
  if (associated(forces%frac_shelf_v)) then
    !$omp target exit data map(release: forces%frac_shelf_v)
  endif
  if (associated(forces%frac_shelf_u)) then
    !$omp target exit data map(release: forces%frac_shelf_u)
  endif
  !$omp end target data

  ! Host debug/diagnostic consumers need only the requested coefficients.
  if (CS%debug .or. (CS%id_h_u > 0)) then
    !$omp target update from(CS%h_u)
  endif
  if (CS%debug .or. (CS%id_h_v > 0)) then
    !$omp target update from(CS%h_v)
  endif
  if (CS%debug .or. (CS%id_au_vv > 0)) then
    !$omp target update from(CS%a_u)
  endif
  if (CS%debug .or. (CS%id_av_vv > 0)) then
    !$omp target update from(CS%a_v)
  endif
  if (CS%id_au_gl90_vv > 0) then
    !$omp target update from(CS%a_u_gl90)
  endif
  if (CS%id_av_gl90_vv > 0) then
    !$omp target update from(CS%a_v_gl90)
  endif

  if (CS%debug) then
    call uvchksum("vertvisc_coef h_[uv]", CS%h_u, CS%h_v, G%HI, haloshift=0, &
                  unscale=GV%H_to_m, scalar_pair=.true.)
    call uvchksum("vertvisc_coef a_[uv]", CS%a_u, CS%a_v, G%HI, haloshift=0, &
                  unscale=GV%H_to_m*US%s_to_T, scalar_pair=.true.)
    if (allocated(hML_u) .and. allocated(hML_v)) &
      call uvchksum("vertvisc_coef hML_[uv]", hML_u, hML_v, G%HI, &
                    haloshift=0, unscale=US%Z_to_m, scalar_pair=.true.)
  endif
  if (query_averaging_enabled(CS%diag)) then
    if (associated(visc%Kv_slow) .and. (CS%id_Kv_slow > 0)) &
        call post_data(CS%id_Kv_slow, visc%Kv_slow, CS%diag)
    if (CS%id_Kv_u > 0) call post_data(CS%id_Kv_u, Kv_u, CS%diag)
    if (CS%id_Kv_v > 0) call post_data(CS%id_Kv_v, Kv_v, CS%diag)
    if (CS%id_Kv_gl90_u > 0) call post_data(CS%id_Kv_gl90_u, Kv_gl90_u, CS%diag)
    if (CS%id_Kv_gl90_v > 0) call post_data(CS%id_Kv_gl90_v, Kv_gl90_v, CS%diag)
    if (CS%id_au_vv > 0) call post_data(CS%id_au_vv, CS%a_u, CS%diag)
    if (CS%id_av_vv > 0) call post_data(CS%id_av_vv, CS%a_v, CS%diag)
    if (CS%id_au_gl90_vv > 0) call post_data(CS%id_au_gl90_vv, CS%a_u_gl90, CS%diag)
    if (CS%id_av_gl90_vv > 0) call post_data(CS%id_av_gl90_vv, CS%a_v_gl90, CS%diag)
    if (CS%id_h_u > 0) call post_data(CS%id_h_u, CS%h_u, CS%diag)
    if (CS%id_h_v > 0) call post_data(CS%id_h_v, CS%h_v, CS%diag)
    if (CS%id_hML_u > 0) call post_data(CS%id_hML_u, hML_u, CS%diag)
    if (CS%id_hML_v > 0) call post_data(CS%id_hML_v, hML_v, CS%diag)
  endif
  if (allocated(hML_u)) deallocate(hML_u)
  if (allocated(hML_v)) deallocate(hML_v)
end subroutine vertvisc_coef


!> Construct coupling coefficients and effective thicknesses in blocks at zonal velocity points.
subroutine vertvisc_coef_u(u, h, dz, forces, visc, tv, dt, G, GV, US, CS, OBC, VarMix, Ustar_2d, &
                           hML_u, Kv_u, Kv_gl90_u, nIIB, njj)
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  type(unit_scale_type), intent(in) :: US
    !< Unit scaling
  real, intent(in) :: u(SZIB_(G),SZJ_(G),SZK_(GV))
    !< Zonal velocity [L T-1 ~> m s-1]
  real, intent(in) :: h(SZI_(G),SZJ_(G),SZK_(GV))
    !< Layer thickness [H ~> m or kg m-2]
  real, intent(in) :: dz(SZI_(G),SZJ_(G),SZK_(GV))
    !< Vertical layer distance [Z ~> m]
  type(mech_forcing), intent(in) :: forces
    !< Mechanical forcing
  type(vertvisc_type), intent(in) :: visc
    !< Viscosities and bottom drag
  type(thermo_var_ptrs), intent(in) :: tv
    !< Thermodynamic fields
  real, intent(in) :: dt
    !< Time step [T ~> s]
  type(vertvisc_CS), intent(inout) :: CS
    !< Vertical viscosity control structure
  type(ocean_OBC_type), pointer :: OBC
    !< Open boundary conditions
  type(VarMix_CS), intent(in) :: VarMix
    !< Variable mixing coefficients
  real, intent(in) :: Ustar_2d(SZI_(G),SZJ_(G))
    !< Wind friction velocity [Z T-1 ~> m s-1]
  real, allocatable, intent(inout) :: hML_u(:,:)
    !< Mixed layer depth diagnostic [Z ~> m]
  real, allocatable, intent(inout) :: Kv_u(:,:,:)
    !< Total vertical viscosity diagnostic [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  real, allocatable, intent(inout) :: Kv_gl90_u(:,:,:)
    !< GL90 vertical viscosity diagnostic [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  integer, intent(in) :: nIIB
    !< Zonal block size at zonal velocity points
  integer, intent(in) :: njj
    !< Meridional block size at zonal velocity points

  real :: hvel(nIIB,njj,SZK_(GV))
    ! Effective open-ocean layer thickness [H ~> m or kg m-2]
  real :: hvel_shelf(nIIB,njj,SZK_(GV))
    ! Effective layer thickness under ice shelves [H ~> m or kg m-2]
  real :: dz_harm(nIIB,njj,SZK_(GV))
    ! Harmonic mean of adjacent layer distances [Z ~> m]
  real :: dz_vel(nIIB,njj,SZK_(GV))
    ! Effective open-ocean layer distance [Z ~> m]
  real :: dz_vel_shelf(nIIB,njj,SZK_(GV))
    ! Effective layer distance under ice shelves [Z ~> m]
  real :: z_i(nIIB,njj,SZK_(GV)+1)
    ! Interface height above the bottom, normalized by bottom boundary layer thickness [nondim]
  real :: z_i_gl90(nIIB,njj,SZK_(GV)+1)
    ! Interface height above the bottom, normalized by GL90 bottom boundary layer thickness [nondim]
  real :: a_cpl(nIIB,njj,SZK_(GV)+1)
    ! Open-ocean coupling coefficient [H T-1 ~> m s-1 or Pa s m-1]
  real :: a_cpl_gl90(nIIB,njj,SZK_(GV)+1)
    ! GL90 coupling coefficient, part of a_cpl [H T-1 ~> m s-1 or Pa s m-1]
  real :: a_shelf(nIIB,njj,SZK_(GV)+1)
    ! Coupling coefficient under ice shelves [H T-1 ~> m s-1 or Pa s m-1]
  real :: kv_bbl(nIIB,njj)
    ! Bottom boundary layer viscosity [H Z T-1 ~> m2 s-1 or Pa s]
  real :: bbl_thick(nIIB,njj)
    ! Bottom boundary layer thickness [Z ~> m]
  real :: I_Hbbl(nIIB,njj)
    ! Inverse bottom boundary layer thickness [Z-1 ~> m-1]
  real :: I_HTbl(nIIB,njj)
    ! Inverse top boundary layer thickness [Z-1 ~> m-1]
  real :: Ztop_min(nIIB,njj)
    ! Deeper of the two adjacent surface heights [Z ~> m]
  real :: Dmin(nIIB,njj)
    ! Shallower of the two adjacent bottom depths [Z ~> m]
  real :: zh(nIIB,njj)
    ! Distance from the bottom based on accumulated harmonic mean layer distances [Z ~> m]
  real :: zcol(nIIB,njj)
    ! Interface height at the adjacent western tracer point [Z ~> m]
  real :: zcol_p1(nIIB,njj)
    ! Interface height at the adjacent eastern tracer point [Z ~> m]
  real :: h_ml(nIIB,njj)
    ! Mixed layer depth [Z ~> m]
  real :: h_harm
    ! Harmonic mean of adjacent layer thicknesses [H ~> m or kg m-2]
  real :: h_arith
    ! Arithmetic mean of adjacent layer thicknesses [H ~> m or kg m-2]
  real :: h_delta
    ! Lateral difference of adjacent layer thicknesses [H ~> m or kg m-2]
  real :: dz_arith
    ! Arithmetic mean of adjacent layer distances [Z ~> m]
  real :: I_Hbbl_gl90
    ! Inverse GL90 bottom boundary layer thickness [Z-1 ~> m-1]
  real :: I_valBL
    ! Inverse scaling factor for the harmonic boundary layer criterion [nondim]
  real :: botfn
    ! Weight equal to 1 at the bottom and tending to 0 in the interior [nondim]
  real :: topfn
    ! Weight equal to 1 at the top and tending to 0 in the interior [nondim]
  real :: z2
    ! Distance from the boundary normalized by its boundary layer thickness [nondim]
  real :: z2_wt
    ! Transition weight used to calculate z2 [nondim]
  real :: z_clear
    ! Clearance of an interface above the surrounding topography [Z ~> m]
  real :: a_cpl_max
    ! Maximum representable coupling coefficient [H T-1 ~> m s-1 or Pa s m-1]
  real :: h_neglect
    ! Negligible layer thickness, usually lost in roundoff [H ~> m or kg m-2]
  real :: dz_neglect
    ! Negligible vertical distance, usually lost in roundoff [Z ~> m]
  logical :: do_ij(nIIB,njj)
    ! Wet-point mask for the active block
  logical :: do_ij_shelf(nIIB,njj)
    ! Mask of active block points with ice shelves
  logical :: do_any_shelf
    ! If true, the block contains active ice-shelf points
  integer :: zi_dir(nIIB,njj)
    ! Ternary selector for the adjacent column used to calculate z_clear
  integer :: i, j, k, js, je, Isq, Ieq, nz
  integer :: IsbB, IebB, jsb, jeb, II, jj, IIe, jje

  js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB
  nz = GV%ke

  h_neglect = GV%H_subroundoff ; dz_neglect = GV%dZ_subroundoff
  a_cpl_max = 1.0e37 * GV%m_to_H * US%T_to_s
  I_valBL = 0.0 ; if (CS%harm_BL_val > 0.0) I_valBL = 1.0 / CS%harm_BL_val
  I_Hbbl_gl90 = 0.0
  if (CS%use_GL90_in_SSW) I_Hbbl_gl90 = 1. / (CS%Hbbl_gl90 + dz_neglect)

  !$omp target enter data &
  !$omp   map(alloc: hvel, dz_harm, dz_vel, hvel_shelf, dz_vel_shelf, z_i, z_i_gl90, a_cpl, a_cpl_gl90, &
  !$omp     a_shelf, kv_bbl, bbl_thick, I_Hbbl, I_HTbl, Ztop_min, Dmin, zh, zcol, zcol_p1, h_ml, zi_dir, &
  !$omp     do_ij, do_ij_shelf)

  do jsb=js,je,njj ; do IsbB=Isq,Ieq,nIIB
    IebB = min(IsbB + nIIB - 1, Ieq)
    jeb = min(jsb + njj - 1, je)
    IIe = IebB - IsbB + 1
    jje = jeb - jsb + 1
    do concurrent (jj=1:jje, II=1:IIe) DO_LOCALITY(local(I,j))
      i = IsbB + II - 1
      j = jsb + jj - 1

      do_ij(II,jj) = G%mask2dCu(I,j) > 0.
      do_ij_shelf(II,jj) = .false.
      if (do_ij(II,jj)) then
        I_Hbbl(II,jj) = 1. / (CS%Hbbl + dz_neglect)

        if (CS%bottomdraglaw) then
          kv_bbl(II,jj) = visc%Kv_bbl_u(I,j)
          bbl_thick(II,jj) = visc%bbl_thick_u(I,j) + dz_neglect
          I_Hbbl(II,jj) = 1. / bbl_thick(II,jj)
        endif

        Dmin(II,jj) = min(G%bathyT(i,j), G%bathyT(i+1,j))
        zi_dir(II,jj) = 0

        ! Project thickness outward across OBCs using a zero-gradient condition.
        if (associated(OBC)) then
          if (OBC%u_E_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) > 0) then
              Dmin(II,jj) = G%bathyT(i,j)
              zi_dir(II,jj) = -1
            endif
          endif

          if (OBC%u_W_OBCs_on_PE) then
            if (OBC%segnum_u(I,j) < 0) then
              Dmin(II,jj) = G%bathyT(i+1,j)
              zi_dir(II,jj) = 1
            endif
          endif
        endif

        ! The following block calculates the thicknesses at velocity grid points for
        ! the vertical viscosity (hvel and dz_vel).  Near the bottom an upwind biased
        ! thickness is used to control the effect of spurious Montgomery potential
        ! gradients at the bottom where nearly massless layers layers ride over the
        ! topography.

        z_i(II,jj,nz+1) = 0.

        if (.not. CS%harmonic_visc) then
          zh(II,jj) = 0.
          zcol(II,jj) = -G%bathyT(i,j)
          zcol_p1(II,jj) = -G%bathyT(i+1,j)
        endif

        if (CS%use_GL90_in_SSW) then
          z_i_gl90(II,jj,nz+1) = 0.
        endif
      endif
    enddo

    do concurrent (jj=1:jje) DO_LOCALITY(local(k))
      do k=nz,1,-1
        do concurrent (II=1:IIe, do_ij(II,jj)) &
            DO_LOCALITY(local(I,j,h_harm,h_arith,h_delta,dz_arith,z2,botfn,z_clear,z2_wt))
          i = IsbB + II - 1
          j = jsb + jj - 1

          h_harm = 2. * h(i,j,k) * h(i+1,j,k) / (h(i,j,k) + h(i+1,j,k) + h_neglect)
          h_arith = 0.5 * (h(i+1,j,k) + h(i,j,k))
          h_delta = h(i+1,j,k) - h(i,j,k)
          dz_harm(II,jj,k) = 2. * dz(i,j,k) * dz(i+1,j,k) / (dz(i,j,k) + dz(i+1,j,k) + dz_neglect)
          dz_arith = 0.5 * (dz(i+1,j,k) + dz(i,j,k))

          ! Project thickness outward across OBCs using a zero-gradient condition.
          if (associated(OBC)) then
            if (OBC%u_E_OBCs_on_PE) then
              if (OBC%segnum_u(I,j) > 0) then
                h_harm = h(i,j,k)
                h_arith = h(i,j,k)
                h_delta = 0.
                dz_harm(II,jj,k) = dz(i,j,k)
                dz_arith = dz(i,j,k)
              endif
            endif

            if (OBC%u_W_OBCs_on_PE) then
              if (OBC%segnum_u(I,j) < 0) then
                h_harm = h(i+1,j,k)
                h_arith = h(i+1,j,k)
                h_delta = 0.
                dz_harm(II,jj,k) = dz(i+1,j,k)
                dz_arith = dz(i+1,j,k)
              endif
            endif
          endif

          if (CS%harmonic_visc) then
            ! The following block calculates the thicknesses at velocity grid points
            ! for the vertical viscosity (hvel and dz_vel).  Near the bottom an
            ! upwind biased thickness is used to control the effect of spurious
            ! Montgomery potential gradients at the bottom where nearly massless
            ! layers ride over the topography.

            hvel(II,jj,k) = h_harm
            dz_vel(II,jj,k) = dz_harm(II,jj,k)

            if (u(I,j,k) * h_delta < 0) then
              z2 = z_i(II,jj,k+1)
              botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

              hvel(II,jj,k) = (1. - botfn) * h_harm + botfn * h_arith
              dz_vel(II,jj,k) = (1. - botfn) * dz_harm(II,jj,k) + botfn * dz_arith
            endif

            z_i(II,jj,k) =  z_i(II,jj,k+1) + dz_harm(II,jj,k) * I_Hbbl(II,jj)
          else
            zcol(II,jj) = zcol(II,jj) + dz(i,j,k)
            zcol_p1(II,jj) = zcol_p1(II,jj) + dz(i+1,j,k)

            zh(II,jj) = zh(II,jj) + dz_harm(II,jj,k)

            z_clear = max(zcol(II,jj), zcol_p1(II,jj)) + Dmin(II,jj)
            if (zi_dir(II,jj) < 0) z_clear = zcol(II,jj) + Dmin(II,jj)
            if (zi_dir(II,jj) > 0) z_clear = zcol_p1(II,jj) + Dmin(II,jj)

            z_i(II,jj,k) = max(zh(II,jj), z_clear) * I_Hbbl(II,jj)

            hvel(II,jj,k) = h_arith
            dz_vel(II,jj,k) = dz_arith

            if (u(I,j,k) * h_delta > 0.) then
              if (zh(II,jj) * I_Hbbl(II,jj) < CS%harm_BL_val) then
                hvel(II,jj,k) = h_harm
                dz_vel(II,jj,k) = dz_harm(II,jj,k)
              else
                z2_wt = 1.
                if (zh(II,jj) * I_Hbbl(II,jj) < 2. * CS%harm_BL_val) &
                  z2_wt = max(0., min(1., zh(II,jj) * I_Hbbl(II,jj) * I_valBL - 1.))

                z2 = z2_wt * (max(zh(II,jj), z_clear) * I_Hbbl(II,jj))
                botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

                hvel(II,jj,k) = (1. - botfn) * h_arith + botfn * h_harm
                dz_vel(II,jj,k) = (1. - botfn) * dz_arith + botfn * dz_harm(II,jj,k)
              endif
            endif
          endif

          if (CS%use_GL90_in_SSW) then
            ! The following block calculates the normalized height above the GL90 BBL
            ! (z_i_gl90), using a harmonic mean between layer thicknesses. For the
            ! GL90 BBL we use simply a constant (Hbbl_gl90). The purpose isthat the
            ! GL90 coupling coefficient is zeroed out within Hbbl_gl90, to ensure
            ! that no momentum gets fluxed into vanished layers. The scheme is not
            ! sensitive to the exact value of Hbbl_gl90, as long as it is in a
            ! reasonable range (~1-20 m): large enough to capture vanished layers
            ! over topography, small enough to not contaminate the interior.

            z_i_gl90(II,jj,k) = z_i_gl90(II,jj,k+1) + dz_harm(II,jj,k) * I_Hbbl_gl90
          endif
        enddo
      enddo
    enddo

    call find_coupling_coef(a_cpl, dz_vel, do_ij, dz_harm, bbl_thick, kv_bbl, z_i, &
        h_ml, dt, G, GV, US, CS, visc, Ustar_2d, tv, .true., OBC, IsbB, IebB, jsb, jeb, nIIB, njj)

    if (allocated(hML_u)) then
      do concurrent (jj=1:jje, II=1:IIe) DO_LOCALITY(local(I,j))
        i = IsbB + II - 1
        j = jsb + jj - 1

        if (do_ij(II,jj)) hML_u(I,j) = h_ml(II,jj)
      enddo
    endif

    if (CS%use_GL90_in_SSW) then
      call find_coupling_coef_gl90(a_cpl_gl90, dz_vel, do_ij, z_i_gl90, G, GV, &
          CS, VarMix, .true., IsbB, IebB, jsb, jeb, nIIB, njj)
    endif

    do_any_shelf = .false.
    if (associated(forces%frac_shelf_u)) then
      do concurrent (jj=1:jje, II=1:IIe, do_ij(II,jj)) DO_LOCALITY(local(I,j) reduce(.or.:do_any_shelf))
        i = IsbB + II - 1
        j = jsb + jj - 1

        CS%a1_shelf_u(I,j) = 0.
        do_ij_shelf(II,jj) = forces%frac_shelf_u(I,j) > 0.
        do_any_shelf = do_any_shelf .or. do_ij_shelf(II,jj)
      enddo
    endif

    if (do_any_shelf) then
      do concurrent (jj=1:jje, II=1:IIe, do_ij_shelf(II,jj)) DO_LOCALITY(local(I,j))
        i = IsbB + II - 1
        j = jsb + jj - 1

        if (.not. CS%harmonic_visc) then
          zh(II,jj) = 0.
          Ztop_min(II,jj) = min(zcol(II,jj), zcol_p1(II,jj))
          I_HTbl(II,jj) = 1. / (visc%tbl_thick_shelf_u(I,j) + dz_neglect)
        endif
      enddo

      do concurrent (jj=1:jje) DO_LOCALITY(local(k))
        do k=1,nz
          do concurrent (II=1:IIe, do_ij_shelf(II,jj)) &
              DO_LOCALITY(local(I,j,h_harm,h_arith,h_delta,dz_arith,z2,z2_wt,topfn))
            i = IsbB + II - 1
            j = jsb + jj - 1

            if (CS%harmonic_visc) then
              hvel_shelf(II,jj,k) = hvel(II,jj,k)
              dz_vel_shelf(II,jj,k) = dz_vel(II,jj,k)
            else
              ! Find upwind-biased thickness near the surface.
              ! (Perhaps this needs to be done more carefully, via find_eta.)

              h_harm = 2. * h(i,j,k) * h(i+1,j,k) &
                  / (h(i,j,k) + h(i+1,j,k) + h_neglect)
              h_arith = 0.5 * (h(i+1,j,k) + h(i,j,k))
              h_delta = h(i+1,j,k) - h(i,j,k)
              dz_arith = 0.5 * (dz(i+1,j,k) + dz(i,j,k))

              if (associated(OBC)) then
                if (OBC%u_E_OBCs_on_PE) then
                  if (OBC%segnum_u(I,j) > 0) then
                    h_harm = h(i,j,k)
                    h_arith = h(i,j,k)
                    h_delta = 0.
                    dz_arith = dz(i,j,k)
                  endif
                endif

                if (OBC%u_W_OBCs_on_PE) then
                  if (OBC%segnum_u(I,j) < 0) then
                    h_harm = h(i+1,j,k)
                    h_arith = h(i+1,j,k)
                    h_delta = 0.
                    dz_arith = dz(i+1,j,k)
                  endif
                endif
              endif

              zcol(II,jj) = zcol(II,jj) - dz(i,j,k)
              zcol_p1(II,jj) = zcol_p1(II,jj) - dz(i+1,j,k)

              zh(II,jj) = zh(II,jj) + dz_harm(II,jj,k)

              hvel_shelf(II,jj,k) = hvel(II,jj,k)
              dz_vel_shelf(II,jj,k) = dz_vel(II,jj,k)

              if (u(I,j,k) * h_delta > 0.) then
                if (zh(II,jj) * I_HTbl(II,jj) < CS%harm_BL_val) then
                  hvel_shelf(II,jj,k) = min(hvel(II,jj,k), h_harm)
                  dz_vel_shelf(II,jj,k) = min(dz_vel(II,jj,k), dz_harm(II,jj,k))
                else
                  z2_wt = 1.
                  if (zh(II,jj) * I_HTbl(II,jj) < 2. * CS%harm_BL_val) then
                    z2_wt = max(0., min(1., zh(II,jj) * I_HTbl(II,jj) * I_valBL - 1.))
                  endif

                  z2 = z2_wt * (max(zh(II,jj), Ztop_min(II,jj) - min(zcol(II,jj), zcol_p1(II,jj))) * I_HTbl(II,jj))
                  ! TODO: replace **6 with multiply
                  topfn = 1. / (1. + 0.09 * z2**6)

                  hvel_shelf(II,jj,k) = min(hvel(II,jj,k), (1. - topfn) * h_arith + topfn * h_harm)
                  dz_vel_shelf(II,jj,k) = min(dz_vel(II,jj,k), (1. - topfn) * dz_arith + topfn * dz_harm(II,jj,k))
                endif
              endif
            endif
          enddo
        enddo
      enddo

      call find_coupling_coef(a_shelf, dz_vel_shelf, do_ij_shelf, dz_harm, &
          bbl_thick, kv_bbl, z_i, h_ml, dt, G, GV, US, CS, visc, Ustar_2d, &
          tv, .true., OBC, IsbB, IebB, jsb, jeb, nIIB, njj, shelf=.true.)

      do concurrent (jj=1:jje, II=1:IIe) DO_LOCALITY(local(I,j))
        i = IsbB + II - 1
        j = jsb + jj - 1

        if (do_ij_shelf(II,jj)) CS%a1_shelf_u(I,j) = a_shelf(II,jj,1)
      enddo
    endif

    do concurrent (K=1:nz+1, jj=1:jje, II=1:IIe, do_ij(II,jj)) DO_LOCALITY(local(I,j))
      I = IsbB + II - 1
      j = jsb + jj - 1

      if (do_ij_shelf(II,jj)) then
        ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
        ! CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j) * max(a_shelf(II,jj,K), a_cpl(II,jj,K)) + &
        !                                (1. - forces%frac_shelf_u(I,j)) * a_cpl(II,jj,K))
        if (CS%use_GL90_in_SSW) then
          CS%a_u(I,j,K) = min(a_cpl_max, (forces%frac_shelf_u(I,j) * a_shelf(II,jj,K) + &
              (1. - forces%frac_shelf_u(I,j)) * a_cpl(II,jj,K)) + a_cpl_gl90(II,jj,K))
          CS%a_u_gl90(I,j,K) = min(a_cpl_max, a_cpl_gl90(II,jj,K))
        else
          CS%a_u(I,j,K) = min(a_cpl_max, (forces%frac_shelf_u(I,j) * a_shelf(II,jj,K) + &
              (1. - forces%frac_shelf_u(I,j)) * a_cpl(II,jj,K)))
        endif
      else
        if (CS%use_GL90_in_SSW) then
          a_cpl(II,jj,K) = a_cpl(II,jj,K) + a_cpl_gl90(II,jj,K)
          CS%a_u_gl90(I,j,K) = min(a_cpl_max, a_cpl_gl90(II,jj,K))
        endif

        CS%a_u(I,j,K) = min(a_cpl_max, a_cpl(II,jj,K))
      endif
    enddo
    do concurrent (k=1:nz, jj=1:jje, II=1:IIe, do_ij(II,jj)) DO_LOCALITY(local(I,j))
      i = IsbB + II - 1
      j = jsb + jj - 1

      if (do_ij_shelf(II,jj)) then
        ! Should we instead take the inverse of the average of the inverses?
        CS%h_u(I,j,k) = forces%frac_shelf_u(I,j) * hvel_shelf(II,jj,k) &
            + (1. - forces%frac_shelf_u(I,j)) * hvel(II,jj,k) + h_neglect
      else
        CS%h_u(I,j,k) = hvel(II,jj,k) + h_neglect
      endif
    enddo

    ! Diagnose total Kv at u-points
    if (allocated(Kv_u)) then
      do concurrent (k=1:nz, jj=1:jje, II=1:IIe, do_ij(II,jj)) DO_LOCALITY(local(I,j))
        i = IsbB + II - 1
        j = jsb + jj - 1

        Kv_u(I,j,k) = 0.5 * (CS%a_u(I,j,K) + CS%a_u(I,j,K+1)) * CS%h_u(I,j,k)
      enddo
    endif

    ! Diagnose GL90 Kv at u-points
    if (allocated(Kv_gl90_u)) then
      do concurrent (k=1:nz, jj=1:jje, II=1:IIe, do_ij(II,jj)) DO_LOCALITY(local(I,j))
        i = IsbB + II - 1
        j = jsb + jj - 1

        Kv_gl90_u(I,j,k) = 0.5 * (CS%a_u_gl90(I,j,K) + CS%a_u_gl90(I,j,K+1)) * CS%h_u(I,j,k)
      enddo
    endif
  enddo ; enddo ! Horizontal blocks
  !$omp target exit data &
  !$omp   map(release: hvel, dz_harm, dz_vel, hvel_shelf, dz_vel_shelf, z_i, z_i_gl90, a_cpl, a_cpl_gl90, &
  !$omp     a_shelf, kv_bbl, bbl_thick, I_Hbbl, I_HTbl, Ztop_min, Dmin, zh, zcol, zcol_p1, h_ml, zi_dir, &
  !$omp     do_ij, do_ij_shelf)
end subroutine vertvisc_coef_u


!> Construct coupling coefficients and effective thicknesses in blocks at meridional velocity points.
subroutine vertvisc_coef_v(v, h, dz, forces, visc, tv, dt, G, GV, US, CS, OBC, VarMix, Ustar_2d, &
                           hML_v, Kv_v, Kv_gl90_v, nii, nJJB)
  type(ocean_grid_type), intent(in) :: G
    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV
    !< Vertical grid structure
  type(unit_scale_type), intent(in) :: US
    !< Unit scaling
  real, intent(in) :: v(SZI_(G),SZJB_(G),SZK_(GV))
    !< Meridional velocity [L T-1 ~> m s-1]
  real, intent(in) :: h(SZI_(G),SZJ_(G),SZK_(GV))
    !< Layer thickness [H ~> m or kg m-2]
  real, intent(in) :: dz(SZI_(G),SZJ_(G),SZK_(GV))
    !< Vertical layer distance [Z ~> m]
  type(mech_forcing), intent(in) :: forces
    !< Mechanical forcing
  type(vertvisc_type), intent(in) :: visc
    !< Viscosities and bottom drag
  type(thermo_var_ptrs), intent(in) :: tv
    !< Thermodynamic fields
  real, intent(in) :: dt
    !< Time step [T ~> s]
  type(vertvisc_CS), intent(inout) :: CS
    !< Vertical viscosity control structure
  type(ocean_OBC_type), pointer :: OBC
    !< Open boundary conditions
  type(VarMix_CS), intent(in) :: VarMix
    !< Variable mixing coefficients
  real, intent(in) :: Ustar_2d(SZI_(G),SZJ_(G))
    !< Wind friction velocity [Z T-1 ~> m s-1]
  real, allocatable, intent(inout) :: hML_v(:,:)
    !< Mixed layer depth diagnostic [Z ~> m]
  real, allocatable, intent(inout) :: Kv_v(:,:,:)
    !< Total vertical viscosity diagnostic [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  real, allocatable, intent(inout) :: Kv_gl90_v(:,:,:)
    !< GL90 vertical viscosity diagnostic [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  integer, intent(in) :: nii
    !< Zonal block size at meridional velocity points
  integer, intent(in) :: nJJB
    !< Meridional block size at meridional velocity points

  real :: hvel(nii,nJJB,SZK_(GV))
    ! Effective open-ocean layer thickness [H ~> m or kg m-2]
  real :: hvel_shelf(nii,nJJB,SZK_(GV))
    ! Effective layer thickness under ice shelves [H ~> m or kg m-2]
  real :: dz_harm(nii,nJJB,SZK_(GV))
    ! Harmonic mean of adjacent layer distances [Z ~> m]
  real :: dz_vel(nii,nJJB,SZK_(GV))
    ! Effective open-ocean layer distance [Z ~> m]
  real :: dz_vel_shelf(nii,nJJB,SZK_(GV))
    ! Effective layer distance under ice shelves [Z ~> m]
  real :: z_i(nii,nJJB,SZK_(GV)+1)
    ! Interface height above the bottom, normalized by bottom boundary layer thickness [nondim]
  real :: z_i_gl90(nii,nJJB,SZK_(GV)+1)
    ! Interface height above the bottom, normalized by GL90 bottom boundary layer thickness [nondim]
  real :: a_cpl(nii,nJJB,SZK_(GV)+1)
    ! Open-ocean coupling coefficient [H T-1 ~> m s-1 or Pa s m-1]
  real :: a_cpl_gl90(nii,nJJB,SZK_(GV)+1)
    ! GL90 coupling coefficient, part of a_cpl [H T-1 ~> m s-1 or Pa s m-1]
  real :: a_shelf(nii,nJJB,SZK_(GV)+1)
    ! Coupling coefficient under ice shelves [H T-1 ~> m s-1 or Pa s m-1]
  real :: kv_bbl(nii,nJJB)
    ! Bottom boundary layer viscosity [H Z T-1 ~> m2 s-1 or Pa s]
  real :: bbl_thick(nii,nJJB)
    ! Bottom boundary layer thickness [Z ~> m]
  real :: I_Hbbl(nii,nJJB)
    ! Inverse bottom boundary layer thickness [Z-1 ~> m-1]
  real :: I_HTbl(nii,nJJB)
    ! Inverse top boundary layer thickness [Z-1 ~> m-1]
  real :: Ztop_min(nii,nJJB)
    ! Deeper of the two adjacent surface heights [Z ~> m]
  real :: Dmin(nii,nJJB)
    ! Shallower of the two adjacent bottom depths [Z ~> m]
  real :: zh(nii,nJJB)
    ! Distance from the bottom based on accumulated harmonic mean layer distances [Z ~> m]
  real :: zcol(nii,nJJB)
    ! Interface height at the adjacent southern tracer point [Z ~> m]
  real :: zcol_p1(nii,nJJB)
    ! Interface height at the adjacent northern tracer point [Z ~> m]
  real :: h_ml(nii,nJJB)
    ! Mixed layer depth [Z ~> m]
  real :: h_harm
    ! Harmonic mean of adjacent layer thicknesses [H ~> m or kg m-2]
  real :: h_arith
    ! Arithmetic mean of adjacent layer thicknesses [H ~> m or kg m-2]
  real :: h_delta
    ! Lateral difference of adjacent layer thicknesses [H ~> m or kg m-2]
  real :: dz_arith
    ! Arithmetic mean of adjacent layer distances [Z ~> m]
  real :: I_Hbbl_gl90
    ! Inverse GL90 bottom boundary layer thickness [Z-1 ~> m-1]
  real :: I_valBL
    ! Inverse scaling factor for the harmonic boundary layer criterion [nondim]
  real :: botfn
    ! Weight equal to 1 at the bottom and tending to 0 in the interior [nondim]
  real :: topfn
    ! Weight equal to 1 at the top and tending to 0 in the interior [nondim]
  real :: z2
    ! Distance from the boundary normalized by its boundary layer thickness [nondim]
  real :: z2_wt
    ! Transition weight used to calculate z2 [nondim]
  real :: z_clear
    ! Clearance of an interface above the surrounding topography [Z ~> m]
  real :: a_cpl_max
    ! Maximum representable coupling coefficient [H T-1 ~> m s-1 or Pa s m-1]
  real :: h_neglect
    ! Negligible layer thickness, usually lost in roundoff [H ~> m or kg m-2]
  real :: dz_neglect
    ! Negligible vertical distance, usually lost in roundoff [Z ~> m]
  logical :: do_ij(nii,nJJB)
    ! Wet-point mask for the active block
  logical :: do_ij_shelf(nii,nJJB)
    ! Mask of active block points with ice shelves
  logical :: do_any_shelf
    ! If true, the block contains active ice-shelf points
  integer :: zi_dir(nii,nJJB)
    ! Ternary selector for the adjacent column used to calculate z_clear
  integer :: i, j, k, is, ie, Jsq, Jeq, nz
  integer :: isb, ieb, JsbB, JebB, ii, JJ, iie, JJe

  is = G%isc
  ie = G%iec
  Jsq = G%JscB
  Jeq = G%JecB
  nz = GV%ke
  h_neglect = GV%H_subroundoff ; dz_neglect = GV%dZ_subroundoff
  a_cpl_max = 1.0e37 * GV%m_to_H * US%T_to_s
  I_valBL = 0.0 ; if (CS%harm_BL_val > 0.0) I_valBL = 1.0 / CS%harm_BL_val
  I_Hbbl_gl90 = 0.0
  if (CS%use_GL90_in_SSW) I_Hbbl_gl90 = 1. / (CS%Hbbl_gl90 + dz_neglect)

  !$omp target enter data &
  !$omp   map(alloc: hvel, dz_harm, dz_vel, hvel_shelf, dz_vel_shelf, z_i, z_i_gl90, a_cpl, a_cpl_gl90, &
  !$omp     a_shelf, kv_bbl, bbl_thick, I_Hbbl, I_HTbl, Ztop_min, Dmin, zh, zcol, zcol_p1, h_ml, zi_dir, &
  !$omp     do_ij, do_ij_shelf)

  do JsbB=Jsq,Jeq,nJJB ; do isb=is,ie,nii
    ieb = min(isb + nii - 1, ie)
    JebB = min(JsbB + nJJB - 1, Jeq)
    iie = ieb - isb + 1
    JJe = JebB - JsbB + 1
    do concurrent (JJ=1:JJe, ii=1:iie) DO_LOCALITY(local(i,J))
      i = isb + ii - 1
      j = JsbB + JJ - 1

      do_ij(ii,JJ) = G%mask2dCv(i,J) > 0.
      do_ij_shelf(ii,JJ) = .false.
      if (do_ij(ii,JJ)) then
        I_Hbbl(ii,JJ) = 1. / (CS%Hbbl + dz_neglect)

        if (CS%bottomdraglaw) then
          kv_bbl(ii,JJ) = visc%Kv_bbl_v(i,J)
          bbl_thick(ii,JJ) = visc%bbl_thick_v(i,J) + dz_neglect
          I_Hbbl(ii,JJ) = 1. / bbl_thick(ii,JJ)
        endif

        Dmin(ii,JJ) = min(G%bathyT(i,j), G%bathyT(i,j+1))
        zi_dir(ii,JJ) = 0

        ! Project thickness outward across OBCs using a zero-gradient condition.
        if (associated(OBC)) then
          if (OBC%v_N_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) > 0) then
              Dmin(ii,JJ) = G%bathyT(i,j)
              zi_dir(ii,JJ) = -1
            endif
          endif

          if (OBC%v_S_OBCs_on_PE) then
            if (OBC%segnum_v(i,J) < 0) then
              Dmin(ii,JJ) = G%bathyT(i,j+1)
              zi_dir(ii,JJ) = 1
            endif
          endif
        endif

        z_i(ii,JJ,nz+1) = 0.

        if (.not. CS%harmonic_visc) then
          zh(ii,JJ) = 0.
          zcol(ii,JJ) = -G%bathyT(i,j)
          zcol_p1(ii,JJ) = -G%bathyT(i,j+1)
        endif

        if (CS%use_GL90_in_SSW) then
          z_i_gl90(ii,JJ,nz+1) = 0.
        endif
      endif
    enddo

    do concurrent (JJ=1:JJe) DO_LOCALITY(local(k))
      do k=nz,1,-1
        do concurrent (ii=1:iie, do_ij(ii,JJ)) &
            DO_LOCALITY(local(i,J,h_harm,h_arith,h_delta,dz_arith,z2,botfn,z_clear,z2_wt))
          i = isb + ii - 1
          j = JsbB + JJ - 1

          h_harm = 2. * h(i,j,k) * h(i,j+1,k) / (h(i,j,k) + h(i,j+1,k) + h_neglect)
          h_arith = 0.5 * (h(i,j+1,k) + h(i,j,k))
          h_delta = h(i,j+1,k) - h(i,j,k)
          dz_harm(ii,JJ,k) = 2. * dz(i,j,k) * dz(i,j+1,k) / (dz(i,j,k) + dz(i,j+1,k) + dz_neglect)
          dz_arith = 0.5 * (dz(i,j+1,k) + dz(i,j,k))

          ! Project thickness outward across OBCs using a zero-gradient condition.
          if (associated(OBC)) then
            if (OBC%v_N_OBCs_on_PE) then
              if (OBC%segnum_v(i,J) > 0) then
                h_harm = h(i,j,k)
                h_arith = h(i,j,k)
                h_delta = 0.
                dz_harm(ii,JJ,k) = dz(i,j,k)
                dz_arith = dz(i,j,k)
              endif
            endif

            if (OBC%v_S_OBCs_on_PE) then
              if (OBC%segnum_v(i,J) < 0) then
                h_harm = h(i,j+1,k)
                h_arith = h(i,j+1,k)
                h_delta = 0.
                dz_harm(ii,JJ,k) = dz(i,j+1,k)
                dz_arith = dz(i,j+1,k)
              endif
            endif
          endif

          if (CS%harmonic_visc) then
            ! The following block calculates the thicknesses at velocity grid points
            ! for the vertical viscosity (hvel and dz_vel).  Near the bottom an
            ! upwind biased thickness is used to control the effect of spurious
            ! Montgomery potential gradients at the bottom where nearly massless
            ! layers ride over the topography.

            hvel(ii,JJ,k) = h_harm
            dz_vel(ii,JJ,k) = dz_harm(ii,JJ,k)

            if (v(i,J,k) * h_delta < 0) then
              z2 = z_i(ii,JJ,k+1)
              botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

              hvel(ii,JJ,k) = (1. - botfn) * h_harm + botfn * h_arith
              dz_vel(ii,JJ,k) = (1. - botfn) * dz_harm(ii,JJ,k) + botfn * dz_arith
            endif

            z_i(ii,JJ,k) = z_i(ii,JJ,k+1) + dz_harm(ii,JJ,k) * I_Hbbl(ii,JJ)
          else
            zcol(ii,JJ) = zcol(ii,JJ) + dz(i,j,k)
            zcol_p1(ii,JJ) = zcol_p1(ii,JJ) + dz(i,j+1,k)

            zh(ii,JJ) = zh(ii,JJ) + dz_harm(ii,JJ,k)

            z_clear = max(zcol(ii,JJ), zcol_p1(ii,JJ)) + Dmin(ii,JJ)
            if (zi_dir(ii,JJ) < 0) z_clear = zcol(ii,JJ) + Dmin(ii,JJ)
            if (zi_dir(ii,JJ) > 0) z_clear = zcol_p1(ii,JJ) + Dmin(ii,JJ)

            z_i(ii,JJ,k) = max(zh(ii,JJ), z_clear) * I_Hbbl(ii,JJ)

            hvel(ii,JJ,k) = h_arith
            dz_vel(ii,JJ,k) = dz_arith

            if (v(i,J,k) * h_delta > 0) then
              if (zh(ii,JJ) * I_Hbbl(ii,JJ) < CS%harm_BL_val) then
                hvel(ii,JJ,k) = h_harm
                dz_vel(ii,JJ,k) = dz_harm(ii,JJ,k)
              else
                z2_wt = 1.
                if (zh(ii,JJ) * I_Hbbl(ii,JJ) < 2. * CS%harm_BL_val) &
                  z2_wt = max(0., min(1., zh(ii,JJ) * I_Hbbl(ii,JJ) * I_valBL - 1.))

                ! TODO: should z_clear be used here?
                z2 = z2_wt * (max(zh(ii,JJ), max(zcol(ii,JJ), zcol_p1(ii,JJ)) + Dmin(ii,JJ)) * I_Hbbl(ii,JJ))
                botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

                hvel(ii,JJ,k) = (1. - botfn) * h_arith + botfn * h_harm
                dz_vel(ii,JJ,k) = (1. - botfn) * dz_arith + botfn * dz_harm(ii,JJ,k)
              endif
            endif
          endif

          if (CS%use_GL90_in_SSW) then
            ! The following block calculates the normalized height above the GL90 BBL
            ! (z_i_gl90), using a harmonic mean between layer thicknesses. For the
            ! GL90 BBL we use simply a constant (Hbbl_gl90). The purpose is that the
            ! GL90 coupling coefficient is zeroed out within Hbbl_gl90, to ensure
            ! that no momentum gets fluxed into vanished layers. The scheme is not
            ! sensitive to the exact value of Hbbl_gl90, as long as it is in a
            ! reasonable range (~1-20 m): large enough to capture vanished layers
            ! over topography, small enough to not contaminate the interior.

            z_i_gl90(ii,JJ,k) = z_i_gl90(ii,JJ,k+1) + dz_harm(ii,JJ,k) * I_Hbbl_gl90
          endif
        enddo
      enddo
    enddo

    call find_coupling_coef(a_cpl, dz_vel, do_ij, dz_harm, bbl_thick, kv_bbl, z_i, &
        h_ml, dt, G, GV, US, CS, visc, Ustar_2d, tv, .false., OBC, isb, ieb, JsbB, JebB, nii, nJJB)

    if (allocated(hML_v)) then
      do concurrent (JJ=1:JJe, ii=1:iie) DO_LOCALITY(local(i,J))
        i = isb + ii - 1
        j = JsbB + JJ - 1

        if (do_ij(ii,JJ)) hML_v(i,J) = h_ml(ii,JJ)
      enddo
    endif

    if (CS%use_GL90_in_SSW) then
      call find_coupling_coef_gl90(a_cpl_gl90, dz_vel, do_ij, z_i_gl90, G, GV, &
          CS, VarMix, .false., isb, ieb, JsbB, JebB, nii, nJJB)
    endif

    do_any_shelf = .false.
    if (associated(forces%frac_shelf_v)) then
      do concurrent (JJ=1:JJe, ii=1:iie, do_ij(ii,JJ)) DO_LOCALITY(local(i,J) reduce(.or.:do_any_shelf))
        i = isb + ii - 1
        j = JsbB + JJ - 1

        CS%a1_shelf_v(i,J) = 0.
        do_ij_shelf(ii,JJ) = forces%frac_shelf_v(i,J) > 0.
        do_any_shelf = do_any_shelf .or. do_ij_shelf(ii,JJ)
      enddo
    endif
    if (do_any_shelf) then
      do concurrent (JJ=1:JJe, ii=1:iie, do_ij_shelf(ii,JJ)) DO_LOCALITY(local(i,J))
        i = isb + ii - 1
        j = JsbB + JJ - 1

        ! Initialize non-harmonic depths
        if (.not. CS%harmonic_visc) then
          zh(ii,JJ) = 0.
          Ztop_min(ii,JJ) = min(zcol(ii,JJ), zcol_p1(ii,JJ))
          I_HTbl(ii,JJ) = 1. / (visc%tbl_thick_shelf_v(i,J) + dz_neglect)
        endif
      enddo
      do concurrent (JJ=1:JJe) DO_LOCALITY(local(k))
        do k=1,nz
          do concurrent (ii=1:iie, do_ij_shelf(ii,JJ)) &
              DO_LOCALITY(local(i,J,h_harm,h_arith,h_delta,dz_arith,z2,z2_wt,topfn))
            i = isb + ii - 1
            j = JsbB + JJ - 1

            if (CS%harmonic_visc) then
              hvel_shelf(ii,JJ,k) = hvel(ii,JJ,k)
              dz_vel_shelf(ii,JJ,k) = dz_vel(ii,JJ,k)
            else
              ! Find upwind-biased thickness near the surface.
              ! Perhaps this needs to be done more carefully, via find_eta.
              h_harm = 2. * h(i,j,k) * h(i,j+1,k) &
                  / (h(i,j,k) + h(i,j+1,k) + h_neglect)
              h_arith = 0.5 * (h(i,j+1,k) + h(i,j,k))
              h_delta = h(i,j+1,k) - h(i,j,k)
              dz_arith = 0.5 * (dz(i,j+1,k) + dz(i,j,k))

              ! Project thickness outward across OBCs using a zero-gradient condition.
              if (associated(OBC)) then
                if (OBC%v_N_OBCs_on_PE) then
                  if (OBC%segnum_v(i,J) > 0) then
                    h_harm = h(i,j,k)
                    h_arith = h(i,j,k)
                    h_delta = 0.
                    dz_arith = dz(i,j,k)
                  endif
                endif

                if (OBC%v_S_OBCs_on_PE) then
                  if (OBC%segnum_v(i,J) < 0) then
                    h_harm = h(i,j+1,k)
                    h_arith = h(i,j+1,k)
                    h_delta = 0.
                    dz_arith = dz(i,j+1,k)
                  endif
                endif
              endif

              zcol(ii,JJ) = zcol(ii,JJ) - dz(i,j,k)
              zcol_p1(ii,JJ) = zcol_p1(ii,JJ) - dz(i,j+1,k)

              zh(ii,JJ) = zh(ii,JJ) + dz_harm(ii,JJ,k)

              hvel_shelf(ii,JJ,k) = hvel(ii,JJ,k)
              dz_vel_shelf(ii,JJ,k) = dz_vel(ii,JJ,k)

              if (v(i,J,k) * h_delta > 0.) then
                if (zh(ii,JJ) * I_HTbl(ii,JJ) < CS%harm_BL_val) then
                  hvel_shelf(ii,JJ,k) = min(hvel(ii,JJ,k), h_harm)
                  dz_vel_shelf(ii,JJ,k) = min(dz_vel(ii,JJ,k), dz_harm(ii,JJ,k))
                else
                  z2_wt = 1.
                  if (zh(ii,JJ) * I_HTbl(ii,JJ) < 2. * CS%harm_BL_val) &
                    z2_wt = max(0., min(1., zh(ii,JJ) * I_HTbl(ii,JJ) * I_valBL - 1.))

                  z2 = z2_wt * (max(zh(ii,JJ), Ztop_min(ii,JJ) - min(zcol(ii,JJ), zcol_p1(ii,JJ))) * I_HTbl(ii,JJ))
                  ! TODO: Replace **6
                  topfn = 1. / (1. + 0.09 * z2**6)

                  hvel_shelf(ii,JJ,k) = min(hvel(ii,JJ,k), (1. - topfn) * h_arith + topfn * h_harm)
                  dz_vel_shelf(ii,JJ,k) = min(dz_vel(ii,JJ,k), (1. - topfn) * dz_arith + topfn * dz_harm(ii,JJ,k))
                endif
              endif
            endif
          enddo
        enddo
      enddo

      call find_coupling_coef(a_shelf, dz_vel_shelf, do_ij_shelf, dz_harm, &
          bbl_thick, kv_bbl, z_i, h_ml, dt, G, GV, US, CS, visc, Ustar_2d, &
          tv, .false., OBC, isb, ieb, JsbB, JebB, nii, nJJB, shelf=.true.)

      do concurrent (JJ=1:JJe, ii=1:iie) DO_LOCALITY(local(i,J))
        i = isb + ii - 1
        j = JsbB + JJ - 1

        if (do_ij_shelf(ii,JJ)) CS%a1_shelf_v(i,J) = a_shelf(ii,JJ,1)
      enddo
    endif

    do concurrent (K=1:nz+1, JJ=1:JJe, ii=1:iie, do_ij(ii,JJ)) DO_LOCALITY(local(i,J))
      i = isb + ii - 1
      j = JsbB + JJ - 1

      if (do_ij_shelf(ii,JJ)) then
        ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
        ! CS%a_v(i,J,K) = min(a_cpl_max, forces%frac_shelf_v(i,J) * max(a_shelf(ii,JJ,K), a_cpl(ii,JJ,K)) + &
        !                                (1. - forces%frac_shelf_v(i,J)) * a_cpl(ii,JJ,K))
        if (CS%use_GL90_in_SSW) then
          CS%a_v(i,J,K) = min(a_cpl_max, (forces%frac_shelf_v(i,J) * a_shelf(ii,JJ,K) + &
              (1. - forces%frac_shelf_v(i,J)) * a_cpl(ii,JJ,K)) + a_cpl_gl90(ii,JJ,K))
          CS%a_v_gl90(i,J,K) = min(a_cpl_max, a_cpl_gl90(ii,JJ,K))
        else
          CS%a_v(i,J,K) = min(a_cpl_max, (forces%frac_shelf_v(i,J) * a_shelf(ii,JJ,K) + &
              (1. - forces%frac_shelf_v(i,J)) * a_cpl(ii,JJ,K)))
        endif
      else
        if (CS%use_GL90_in_SSW) then
          a_cpl(ii,JJ,K) = a_cpl(ii,JJ,K) + a_cpl_gl90(ii,JJ,K)
          CS%a_v_gl90(i,J,K) = min(a_cpl_max, a_cpl_gl90(ii,JJ,K))
        endif
        CS%a_v(i,J,K) = min(a_cpl_max, a_cpl(ii,JJ,K))
      endif
    enddo
    do concurrent (k=1:nz, JJ=1:JJe, ii=1:iie, do_ij(ii,JJ)) DO_LOCALITY(local(i,J))
      i = isb + ii - 1
      j = JsbB + JJ - 1

      if (do_ij_shelf(ii,JJ)) then
        ! Should we instead take the inverse of the average of the inverses?
        CS%h_v(i,J,k) = forces%frac_shelf_v(i,J) * hvel_shelf(ii,JJ,k) &
            + (1. - forces%frac_shelf_v(i,J)) * hvel(ii,JJ,k) + h_neglect
      else
        CS%h_v(i,J,k) = hvel(ii,JJ,k) + h_neglect
      endif
    enddo

    ! Diagnose total Kv at v-points
    if (allocated(Kv_v)) then
      do concurrent (k=1:nz, JJ=1:JJe, ii=1:iie, do_ij(ii,JJ)) DO_LOCALITY(local(i,J))
        i = isb + ii - 1
        j = JsbB + JJ - 1

        Kv_v(i,J,k) = 0.5 * (CS%a_v(i,J,K) + CS%a_v(i,J,K+1)) * CS%h_v(i,J,k)
      enddo
    endif

    ! Diagnose GL90 Kv at v-points
    if (allocated(Kv_gl90_v)) then
      do concurrent (k=1:nz, JJ=1:JJe, ii=1:iie, do_ij(ii,JJ)) DO_LOCALITY(local(i,J))
        i = isb + ii - 1
        j = JsbB + JJ - 1

        Kv_gl90_v(i,J,k) = 0.5 * (CS%a_v_gl90(i,J,K) + CS%a_v_gl90(i,J,K+1)) * CS%h_v(i,J,k)
      enddo
    endif
  enddo ; enddo ! Horizontal blocks
  !$omp target exit data &
  !$omp   map(release: hvel, dz_harm, dz_vel, hvel_shelf, dz_vel_shelf, z_i, z_i_gl90, a_cpl, a_cpl_gl90, &
  !$omp     a_shelf, kv_bbl, bbl_thick, I_Hbbl, I_HTbl, Ztop_min, Dmin, zh, zcol, zcol_p1, h_ml, zi_dir, &
  !$omp     do_ij, do_ij_shelf)
end subroutine vertvisc_coef_v


!> Calculate the 'coupling coefficient' (a_cpl) at the interfaces.
!! If BOTTOMDRAGLAW is defined, the minimum of Hbbl and half the adjacent
!! layer thicknesses are used to calculate a_cpl near the bottom.
subroutine find_coupling_coef(a_cpl, hvel, do_ij, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                              dt, G, GV, US, CS, visc, Ustar_2d, tv, work_on_u, OBC, &
                              isb, ieb, jsb, jeb, nii, njj, shelf)
  type(ocean_grid_type),     intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),   intent(in)  :: GV !< Ocean vertical grid structure
  type(unit_scale_type),     intent(in)  :: US !< A dimensional unit scaling type
  integer,                  intent(in)  :: nii !< Zonal block storage extent
  integer,                  intent(in)  :: njj !< Meridional block storage extent
  integer,                  intent(in)  :: isb !< Start of the block's global i index range
  integer,                  intent(in)  :: ieb !< End of the block's global i index range
  integer,                  intent(in)  :: jsb !< Start of the block's global j index range
  integer,                  intent(in)  :: jeb !< End of the block's global j index range
  real, dimension(nii,njj,SZK_(GV)+1), &
                             intent(inout) :: a_cpl !< Coupling coefficient [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(nii,njj,SZK_(GV)), &
                             intent(in)  :: hvel !< Distance between interfaces at velocity points [Z ~> m]
  logical, intent(in) :: do_ij(nii,njj) !< If true, calculate the column's coupling coefficients
  real, dimension(nii,njj,SZK_(GV)), &
                             intent(in)  :: h_harm !< Harmonic mean of thicknesses around a velocity
                                                   !! grid point [Z ~> m]
  real, intent(in)  :: bbl_thick(nii,njj) !< Bottom boundary layer thickness [Z ~> m]
  real, intent(in)  :: kv_bbl(nii,njj) !< Bottom boundary layer viscosity, exclusive of
                                                   !! any depth-dependent contributions from
                                                   !! visc%Kv_shear [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(nii,njj,SZK_(GV)+1), &
                             intent(in)  :: z_i  !< Estimate of interface heights above the bottom,
                                                 !! normalized by the bottom boundary layer thickness [nondim]
  real, intent(inout) :: h_ml(nii,njj) !< Mixed layer depth [Z ~> m]
  real,                      intent(in)  :: dt   !< Time increment [T ~> s]
  type(vertvisc_CS),         intent(in)  :: CS   !< Vertical viscosity control structure
  type(vertvisc_type),       intent(in)  :: visc !< Structure containing viscosities and bottom drag
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)  :: Ustar_2d !< The wind friction velocity, calculated using
                                                 !! the Boussinesq reference density or the
                                                 !! time-evolving surface density in non-Boussinesq
                                                 !! mode [Z T-1 ~> m s-1]
  type(thermo_var_ptrs),     intent(in)  :: tv   !< A structure containing pointers to any available
                                                 !! thermodynamic fields.
  logical,                   intent(in)  :: work_on_u !< If true, u-points are being calculated,
                                                  !! otherwise they are v-points
  type(ocean_OBC_type),      pointer     :: OBC   !< Open boundary condition structure
  logical,         optional, intent(in)  :: shelf !< If present and true, use a surface boundary
                                                  !! condition appropriate for an ice shelf.

  ! Local variables

  real :: u_star(nii,njj)
    ! ustar at a velocity point [Z T-1 ~> m s-1]
  real :: tau_mag(nii,njj)
    ! Wind stress magnitude including gustiness [H Z T-2 ~> m2 s-2 or Pa]
  real :: absf(nii,njj)
    ! The average of the neighboring absolute values of f [T-1 ~> s-1].
  real :: z_t(nii,njj)
    ! The distance from the top, sometimes normalized by Hmix [Z ~> m] or [nondim].
  real :: kv_TBL(nii,njj)
    ! The viscosity in a top boundary layer under ice [H Z T-1 ~> m2 s-1 or Pa s]
  real :: tbl_thick(nii,njj)
    ! The thickness of the top boundary layer [Z ~> m]
  real :: Kv_add
    ! A viscosity to add [H Z T-1 ~> m2 s-1 or Pa s]
  real :: Kv_tot
    ! The total viscosity at an interface [H Z T-1 ~> m2 s-1 or Pa s]
  integer :: nk_in_ml(nii,njj)
    ! The index of the deepest interface in the mixed layer.
  real :: rho_av1 ! The harmonic mean surface layer density at velocity points [R ~> kg m-3]
  real :: h_shear ! The distance over which shears occur [Z ~> m].
  real :: dhc     ! The distance between the center of adjacent layers [Z ~> m].
  real :: visc_ml ! The mixed layer viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real :: I_Hmix  ! The inverse of the mixed layer thickness [Z-1 ~> m-1].
  real :: a_ml    ! The layer coupling coefficient across an interface in
                  ! the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
  real :: a_floor ! A lower bound on the layer coupling coefficient across an interface in
                  ! the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
  real :: I_amax  ! The inverse of the maximum coupling coefficient [T H-1 ~> s m-1 or s m2 kg-1].
  real :: temp1   ! A temporary variable [Z2 ~> m2]
  real :: ustar2_denom ! A temporary variable in the surface boundary layer turbulence
                  ! calculations [H Z-1 T-1 ~> s-1 or kg m-3 s-1]
  real :: h_neglect ! A vertical distance that is so small it is usually lost
                  ! in roundoff and can be neglected [Z ~> m].
  real :: z2      ! A copy of z_i [nondim]
  real :: botfn   ! A function that is 1 at the bottom and small far from it [nondim]
  real :: topfn   ! A function that is 1 at the top and small far from it [nondim]
  real :: kv_top  ! A viscosity associated with the top boundary layer [H Z T-1 ~> m2 s-1 or Pa s]
  logical :: do_shelf, do_OBCs
  integer :: i, j, k, ii, jj, iie, jje
  integer :: nz, max_nk

  nz = GV%ke
  iie = ieb - isb + 1
  jje = jeb - jsb + 1

  !$omp target enter data map(alloc: u_star, tau_mag, absf, z_t, kv_TBL, tbl_thick, nk_in_ml)

  h_neglect = GV%dZ_subroundoff

  if (CS%answer_date < 20190101) then
    !   The maximum coupling coefficient was originally introduced to avoid
    ! truncation error problems in the tridiagonal solver. Effectively, the 1e-10
    ! sets the maximum coupling coefficient increment to 1e10 m per timestep.
    I_amax = (1.0e-10*GV%H_to_m) * dt
  else
    I_amax = 0.0
  endif

  do_shelf = .false. ; if (present(shelf)) do_shelf = shelf

  do_OBCs = .false.
  if (associated(OBC)) then
    if (work_on_u) then
      do_OBCS = OBC%u_E_OBCs_on_PE .or. OBC%u_W_OBCs_on_PE
    else
      do_OBCS = OBC%v_N_OBCs_on_PE .or. OBC%v_S_OBCs_on_PE
    endif
  endif

  do concurrent (K=1:nz+1, jj=1:jje, ii=1:iie, do_ij(ii,jj))
    a_cpl(ii,jj,K) = 0.
  enddo

  I_Hmix = 0.
  if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) I_Hmix = 1. / (CS%Hmix + h_neglect)

  do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj))
    h_ml(ii,jj) = 0.
    if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) z_t(ii,jj) = h_neglect * I_Hmix
  enddo

  do concurrent (jj=1:jje) DO_LOCALITY(local(k))
    do K=2,nz
      do concurrent (ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(i,j,Kv_tot,Kv_add,z2,botfn,dhc,h_shear))
        i = isb + ii - 1
        j = jsb + jj - 1

        Kv_tot = CS%Kv
        if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) then
          ! This is an older (vintage ~1997) way to prevent wind stresses from driving very
          ! large flows in nearly massless near-surface layers when there is not a physically-
          ! based surface boundary layer parameterization.  It does not have a plausible
          ! physical basis, and probably should not be used.
          z_t(ii,jj) = z_t(ii,jj) + h_harm(ii,jj,k-1) * I_Hmix
          Kv_tot = CS%Kv + CS%Kvml_invZ2 / ((z_t(ii,jj) * z_t(ii,jj)) * &
              (1. + 0.09 * z_t(ii,jj) * z_t(ii,jj) * z_t(ii,jj) * z_t(ii,jj) * z_t(ii,jj) * z_t(ii,jj)))
        endif

        if (associated(visc%Kv_shear)) then
          ! Add in viscosities that are determined by physical processes that are handled in
          ! other modules, and which do not respond immediately to the changing layer thicknesses.
          ! These processes may include shear-driven mixing or contributions from some boundary
          ! layer turbulence schemes.  Other viscosity contributions that respond to the evolving
          ! layer thicknesses or the surface wind stresses are added later.
          if (work_on_u) then
            Kv_add = 0.5 * (visc%Kv_shear(i,j,k) + visc%Kv_shear(i+1,j,k))
            if (do_OBCs) then
              if (OBC%u_E_OBCs_on_PE) then
                if (OBC%segnum_u(I,j) > 0) Kv_add = visc%Kv_shear(i,j,k)
              endif
              if (OBC%u_W_OBCs_on_PE) then
                if (OBC%segnum_u(I,j) < 0) Kv_add = visc%Kv_shear(i+1,j,k)
              endif
            endif
            Kv_tot = Kv_tot + Kv_add
          else
            Kv_add = 0.5 * (visc%Kv_shear(i,j,k) + visc%Kv_shear(i,j+1,k))
            if (do_OBCs) then
              if (OBC%v_N_OBCs_on_PE) then
                if (OBC%segnum_v(i,J) > 0) Kv_add = visc%Kv_shear(i,j,k)
              endif
              if (OBC%v_S_OBCs_on_PE) then
                if (OBC%segnum_v(i,J) < 0) Kv_add = visc%Kv_shear(i,j+1,k)
              endif
            endif
            Kv_tot = Kv_tot + Kv_add
          endif
        endif

        if (associated(visc%Kv_shear_Bu)) then
          ! This is similar to what was done above, but for contributions coming from the corner
          ! (vorticity) points.  Because OBCs run through the faces and corners there is no need
          ! to further modify these viscosities here to take OBCs into account.
          if (work_on_u) then
            Kv_tot = Kv_tot + 0.5 * (visc%Kv_shear_Bu(I,J-1,k) + visc%Kv_shear_Bu(I,J,k))
          else
            Kv_tot = Kv_tot + 0.5 * (visc%Kv_shear_Bu(I-1,J,k) + visc%Kv_shear_Bu(I,J,k))
          endif
        endif

        ! Set interior and bottom-enhanced coupling, excluding surface mixed layer contributions.
        if (CS%bottomdraglaw) then
          z2 = z_i(ii,jj,k)
          botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)
          Kv_tot = Kv_tot + (kv_bbl(ii,jj) - CS%Kv) * botfn
          dhc = 0.5 * (hvel(ii,jj,k) + hvel(ii,jj,k-1))
          if (dhc > bbl_thick(ii,jj)) then
            h_shear = ((1. - botfn) * dhc + botfn * bbl_thick(ii,jj)) + h_neglect
          else
            h_shear = dhc + h_neglect
          endif
          a_cpl(ii,jj,K) = Kv_tot / (h_shear + (I_amax * Kv_tot))
        elseif (abs(CS%Kv_extra_bbl) > 0.0) then
          z2 = z_i(ii,jj,k)
          botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)
          Kv_tot = Kv_tot + CS%Kv_extra_bbl * botfn
          h_shear = 0.5 * (hvel(ii,jj,k) + hvel(ii,jj,k-1) + h_neglect)
          a_cpl(ii,jj,K) = Kv_tot / (h_shear + I_amax * Kv_tot)
        else
          h_shear = 0.5 * (hvel(ii,jj,k) + hvel(ii,jj,k-1) + h_neglect)
          a_cpl(ii,jj,K) = Kv_tot / (h_shear + I_amax * Kv_tot)
        endif
      enddo
    enddo
  enddo ! jj-loop over rows in the block

  ! Assign the bottom coupling coefficients.
  do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(dhc))
    if (CS%bottomdraglaw) then
      dhc = hvel(ii,jj,nz) * 0.5
      a_cpl(ii,jj,nz+1) = kv_bbl(ii,jj) / ((min(dhc, bbl_thick(ii,jj)) + h_neglect) + I_amax * kv_bbl(ii,jj))
    elseif (abs(CS%Kv_extra_bbl) > 0.0) then
      a_cpl(ii,jj,nz+1) = (CS%Kv + CS%Kv_extra_bbl) &
          / ((0.5 * hvel(ii,jj,nz) + h_neglect) + I_amax * (CS%Kv + CS%Kv_extra_bbl))
    else
      a_cpl(ii,jj,nz+1) = CS%Kv / ((0.5 * hvel(ii,jj,nz) + h_neglect) + I_amax * CS%Kv)
    endif
  enddo

  ! Add surface intensified viscous coupling, either as a no-slip boundary condition under a
  ! rigid ice-shelf, or due to wind-stress driven surface boundary layer mixing that has not
  ! already been added via visc%Kv_shear.
  if (do_shelf) then
    ! Set the coefficients to include the no-slip surface stress.
    do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(i,j))
      i = isb + ii - 1
      j = jsb + jj - 1

      ! Set the coefficients to include the no-slip surface stress.
      if (work_on_u) then
        kv_TBL(ii,jj) = visc%Kv_tbl_shelf_u(I,j)
        tbl_thick(ii,jj) = visc%tbl_thick_shelf_u(I,j) + h_neglect
      else
        kv_TBL(ii,jj) = visc%Kv_tbl_shelf_v(i,J)
        tbl_thick(ii,jj) = visc%tbl_thick_shelf_v(i,J) + h_neglect
      endif
      z_t(ii,jj) = 0.0

      ! If a_cpl(1) were not already 0, it would be added here.
      if (0.5 * hvel(ii,jj,1) > tbl_thick(ii,jj)) then
        a_cpl(ii,jj,1) = kv_TBL(ii,jj) / (tbl_thick(ii,jj) + I_amax * kv_TBL(ii,jj))
      else
        a_cpl(ii,jj,1) = kv_TBL(ii,jj) / ((0.5 * hvel(ii,jj,1) + h_neglect) + I_amax * kv_TBL(ii,jj))
      endif
    enddo

    do concurrent (jj=1:jje) DO_LOCALITY(local(k))
      do K=2,nz
        do concurrent (ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(topfn,dhc,h_shear,kv_top))
          z_t(ii,jj) = z_t(ii,jj) + hvel(ii,jj,k-1) / tbl_thick(ii,jj)
          topfn = 1. / (1. + 0.09 * z_t(ii,jj)**6)

          dhc = 0.5 * (hvel(ii,jj,k) + hvel(ii,jj,k-1))

          if (dhc > tbl_thick(ii,jj)) then
            h_shear = ((1. - topfn) * dhc + topfn * tbl_thick(ii,jj)) + h_neglect
          else
            h_shear = dhc + h_neglect
          endif

          kv_top = topfn * kv_TBL(ii,jj)
          a_cpl(ii,jj,K) = a_cpl(ii,jj,K) + kv_top / (h_shear + I_amax * kv_top)
        enddo
      enddo
    enddo
  elseif (CS%dynamic_viscous_ML .or. (GV%nkml>0) .or. CS%fixed_LOTW_ML .or. CS%apply_LOTW_floor) then

    ! Find friction velocity and absolute Coriolis parameter at each active block point.
    do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(i,j,rho_av1))
      i = isb + ii - 1
      j = jsb + jj - 1

      u_star(ii,jj) = 0.
      tau_mag(ii,jj) = 0.

      if (allocated(tv%SpV_avg)) then
        rho_av1 = 0.

        if (work_on_u) then
          u_star(ii,jj) = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i+1,j))
          rho_av1 = 2. / (tv%SpV_avg(i,j,1) + tv%SpV_avg(i+1,j,1))
          absf(ii,jj) = 0.5 * (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))

          if (do_OBCs) then
            if (OBC%u_E_OBCs_on_PE) then
              if (OBC%segnum_u(I,j) > 0) then
                u_star(ii,jj) = Ustar_2d(i,j)
                rho_av1 = 1. / tv%SpV_avg(i,j,1)
              endif
            endif

            if (OBC%u_W_OBCs_on_PE) then
              if (OBC%segnum_u(I,j) < 0) then
                u_star(ii,jj) = Ustar_2d(i+1,j)
                rho_av1 = 1. / tv%SpV_avg(i+1,j,1)
              endif
            endif
          endif
        else
          u_star(ii,jj) = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i,j+1))
          rho_av1 = 2. / (tv%SpV_avg(i,j,1) + tv%SpV_avg(i,j+1,1))
          absf(ii,jj) = 0.5 * (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))

          if (do_OBCs) then
            if (OBC%v_N_OBCs_on_PE) then
              if (OBC%segnum_v(i,J) > 0) then
                u_star(ii,jj) = Ustar_2d(i,j)
                rho_av1 = 1. / tv%SpV_avg(i,j,1)
              endif
            endif

            if (OBC%v_S_OBCs_on_PE) then
              if (OBC%segnum_v(i,J) < 0) then
                u_star(ii,jj) = Ustar_2d(i,j+1)
                rho_av1 = 1. / tv%SpV_avg(i,j+1,1)
              endif
            endif
          endif
        endif

        tau_mag(ii,jj) = GV%RZ_to_H * rho_av1 * u_star(ii,jj)**2
      else
        if (work_on_u) then
          u_star(ii,jj) = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i+1,j))
          absf(ii,jj) = 0.5 * (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))

          if (do_OBCs) then
            if (OBC%u_E_OBCs_on_PE) then
              if (OBC%segnum_u(I,j) > 0) u_star(ii,jj) = Ustar_2d(i,j)
            endif

            if (OBC%u_W_OBCs_on_PE) then
              if (OBC%segnum_u(I,j) < 0) u_star(ii,jj) = Ustar_2d(i+1,j)
            endif
          endif
        else
          u_star(ii,jj) = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i,j+1))
          absf(ii,jj) = 0.5 * (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))

          if (do_OBCs) then
            if (OBC%v_N_OBCs_on_PE) then
              if (OBC%segnum_v(i,J) > 0) u_star(ii,jj) = Ustar_2d(i,j)
            endif

            if (OBC%v_S_OBCs_on_PE) then
              if (OBC%segnum_v(i,J) < 0) u_star(ii,jj) = Ustar_2d(i,j+1)
            endif
          endif
        endif

        tau_mag(ii,jj) = GV%Z_to_H * u_star(ii,jj)**2
      endif

      nk_in_ml(ii,jj) = 0
    enddo

    ! Determine the thickness of the surface ocean boundary layer and its extent in index space.
    max_nk = 0

    if (CS%dynamic_viscous_ML) then
      ! The fractional number of layers that are within the viscous boundary layer were
      ! previously stored in visc%nkml_visc_[uv].
      do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(i,j) reduce(max:max_nk))
        i = isb + ii - 1
        j = jsb + jj - 1

        h_ml(ii,jj) = h_neglect
        if (work_on_u) then
          nk_in_ml(ii,jj) = ceiling(visc%nkml_visc_u(I,j))
        else
          nk_in_ml(ii,jj) = ceiling(visc%nkml_visc_v(i,J))
        endif
        max_nk = max(max_nk, nk_in_ml(ii,jj))
      enddo

      do concurrent (jj=1:jje) DO_LOCALITY(local(k))
        do k=1,max_nk
          do concurrent (ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(i,j))
            i = isb + ii - 1
            j = jsb + jj - 1

            if (work_on_u) then
              if (k <= visc%nkml_visc_u(I,j)) then ! This layer is all in the ML.
                h_ml(ii,jj) = h_ml(ii,jj) + hvel(ii,jj,k)
              elseif (k < visc%nkml_visc_u(I,j) + 1.) then ! Part of this layer is in the ML.
                h_ml(ii,jj) = h_ml(ii,jj) + ((visc%nkml_visc_u(I,j) + 1.) - k) * hvel(ii,jj,k)
              endif
            else
              if (k <= visc%nkml_visc_v(i,J)) then
                h_ml(ii,jj) = h_ml(ii,jj) + hvel(ii,jj,k)
              elseif (k < visc%nkml_visc_v(i,J) + 1.) then
                h_ml(ii,jj) = h_ml(ii,jj) + ((visc%nkml_visc_v(i,J) + 1.) - k) * hvel(ii,jj,k)
              endif
            endif
          enddo
        enddo
      enddo
    elseif (GV%nkml>0) then
      ! This is a simple application of a refined-bulk mixed layer with GV%nkml sublayers.
      max_nk = GV%nkml

      do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj))
        nk_in_ml(ii,jj) = GV%nkml
        h_ml(ii,jj) = h_neglect
      enddo

      do concurrent (jj=1:jje) DO_LOCALITY(local(k))
        do k=1,max_nk
          do concurrent (ii=1:iie)
            if (do_ij(ii,jj)) h_ml(ii,jj) = h_ml(ii,jj) + hvel(ii,jj,k)
          enddo
        enddo
      enddo
    elseif (CS%fixed_LOTW_ML .or. CS%apply_LOTW_floor) then
      ! Determine which interfaces are within CS%Hmix of the surface, and set the viscous
      ! boundary layer thickness to the smaller of CS%Hmix and the depth of the ocean.
      do concurrent (jj=1:jje, ii=1:iie)
        if (do_ij(ii,jj)) h_ml(ii,jj) = 0.0
      enddo

      ! Once a column reaches Hmix, its state is unchanged by the remaining layers.
      do concurrent (jj=1:jje) DO_LOCALITY(local(k))
        do k=1,nz
          do concurrent (ii=1:iie, do_ij(ii,jj))
            if (h_ml(ii,jj) < CS%Hmix) then
              nk_in_ml(ii,jj) = k

              if (h_ml(ii,jj) + hvel(ii,jj,k) < CS%Hmix) then
                h_ml(ii,jj) = h_ml(ii,jj) + hvel(ii,jj,k)
              else
                h_ml(ii,jj) = CS%Hmix
              endif
            endif
          enddo
        enddo
      enddo

      do concurrent (jj=1:jje, ii=1:iie) DO_LOCALITY(reduce(max:max_nk))
        if (do_ij(ii,jj)) max_nk = max(max_nk, nk_in_ml(ii,jj))
      enddo
    endif

    ! Avoid working on columns where the viscous coupling could not be increased.
    do concurrent (jj=1:jje, ii=1:iie, do_ij(ii,jj))
      if (u_star(ii,jj) <= 0.) nk_in_ml(ii,jj) = 0

      z_t(ii,jj) = 0.
    enddo

    ! Set the viscous coupling at the interfaces as the larger of what was previously
    ! set and the contributions from the surface boundary layer.
    if (CS%apply_LOTW_floor .and. &
        (CS%dynamic_viscous_ML .or. GV%nkml > 0 .or. CS%fixed_LOTW_ML)) then
      do concurrent (jj=1:jje) DO_LOCALITY(local(k))
        do K=2,max_nk
          do concurrent (ii=1:iie, do_ij(ii,jj)) &
              DO_LOCALITY(local(temp1,ustar2_denom,visc_ml,a_ml,a_floor))
            if (k <= nk_in_ml(ii,jj)) then
              z_t(ii,jj) = z_t(ii,jj) + hvel(ii,jj,k-1)
              !   The viscosity in visc_ml is set to go to 0 at the mixed layer top and bottom
              ! (in a log-layer) and be further limited by rotation to give the natural Ekman length.
              temp1 = (z_t(ii,jj) * h_ml(ii,jj) - z_t(ii,jj) * z_t(ii,jj))
              if (GV%Boussinesq) then
                ustar2_denom = (CS%vonKar * GV%Z_to_H * u_star(ii,jj)**2) &
                    / (absf(ii,jj) * temp1 + (h_ml(ii,jj) + h_neglect) * u_star(ii,jj))
              else
                ustar2_denom = (CS%vonKar * tau_mag(ii,jj)) &
                    / (absf(ii,jj) * temp1 + (h_ml(ii,jj) + h_neglect) * u_star(ii,jj))
              endif
              visc_ml = temp1 * ustar2_denom
              ! Set the viscous coupling based on the model's vertical resolution.  The omission of
              ! the I_amax factor here is consistent with answer dates above 20190101.
              a_ml = visc_ml / (0.25 * (hvel(ii,jj,k) + hvel(ii,jj,k-1) + h_neglect))
              ! As a floor on the viscous coupling, assume that the length scale in the denominator can
              ! not be larger than the distance from the surface, consistent with a logarithmic velocity
              ! profile.  This is consistent with visc_ml, but cancels out common factors of z_t.
              a_floor = (h_ml(ii,jj) - z_t(ii,jj)) * ustar2_denom
              ! Choose the largest estimate of a_cpl.
              a_cpl(ii,jj,K) = max(a_cpl(ii,jj,K), a_ml, a_floor)
              ! An option could be added to change this to:
              ! a_cpl(ii,jj,K) = max(a_cpl(ii,jj,K) + a_ml, a_floor)
            endif
          enddo
        enddo
      enddo
    elseif (CS%apply_LOTW_floor) then
      do concurrent (jj=1:jje) DO_LOCALITY(local(k))
        do K=2,max_nk
          do concurrent (ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(temp1,ustar2_denom))
            if (k <= nk_in_ml(ii,jj)) then
              z_t(ii,jj) = z_t(ii,jj) + hvel(ii,jj,k-1)
              temp1 = (z_t(ii,jj) * h_ml(ii,jj) - z_t(ii,jj) * z_t(ii,jj))
              if (GV%Boussinesq) then
                ustar2_denom = (CS%vonKar * GV%Z_to_H * u_star(ii,jj)**2) &
                    / (absf(ii,jj) * temp1 + (h_ml(ii,jj) + h_neglect) * u_star(ii,jj))
              else
                ustar2_denom = (CS%vonKar * tau_mag(ii,jj)) &
                    / (absf(ii,jj) * temp1 + (h_ml(ii,jj) + h_neglect) * u_star(ii,jj))
              endif
              ! As a floor on the viscous coupling, assume that the length scale in the denominator can not
              ! be larger than the distance from the surface, consistent with a logarithmic velocity profile.
              a_cpl(ii,jj,K) = max(a_cpl(ii,jj,K), (h_ml(ii,jj) - z_t(ii,jj)) * ustar2_denom)
            endif
          enddo
        enddo
      enddo
    else
      do concurrent (jj=1:jje) DO_LOCALITY(local(k))
        do K=2,max_nk
          do concurrent (ii=1:iie, do_ij(ii,jj)) DO_LOCALITY(local(temp1,visc_ml,a_ml))
            if (k <= nk_in_ml(ii,jj)) then
              z_t(ii,jj) = z_t(ii,jj) + hvel(ii,jj,k-1)

              temp1 = (z_t(ii,jj) * h_ml(ii,jj) - z_t(ii,jj) * z_t(ii,jj))
              !   This viscosity is set to go to 0 at the mixed layer top and bottom (in a log-layer)
              ! and be further limited by rotation to give the natural Ekman length.
              ! The following expressions are mathematically equivalent.
              if (GV%Boussinesq .or. (CS%answer_date < 20230601)) then
                visc_ml = u_star(ii,jj) * CS%vonKar * (GV%Z_to_H * temp1 * u_star(ii,jj)) &
                    / (absf(ii,jj) * temp1 + (h_ml(ii,jj) + h_neglect) * u_star(ii,jj))
              else
                visc_ml = CS%vonKar * (temp1 * tau_mag(ii,jj)) &
                    / (absf(ii,jj) * temp1 + (h_ml(ii,jj) + h_neglect) * u_star(ii,jj))
              endif
              a_ml = visc_ml / (0.25 * (hvel(ii,jj,k) + hvel(ii,jj,k-1) + h_neglect) + 0.5 * I_amax * visc_ml)

              ! Choose the largest estimate of a_cpl, but these could be changed to be additive.
              a_cpl(ii,jj,K) = max(a_cpl(ii,jj,K), a_ml)
              ! An option could be added to change this to: a_cpl(i,K) = a_cpl(i,K) + a_ml
            endif
          enddo
        enddo
      enddo
    endif
  endif
  !$omp target exit data map(release: u_star, tau_mag, absf, z_t, kv_TBL, tbl_thick, nk_in_ml)
end subroutine find_coupling_coef


!> Velocity components which exceed a threshold for physically reasonable values are truncated,
!! and the running sum of the number of trunctionas within the non-symmetric memory computational
!! domain is incremented.  Optionally, any column with excessive velocities may be sent
!! to a diagnostic reporting subroutine.
subroutine vertvisc_limit_vel(u, v, h, ADp, CDp, forces, visc, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(accel_diag_ptrs),   intent(in)    :: ADp    !< Acceleration diagnostic pointers
  type(cont_diag_ptrs),    intent(in)    :: CDp    !< Continuity diagnostic pointers
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       intent(inout) :: CS     !< Vertical viscosity control structure

  ! Local variables
  real :: CFL              ! The local CFL number [nondim]
  real :: H_report         ! A thickness below which not to report truncations [H ~> m or kg m-2]
  real :: vel_report(SZIB_(G),SZJB_(G))   ! The velocity to report [L T-1 ~> m s-1]
  real :: u_old(SZIB_(G),SZJ_(G),SZK_(GV)) ! The previous u-velocity [L T-1 ~> m s-1]
  real :: v_old(SZI_(G),SZJB_(G),SZK_(GV)) ! The previous v-velocity [L T-1 ~> m s-1]
  logical :: trunc_any, dowrite(SZIB_(G),SZJB_(G))
  logical :: trunc_col(SZIB_(G),SZJB_(G)) ! True if any level in a column exceeds the truncation threshold.
  logical :: do_any_write
  integer :: ntrunc         ! Thread-safe truncation count for one velocity component.
  integer :: nIIB           ! Zonal block size at u points.
  integer :: njj            ! Meridional block size at u points.
  integer :: nii            ! Zonal block size at v points.
  integer :: nJJB           ! Meridional block size at v points.
  integer :: IsbB, IebB, jsb, jeb, isb, ieb, JsbB, JebB ! Horizontal block bounds.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  H_report = 3.0 * GV%Angstrom_H

  !$omp target enter data map(alloc: dowrite, trunc_col, vel_report, u_old, v_old)

  nIIB = CS%niblock
  if (nIIB == 0) nIIB = Ieq - Isq + 1
  njj = CS%njblock
  if (njj == 0) njj = je - js + 1

  ntrunc = 0
  if (len_trim(CS%u_trunc_file) > 0) then
    do_any_write = .false.
    trunc_any = .false.

    do jsb=js,je,njj ; do IsbB=Isq,Ieq,nIIB
      IebB = min(IsbB + nIIB - 1, Ieq)
      jeb = min(jsb + njj - 1, je)

      do concurrent (j=jsb:jeb, I=IsbB:IebB)
        dowrite(I,j) = .false.
        trunc_col(I,j) = .false.
        vel_report(I,j) = 3.0e8 * US%m_s_to_L_T
      enddo

      do concurrent (j=jsb:jeb) DO_LOCALITY(local(k))
        do k=1,nz
          do concurrent (I=IsbB:IebB) DO_LOCALITY(local(CFL))
            if (abs(u(I,j,k)) < CS%vel_underflow) u(I,j,k) = 0.0

            if (u(I,j,k) < 0.0) then
              CFL = (-u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
            else
              CFL = (u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
            endif

            if (CFL > CS%CFL_trunc) trunc_col(I,j) = .true.

            if (CFL > CS%CFL_report) then
              dowrite(I,j) = .true.
              vel_report(I,j) = min(vel_report(I,j), abs(u(I,j,k)))
            endif
          enddo
        enddo
      enddo
    enddo ; enddo

    ! Main's truncation decision covers the whole domain, not individual blocks.
    do concurrent (j=js:je, I=Isq:Ieq) DO_LOCALITY(reduce(.or.: trunc_any, do_any_write))
      trunc_any = trunc_any .or. trunc_col(I,j)
      do_any_write = do_any_write .or. dowrite(I,j)
    enddo

    do jsb=js,je,njj ; do IsbB=Isq,Ieq,nIIB
      IebB = min(IsbB + nIIB - 1, Ieq)
      jeb = min(jsb + njj - 1, je)

      do concurrent (j=jsb:jeb, I=IsbB:IebB, dowrite(I,j))
        u_old(I,j,:) = u(I,j,:)
      enddo
    enddo ; enddo

    if (trunc_any) then
      do jsb=js,je,njj ; do IsbB=Isq,Ieq,nIIB
        IebB = min(IsbB + nIIB - 1, Ieq)
        jeb = min(jsb + njj - 1, je)

        do concurrent (k=1:nz, j=jsb:jeb, I=IsbB:IebB) DO_LOCALITY(reduce(+: ntrunc))
          if ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i+1,j) < -CS%CFL_trunc) then
            u(I,j,k) = (-0.9 * CS%CFL_trunc) * (G%areaT(i+1,j) / (dt * G%dy_Cu(I,j)))

            if (I >= G%isc .and. I <= G%iec .and. j >= G%jsc .and. j <= G%jec &
                .and. CS%h_u(I,j,k) > H_report) &
              ntrunc = ntrunc + 1
          elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i,j) > CS%CFL_trunc) then
            u(I,j,k) = (0.9 * CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dy_Cu(I,j)))

            if (I >= G%isc .and. I <= G%iec .and. j >= G%jsc .and. j <= G%jec &
                .and. CS%h_u(I,j,k) > H_report) &
              ntrunc = ntrunc + 1
          endif
        enddo
      enddo ; enddo
    endif

    CS%ntrunc = CS%ntrunc + ntrunc

    if (do_any_write) then
      !$omp target update from(u_old, vel_report, dowrite, CS%a_u, CS%h_u)
      do j=js,je ; do I=Isq,Ieq ; if (dowrite(I,j)) then
        ! Call a diagnostic reporting subroutines are called if unphysically large values are found.
        call write_u_accel(I, j, u_old, h, ADp, CDp, dt, G, GV, US, CS%PointAccel_CSp, &
                           vel_report(I,j), forces%taux(I,j), a=CS%a_u, hv=CS%h_u)
      endif ; enddo ; enddo
    endif
  else  ! Do not report accelerations leading to large velocities.
    do jsb=js,je,njj ; do IsbB=Isq,Ieq,nIIB
      IebB = min(IsbB + nIIB - 1, Ieq)
      jeb = min(jsb + njj - 1, je)

      do concurrent (k=1:nz, j=jsb:jeb, I=IsbB:IebB) DO_LOCALITY(reduce(+: ntrunc))
        if (abs(u(I,j,k)) < CS%vel_underflow) then
          u(I,j,k) = 0.0
        elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i+1,j) < -CS%CFL_trunc) then
          u(I,j,k) = (-0.9 * CS%CFL_trunc) * (G%areaT(i+1,j) / (dt * G%dy_Cu(I,j)))

          if (I >= G%isc .and. I <= G%iec .and. j >= G%jsc .and. j <= G%jec &
              .and. CS%h_u(I,j,k) > H_report) &
            ntrunc = ntrunc + 1
        elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i,j) > CS%CFL_trunc) then
          u(I,j,k) = (0.9 * CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dy_Cu(I,j)))

          if (I >= G%isc .and. I <= G%iec .and. j >= G%jsc .and. j <= G%jec &
              .and. CS%h_u(I,j,k) > H_report) &
            ntrunc = ntrunc + 1
        endif
      enddo
    enddo ; enddo

    CS%ntrunc = CS%ntrunc + ntrunc
  endif

  nii = CS%niblock
  if (nii == 0) nii = ie - is + 1
  nJJB = CS%njblock
  if (nJJB == 0) nJJB = Jeq - Jsq + 1

  ntrunc = 0
  if (len_trim(CS%v_trunc_file) > 0) then
    do_any_write =.false.
    trunc_any = .false.

    do JsbB=Jsq,Jeq,nJJB ; do isb=is,ie,nii
      ieb = min(isb + nii - 1, ie)
      JebB = min(JsbB + nJJB - 1, Jeq)

      do concurrent (J=JsbB:JebB, i=isb:ieb)
        dowrite(i,J) = .false.
        trunc_col(i,J) = .false.
        vel_report(i,J) = 3.0e8 * US%m_s_to_L_T
      enddo

      do concurrent (J=JsbB:JebB) DO_LOCALITY(local(k))
        do k=1,nz
          do concurrent (i=isb:ieb) DO_LOCALITY(local(CFL))
            if (abs(v(i,J,k)) < CS%vel_underflow) v(i,J,k) = 0.0

            if (v(i,J,k) < 0.0) then
              CFL = (-v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
            else
              CFL = (v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
            endif

            if (CFL > CS%CFL_trunc) trunc_col(i,J) = .true.

            if (CFL > CS%CFL_report) then
              dowrite(i,J) = .true.
              vel_report(i,J) = min(vel_report(i,J), abs(v(i,J,k)))
            endif
          enddo
        enddo
      enddo
    enddo ; enddo

    do concurrent (J=Jsq:Jeq, i=is:ie) DO_LOCALITY(reduce(.or.: trunc_any, do_any_write))
      trunc_any = trunc_any .or. trunc_col(i,J)
      do_any_write = do_any_write .or. dowrite(i,J)
    enddo

    do JsbB=Jsq,Jeq,nJJB ; do isb=is,ie,nii
      ieb = min(isb + nii - 1, ie)
      JebB = min(JsbB + nJJB - 1, Jeq)

      do concurrent (J=JsbB:JebB, i=isb:ieb, dowrite(i,J))
        v_old(i,J,:) = v(i,J,:)
      enddo
    enddo ; enddo

    if (trunc_any) then
      do JsbB=Jsq,Jeq,nJJB ; do isb=is,ie,nii
        ieb = min(isb + nii - 1, ie)
        JebB = min(JsbB + nJJB - 1, Jeq)

        do concurrent (k=1:nz, J=JsbB:JebB, i=isb:ieb) DO_LOCALITY(reduce(+: ntrunc))
          if ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j+1) < -CS%CFL_trunc) then
            v(i,J,k) = (-0.9 * CS%CFL_trunc) * (G%areaT(i,j+1) / (dt * G%dx_Cv(i,J)))

            if (i >= G%isc .and. i <= G%iec .and. J >= G%jsc .and. J <= G%jec &
                .and. CS%h_v(i,J,k) > H_report) &
              ntrunc = ntrunc + 1
          elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j) > CS%CFL_trunc) then
            v(i,J,k) = (0.9 * CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dx_Cv(i,J)))

            if (i >= G%isc .and. i <= G%iec .and. J >= G%jsc .and. J <= G%jec &
                .and. CS%h_v(i,J,k) > H_report) &
              ntrunc = ntrunc + 1
          endif
        enddo
      enddo ; enddo
    endif
    CS%ntrunc = CS%ntrunc + ntrunc

    if (do_any_write) then
      !$omp target update from(v_old, vel_report, dowrite, CS%a_v, CS%h_v)
      do J=Jsq,Jeq ; do i=is,ie ; if (dowrite(i,J)) then
        ! Call a diagnostic reporting subroutines are called if unphysically
        ! large values are found.
        call write_v_accel(i, J, v_old, h, ADp, CDp, dt, G, GV, US, &
            CS%PointAccel_CSp, vel_report(i,J), forces%tauy(i,J), &
            a=CS%a_v, hv=CS%h_v)
      endif ; enddo ; enddo
    endif
  else  ! Do not report accelerations leading to large velocities.
    do JsbB=Jsq,Jeq,nJJB ; do isb=is,ie,nii
      ieb = min(isb + nii - 1, ie)
      JebB = min(JsbB + nJJB - 1, Jeq)

      do concurrent (k=1:nz, J=JsbB:JebB, i=isb:ieb) DO_LOCALITY(reduce(+: ntrunc))
        if (abs(v(i,J,k)) < CS%vel_underflow) then
          v(i,J,k) = 0.0
        elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j+1) < -CS%CFL_trunc) then
          v(i,J,k) = (-0.9 * CS%CFL_trunc) * (G%areaT(i,j+1) / (dt * G%dx_Cv(i,J)))

          if (i >= G%isc .and. i <= G%iec .and. J >= G%jsc .and. J <= G%jec &
              .and. CS%h_v(i,J,k) > H_report) &
            ntrunc = ntrunc + 1
        elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j) > CS%CFL_trunc) then
          v(i,J,k) = (0.9 * CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dx_Cv(i,J)))

          if (i >= G%isc .and. i <= G%iec .and. J >= G%jsc .and. J <= G%jec &
              .and. CS%h_v(i,J,k) > H_report) &
            ntrunc = ntrunc + 1
        endif
      enddo
    enddo ; enddo

    CS%ntrunc = CS%ntrunc + ntrunc
  endif

  !$omp target exit data map(release: u_old, v_old, dowrite, trunc_col, vel_report)
end subroutine vertvisc_limit_vel


!> Initialize the vertical friction module
subroutine vertvisc_init(MIS, Time, G, GV, US, param_file, diag, ADp, dirs, &
                          ntrunc, CS, fpmix)
  type(ocean_internal_state), &
                   target, intent(in)    :: MIS    !< The "MOM Internal State", a set of pointers
                                                   !! to the fields and accelerations that make
                                                   !! up the ocean's physical state
  type(time_type), target, intent(in)    :: Time   !< Current model time
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< File to parse for parameters
  type(diag_ctrl), target, intent(inout) :: diag   !< Diagnostic control structure
  type(accel_diag_ptrs),   intent(inout) :: ADp    !< Acceleration diagnostic pointers
  type(directories),       intent(in)    :: dirs   !< Relevant directory paths
  integer, target,         intent(inout) :: ntrunc !< Number of velocity truncations
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure
  logical, optional,       intent(in)    :: fpmix  !< Nonlocal momentum mixing

  ! Local variables

  real :: Kv_BBL  ! A viscosity in the bottom boundary layer with a simple scheme [H Z T-1 ~> m2 s-1 or Pa s]
  real :: Kv_back_z  ! A background kinematic viscosity [Z2 T-1 ~> m2 s-1]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  logical :: lfpmix
  character(len=200) :: kappa_gl90_file, inputdir, kdgl90_varname
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_vert_friction" ! This module's name.
  character(len=40)  :: thickness_units
  real :: Kv_mks ! KVML in MKS [m2 s-1]

  if (associated(CS)) then
    call MOM_error(WARNING, "vertvisc_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%initialized = .true.

  if (GV%Boussinesq) then ; thickness_units = "m"
  else ; thickness_units = "kg m-2" ; endif

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag ; CS%ntrunc => ntrunc ; ntrunc = 0

  lfpmix = .false.
  if (present(fpmix)) lfpmix = fpmix

! Default, read and log parameters
  call log_version(param_file, mdl, version, "", log_to_all=.true., debugging=.true.)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "VERT_FRICTION_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the viscous "//&
                 "calculations.  Values below 20190101 recover the answers from the end of 2018, "//&
                 "while higher values use expressions that do not use an arbitrary hard-coded "//&
                 "maximum viscous coupling coefficient between layers.  Values below 20230601 "//&
                 "recover a form of the viscosity within the mixed layer that breaks up the "//&
                 "magnitude of the wind stress in some non-Boussinesq cases.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  call get_param(param_file, mdl, "BOTTOMDRAGLAW", CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag "//&
                 "law of the form c_drag*|u|*u. The velocity magnitude "//&
                 "may be an assumed value or it may be based on the "//&
                 "actual velocity in the bottommost HBBL, depending on "//&
                 "LINEAR_DRAG.", default=.true.)
  call get_param(param_file, mdl, "DIRECT_STRESS", CS%direct_stress, &
                 "If true, the wind stress is distributed over the topmost HMIX_STRESS of fluid "//&
                 "(like in HYCOM), and an added mixed layer viscosity or a physically based "//&
                 "boundary layer turbulence parameterization is not needed for stability.", &
                 default=.false.)
  call get_param(param_file, mdl, "DYNAMIC_VISCOUS_ML", CS%dynamic_viscous_ML, &
                 "If true, use a bulk Richardson number criterion to "//&
                 "determine the mixed layer thickness for viscosity.", &
                 default=.false.)
  call get_param(param_file, mdl, "FIXED_DEPTH_LOTW_ML", CS%fixed_LOTW_ML, &
                 "If true, use a Law-of-the-wall prescription for the mixed layer viscosity "//&
                 "within a boundary layer that is the lesser of HMIX_FIXED and the total "//&
                 "depth of the ocean in a column.", default=.false.)
  call get_param(param_file, mdl, "LOTW_VISCOUS_ML_FLOOR", CS%apply_LOTW_floor, &
                 "If true, use a Law-of-the-wall prescription to set a lower bound on the "//&
                 "viscous coupling between layers within the surface boundary layer, based "//&
                 "the distance of interfaces from the surface.  This only acts when there "//&
                 "are large changes in the thicknesses of successive layers or when the "//&
                 "viscosity is set externally and the wind stress has subsequently increased.", &
                 default=.false.)
  call get_param(param_file, mdl, 'VON_KARMAN_CONST', CS%vonKar, &
                 'The value the von Karman constant as used for mixed layer viscosity.', &
                 units='nondim', default=0.41)
  call get_param(param_file, mdl, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to a file into which the accelerations "//&
                 "leading to zonal velocity truncations are written. "//&
                 "Undefine this for efficiency if this diagnostic is not needed.", &
                 default=" ", debuggingParam=.true.)
  call get_param(param_file, mdl, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to a file into which the accelerations "//&
                 "leading to meridional velocity truncations are written. "//&
                 "Undefine this for efficiency if this diagnostic is not needed.", &
                 default=" ", debuggingParam=.true.)
  call get_param(param_file, mdl, "HARMONIC_VISC", CS%harmonic_visc, &
                 "If true, use the harmonic mean thicknesses for "//&
                 "calculating the vertical viscosity.", default=.false.)
  call get_param(param_file, mdl, "HARMONIC_BL_SCALE", CS%harm_BL_val, &
                 "A scale to determine when water is in the boundary "//&
                 "layers based solely on harmonic mean thicknesses for "//&
                 "the purpose of determining the extent to which the "//&
                 "thicknesses used in the viscosities are upwinded.", &
                 default=0.0, units="nondim")
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)

  if (GV%nkml < 1) then
    call get_param(param_file, mdl, "HMIX_FIXED", CS%Hmix, &
                 "The prescribed depth over which the near-surface viscosity and "//&
                 "diffusivity are elevated when the bulk mixed layer is not used.", &
                 units="m", scale=US%m_to_Z, fail_if_missing=.true.)
  endif
  if (CS%direct_stress) then
    if (GV%nkml < 1) then
      call get_param(param_file, mdl, "HMIX_STRESS", CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if DIRECT_STRESS is true.", &
                 units="m", default=US%Z_to_m*CS%Hmix, scale=GV%m_to_H)
    else
      call get_param(param_file, mdl, "HMIX_STRESS", CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if DIRECT_STRESS is true.", &
                 units="m", fail_if_missing=.true., scale=GV%m_to_H)
    endif
    if (CS%Hmix_stress <= 0.0) call MOM_error(FATAL, "vertvisc_init: " // &
       "HMIX_STRESS must be set to a positive value if DIRECT_STRESS is true.")
  endif
  call get_param(param_file, mdl, "KV", Kv_back_z, &
                 "The background kinematic viscosity in the interior. "//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", fail_if_missing=.true., scale=US%m2_s_to_Z2_T)
  ! Convert input kinematic viscosity to dynamic viscosity when non-Boussinesq.
  CS%Kv = (US%Z2_T_to_m2_s*GV%m2_s_to_HZ_T) * Kv_back_z

  call get_param(param_file, mdl, "USE_GL90_IN_SSW", CS%use_GL90_in_SSW, &
                 "If true, use simpler method to calculate 1/N^2 in GL90 vertical "// &
                 "viscosity coefficient. This method is valid in stacked shallow water mode.", &
                 default=.false.)
  call get_param(param_file, mdl, "KD_GL90", CS%kappa_gl90, &
                 "The scalar diffusivity used in GL90 vertical viscosity scheme.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L*US%Z_to_L*GV%m_to_H*US%T_to_s, &
                 do_not_log=.not.CS%use_GL90_in_SSW)
  call get_param(param_file, mdl, "READ_KD_GL90", CS%read_kappa_gl90, &
                 "If true, read a file (given by KD_GL90_FILE) containing the "//&
                 "spatially varying diffusivity KD_GL90 used in the GL90 scheme.", default=.false., &
                 do_not_log=.not.CS%use_GL90_in_SSW)
  if (CS%read_kappa_gl90) then
    if (CS%kappa_gl90 > 0) then
        call MOM_error(FATAL, "MOM_vert_friction.F90, vertvisc_init: KD_GL90 > 0 "// &
              "is not compatible with READ_KD_GL90 = .TRUE. ")
    endif
    call get_param(param_file, mdl, "INPUTDIR", inputdir, &
                 "The directory in which all input files are found.", &
                 default=".", do_not_log=.true.)
    inputdir = slasher(inputdir)
    call get_param(param_file, mdl, "KD_GL90_FILE", kappa_gl90_file, &
                 "The file containing the spatially varying diffusivity used in the "// &
                 "GL90 scheme.", default="kd_gl90.nc", do_not_log=.not.CS%use_GL90_in_SSW)
    call get_param(param_file, mdl, "KD_GL90_VARIABLE", kdgl90_varname, &
                 "The name of the GL90 diffusivity variable to read "//&
                 "from KD_GL90_FILE.", default="kd_gl90", do_not_log=.not.CS%use_GL90_in_SSW)
    kappa_gl90_file = trim(inputdir) // trim(kappa_gl90_file)

    allocate(CS%kappa_gl90_2d(G%isd:G%ied, G%jsd:G%jed), source=0.0)
    call MOM_read_data(kappa_gl90_file, kdgl90_varname, CS%kappa_gl90_2d(:,:), G%domain, &
                       scale=US%m_to_L*US%Z_to_L*GV%m_to_H*US%T_to_s)
    call pass_var(CS%kappa_gl90_2d, G%domain)
  endif
  call get_param(param_file, mdl, "USE_GL90_N2", CS%use_GL90_N2, &
                 "If true, use GL90 vertical viscosity coefficient that is depth-independent; "// &
                 "this corresponds to a kappa_GM that scales as N^2 with depth.", &
                 default=.false., do_not_log=.not.CS%use_GL90_in_SSW)
  if (CS%use_GL90_N2) then
    if (.not. CS%use_GL90_in_SSW) call MOM_error(FATAL, &
           "MOM_vert_friction.F90, vertvisc_init: "//&
           "When USE_GL90_N2=True, USE_GL90_in_SSW must also be True.")
    if (CS%kappa_gl90 > 0) then
        call MOM_error(FATAL, "MOM_vert_friction.F90, vertvisc_init: KD_GL90 > 0 "// &
              "is not compatible with USE_GL90_N2 = .TRUE. ")
    endif
    if (CS%read_kappa_gl90) call MOM_error(FATAL, &
           "MOM_vert_friction.F90, vertvisc_init: "//&
           "READ_KD_GL90 = .TRUE. is not compatible with USE_GL90_N2 = .TRUE.")
    call get_param(param_file, mdl, "alpha_GL90", CS%alpha_gl90, &
                   "Coefficient used to compute a depth-independent GL90 vertical "//&
                   "viscosity via Kv_GL90 = alpha_GL90 * f2. Is only used "// &
                   "if USE_GL90_N2 is true. Note that the implied Kv_GL90 "// &
                   "corresponds to a KD_GL90 that scales as N^2 with depth.", &
                   units="m2 s", default=0.0, scale=GV%m_to_H*US%m_to_Z*US%s_to_T, &
                   do_not_log=.not.CS%use_GL90_in_SSW)
  endif
  call get_param(param_file, mdl, "HBBL_GL90", CS%Hbbl_gl90, &
                 "The thickness of the GL90 bottom boundary layer, "//&
                 "which defines the range over which the GL90 coupling "//&
                 "coefficient is zeroed out, in order to avoid fluxing "//&
                 "momentum into vanished layers over steep topography.", &
                 units="m", default=5.0, scale=US%m_to_Z, do_not_log=.not.CS%use_GL90_in_SSW)

  CS%Kvml_invZ2 = 0.0
  if (GV%nkml < 1) then
    call get_param(param_file, mdl, "KVML", Kv_mks, &
                 "The scale for an extra kinematic viscosity in the mixed layer", &
                 units="m2 s-1", default=-1.0, do_not_log=.true.)
    if (Kv_mks >= 0.0) then
      call MOM_error(WARNING, "KVML is a deprecated parameter. Use KV_ML_INVZ2 instead.")
    else
      Kv_mks = 0.0
    endif
    call get_param(param_file, mdl, "KV_ML_INVZ2", CS%Kvml_invZ2, &
                 "An extra kinematic viscosity in a mixed layer of thickness HMIX_FIXED, "//&
                 "with the actual viscosity scaling as 1/(z*HMIX_FIXED)^2, where z is the "//&
                 "distance from the surface, to allow for finite wind stresses to be "//&
                 "transmitted through infinitesimally thin surface layers.  This is an "//&
                 "older option for numerical convenience without a strong physical basis, "//&
                 "and its use is now discouraged.", &
                 units="m2 s-1", default=Kv_mks, scale=GV%m2_s_to_HZ_T)
  endif

  if (.not.CS%bottomdraglaw) then
    call get_param(param_file, mdl, "KV_EXTRA_BBL", CS%Kv_extra_bbl, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=0.0, scale=GV%m2_s_to_HZ_T, do_not_log=.true.)
    if (CS%Kv_extra_bbl == 0.0) then
      call get_param(param_file, mdl, "KVBBL", Kv_BBL, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=US%Z2_T_to_m2_s*Kv_back_z, scale=GV%m2_s_to_HZ_T, &
                 do_not_log=.true.)
      if (abs(Kv_BBL - CS%Kv) > 1.0e-15*abs(CS%Kv)) then
        call MOM_error(WARNING, "KVBBL is a deprecated parameter. Use KV_EXTRA_BBL instead.")
        CS%Kv_extra_bbl = Kv_BBL - CS%Kv
      endif
    endif
    call log_param(param_file, mdl, "KV_EXTRA_BBL", CS%Kv_extra_bbl, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=0.0, unscale=GV%HZ_T_to_m2_s)
  endif
  call get_param(param_file, mdl, "HBBL", CS%Hbbl, &
                 "The thickness of a bottom boundary layer with a viscosity increased by "//&
                 "KV_EXTRA_BBL if BOTTOMDRAGLAW is not defined, or the thickness over which "//&
                 "near-bottom velocities are averaged for the drag law if BOTTOMDRAGLAW is "//&
                 "defined but LINEAR_DRAG is not.", &
                 units="m", fail_if_missing=.true., scale=US%m_to_Z)
  call get_param(param_file, mdl, "CFL_TRUNCATE", CS%CFL_trunc, &
                 "The value of the CFL number that will cause velocity "//&
                 "components to be truncated; instability can occur past 0.5.", &
                 units="nondim", default=0.5)
  call get_param(param_file, mdl, "CFL_REPORT", CS%CFL_report, &
                 "The value of the CFL number that causes accelerations "//&
                 "to be reported; the default is CFL_TRUNCATE.", &
                 units="nondim", default=CS%CFL_trunc)
  call get_param(param_file, mdl, "CFL_TRUNCATE_RAMP_TIME", CS%truncRampTime, &
                 "The time over which the CFL truncation value is ramped "//&
                 "up at the beginning of the run.", &
                 units="s", default=0., scale=US%s_to_T)
  CS%CFL_truncE = CS%CFL_trunc
  call get_param(param_file, mdl, "CFL_TRUNCATE_START", CS%CFL_truncS, &
                 "The start value of the truncation CFL number used when "//&
                 "ramping up CFL_TRUNC.", &
                 units="nondim", default=0.)
  call get_param(param_file, mdl, "STOKES_MIXING_COMBINED", CS%StokesMixing, &
                 "Flag to use Stokes drift Mixing via the Lagrangian "//&
                 " current (Eulerian plus Stokes drift). "//&
                 " Still needs work and testing, so not recommended for use.",&
                 default=.false.)
  !BGR 04/04/2018{
  ! StokesMixing is required for MOM6 for some Langmuir mixing parameterization.
  !   The code used here has not been developed for vanishing layers or in
  !   conjunction with any bottom friction.  Therefore, the following line is
  !   added so this functionality cannot be used without user intervention in
  !   the code.  This will prevent general use of this functionality until proper
  !   care is given to the previously mentioned issues.  Comment out the following
  !   MOM_error to use, but do so at your own risk and with these points in mind.
  !}
  if (CS%StokesMixing) then
    call MOM_error(FATAL, "Stokes mixing requires user intervention in the code.\n"//&
                          "  Model now exiting.  See MOM_vert_friction.F90 for \n"//&
                          "  details (search 'BGR 04/04/2018' to locate comment).")
  endif
  call get_param(param_file, mdl, "VEL_UNDERFLOW", CS%vel_underflow, &
                 "A negligibly small velocity magnitude below which velocity "//&
                 "components are set to 0.  A reasonable value might be "//&
                 "1e-30 m/s, which is less than an Angstrom divided by "//&
                 "the age of the universe.", units="m s-1", default=0.0, scale=US%m_s_to_L_T)

  call get_param(param_file, mdl, "VERT_FRICTION_NIBLOCK", CS%niblock, &
      "The i-direction block size for vertical viscosity calculations. " // &
      "A value of 0 uses the full computational i-domain extent.", &
      default=0, layoutParam=.true.)
  call get_param(param_file, mdl, "VERT_FRICTION_NJBLOCK", CS%njblock, &
      "The j-direction block size for vertical viscosity calculations. " // &
      "A value of 0 uses the full computational j-domain extent.", &
      default=1, layoutParam=.true.)
  if (CS%niblock < 0) &
    call MOM_error(FATAL, &
        "MOM_vert_friction_init: VERT_FRICTION_NIBLOCK must be nonnegative.")
  if (CS%njblock < 0) &
    call MOM_error(FATAL, &
        "MOM_vert_friction_init: VERT_FRICTION_NJBLOCK must be nonnegative.")

  ALLOC_(CS%a_u(IsdB:IedB,jsd:jed,nz+1)) ; CS%a_u(:,:,:) = 0.0
  ALLOC_(CS%a_u_gl90(IsdB:IedB,jsd:jed,nz+1)) ; CS%a_u_gl90(:,:,:) = 0.0
  ALLOC_(CS%h_u(IsdB:IedB,jsd:jed,nz))   ; CS%h_u(:,:,:) = 0.0
  ALLOC_(CS%a_v(isd:ied,JsdB:JedB,nz+1)) ; CS%a_v(:,:,:) = 0.0
  ALLOC_(CS%a_v_gl90(isd:ied,JsdB:JedB,nz+1)) ; CS%a_v_gl90(:,:,:) = 0.0
  ALLOC_(CS%h_v(isd:ied,JsdB:JedB,nz))   ; CS%h_v(:,:,:) = 0.0

  CS%id_Kv_slow = register_diag_field('ocean_model', 'Kv_slow', diag%axesTi, Time, &
      'Slow varying vertical viscosity', 'm2 s-1', conversion=GV%HZ_T_to_m2_s)

  CS%id_Kv_u = register_diag_field('ocean_model', 'Kv_u', diag%axesCuL, Time, &
      'Total vertical viscosity at u-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_v = register_diag_field('ocean_model', 'Kv_v', diag%axesCvL, Time, &
      'Total vertical viscosity at v-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_gl90_u = register_diag_field('ocean_model', 'Kv_gl90_u', diag%axesCuL, Time, &
      'GL90 vertical viscosity at u-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_gl90_v = register_diag_field('ocean_model', 'Kv_gl90_v', diag%axesCvL, Time, &
      'GL90 vertical viscosity at v-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_au_vv = register_diag_field('ocean_model', 'au_visc', diag%axesCui, Time, &
      'Zonal Viscous Vertical Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_av_vv = register_diag_field('ocean_model', 'av_visc', diag%axesCvi, Time, &
      'Meridional Viscous Vertical Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_au_gl90_vv = register_diag_field('ocean_model', 'au_gl90_visc', diag%axesCui, Time, &
      'Zonal Viscous Vertical GL90 Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_av_gl90_vv = register_diag_field('ocean_model', 'av_gl90_visc', diag%axesCvi, Time, &
      'Meridional Viscous Vertical GL90 Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_h_u = register_diag_field('ocean_model', 'Hu_visc', diag%axesCuL, Time, &
      'Thickness at Zonal Velocity Points for Viscosity', &
      thickness_units, conversion=GV%H_to_MKS)
      ! Alternately, to always give this variable in 'm' use the following line instead:
      ! 'm', conversion=GV%H_to_m)

  CS%id_h_v = register_diag_field('ocean_model', 'Hv_visc', diag%axesCvL, Time, &
      'Thickness at Meridional Velocity Points for Viscosity', &
      thickness_units, conversion=GV%H_to_MKS)

  CS%id_hML_u = register_diag_field('ocean_model', 'HMLu_visc', diag%axesCu1, Time, &
      'Mixed Layer Thickness at Zonal Velocity Points for Viscosity', &
      thickness_units, conversion=US%Z_to_m)

  CS%id_hML_v = register_diag_field('ocean_model', 'HMLv_visc', diag%axesCv1, Time, &
      'Mixed Layer Thickness at Meridional Velocity Points for Viscosity', &
      thickness_units, conversion=US%Z_to_m)

 if (lfpmix) then
  CS%id_uE_h = register_diag_field('ocean_model', 'uE_h' , CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vE_h = register_diag_field('ocean_model', 'vE_h' , CS%diag%axesTL, &
      Time, 'y-merid Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_uInc_h = register_diag_field('ocean_model','uInc_h',CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vInc_h = register_diag_field('ocean_model','vInc_h',CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_uStk = register_diag_field('ocean_model', 'uStk' , CS%diag%axesTL, &
      Time, 'x-FP du increment' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vStk = register_diag_field('ocean_model', 'vStk' , CS%diag%axesTL, &
      Time, 'y-FP dv increment' , 'm s-1', conversion=US%L_T_to_m_s)

  CS%id_FPtau2s = register_diag_field('ocean_model','Omega_tau2s',CS%diag%axesTi, &
      Time, 'Stress direction from shear','radians')
  CS%id_FPtau2w = register_diag_field('ocean_model','Omega_tau2w',CS%diag%axesTi, &
      Time, 'Stress direction from wind','radians')
  CS%id_uStk0 = register_diag_field('ocean_model', 'uStk0' , diag%axesT1, &
      Time, 'Zonal Surface Stokes', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vStk0 = register_diag_field('ocean_model', 'vStk0' , diag%axesT1, &
      Time, 'Merid Surface Stokes', 'm s-1', conversion=US%L_T_to_m_s)
  endif

  CS%id_du_dt_visc = register_diag_field('ocean_model', 'du_dt_visc', diag%axesCuL, Time, &
      'Zonal Acceleration from Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_visc > 0) call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  CS%id_dv_dt_visc = register_diag_field('ocean_model', 'dv_dt_visc', diag%axesCvL, Time, &
      'Meridional Acceleration from Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_visc > 0) call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  CS%id_GLwork = register_diag_field('ocean_model', 'GLwork', diag%axesTL, Time, &
      'Sign-definite Kinetic Energy Source from GL90 Vertical Viscosity', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_du_dt_visc_gl90 = register_diag_field('ocean_model', 'du_dt_visc_gl90', diag%axesCuL, Time, &
      'Zonal Acceleration from GL90 Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_du_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    call safe_alloc_ptr(ADp%du_dt_visc_gl90,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  endif
  CS%id_dv_dt_visc_gl90 = register_diag_field('ocean_model', 'dv_dt_visc_gl90', diag%axesCvL, Time, &
      'Meridional Acceleration from GL90 Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_dv_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    call safe_alloc_ptr(ADp%dv_dt_visc_gl90,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  endif
  CS%id_du_dt_str = register_diag_field('ocean_model', 'du_dt_str', diag%axesCuL, Time, &
      'Zonal Acceleration from Surface Wind Stresses', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_str > 0) call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
  CS%id_dv_dt_str = register_diag_field('ocean_model', 'dv_dt_str', diag%axesCvL, Time, &
      'Meridional Acceleration from Surface Wind Stresses', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_str > 0) call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)

  CS%id_taux_bot = register_diag_field('ocean_model', 'taux_bot', diag%axesCu1, &
      Time, 'Zonal Bottom Stress from Ocean to Earth', &
      'Pa', conversion=US%RZ_to_kg_m2*US%L_T2_to_m_s2)
  CS%id_tauy_bot = register_diag_field('ocean_model', 'tauy_bot', diag%axesCv1, &
      Time, 'Meridional Bottom Stress from Ocean to Earth', &
      'Pa', conversion=US%RZ_to_kg_m2*US%L_T2_to_m_s2)

  !CS%id_hf_du_dt_visc = register_diag_field('ocean_model', 'hf_du_dt_visc', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Acceleration from Vertical Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_du_dt_visc > 0) then
  !  call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  !  call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  !endif

  !CS%id_hf_dv_dt_visc = register_diag_field('ocean_model', 'hf_dv_dt_visc', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration from Vertical Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_dv_dt_visc > 0) then
  !  call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  !  call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  !endif

  CS%id_hf_du_dt_visc_2d = register_diag_field('ocean_model', 'hf_du_dt_visc_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Vertical Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_du_dt_visc_2d > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_hf_dv_dt_visc_2d = register_diag_field('ocean_model', 'hf_dv_dt_visc_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Vertical Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_dv_dt_visc_2d > 0) then
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_h_du_dt_visc = register_diag_field('ocean_model', 'h_du_dt_visc', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_du_dt_visc > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hu,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_h_dv_dt_visc = register_diag_field('ocean_model', 'h_dv_dt_visc', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_dv_dt_visc > 0) then
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hv,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_h_du_dt_str = register_diag_field('ocean_model', 'h_du_dt_str', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration from Surface Wind Stresses', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_du_dt_str > 0) then
    call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hu,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_h_dv_dt_str = register_diag_field('ocean_model', 'h_dv_dt_str', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration from Surface Wind Stresses', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_dv_dt_str > 0) then
    call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hv,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_du_dt_str_visc_rem = register_diag_field('ocean_model', 'du_dt_str_visc_rem', diag%axesCuL, Time, &
      'Zonal Acceleration from Surface Wind Stresses multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_str_visc_rem > 0) then
    call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%visc_rem_u,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_dv_dt_str_visc_rem = register_diag_field('ocean_model', 'dv_dt_str_visc_rem', diag%axesCvL, Time, &
      'Meridional Acceleration from Surface Wind Stresses multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_str_visc_rem > 0) then
    call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%visc_rem_v,isd,ied,JsdB,JedB,nz)
  endif

  if ((len_trim(CS%u_trunc_file) > 0) .or. (len_trim(CS%v_trunc_file) > 0)) &
    call PointAccel_init(MIS, Time, G, param_file, diag, dirs, CS%PointAccel_CSp)

  !$omp target enter data map(to: CS)
  !$omp target enter data map(to: CS%a_u, CS%a_v, CS%h_u, CS%h_v)
  !$omp target enter data map(to: CS%a_u_gl90, CS%a_v_gl90)
  if (allocated(CS%kappa_gl90_2d)) then
    !$omp target enter data map(to: CS%kappa_gl90_2d)
  endif

end subroutine vertvisc_init

!> Update the CFL truncation value as a function of time.
!! If called with the optional argument activate=.true., record the
!! value of Time as the beginning of the ramp period.
subroutine updateCFLtruncationValue(Time, CS, US, activate)
  type(time_type), target, intent(in)    :: Time     !< Current model time
  type(vertvisc_CS),       intent(inout) :: CS       !< Vertical viscosity control structure
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  logical, optional,       intent(in)    :: activate !< Specify whether to record the value of
                                                     !! Time as the beginning of the ramp period

  ! Local variables
  real :: deltaTime ! The time since CS%rampStartTime [T ~> s], which may be negative.
  real :: wghtA     ! The relative weight of the final value [nondim]
  character(len=12) :: msg

  if (CS%truncRampTime==0.) return ! This indicates to ramping is turned off

  ! We use the optional argument to indicate this Time should be recorded as the
  ! beginning of the ramp-up period.
  if (present(activate)) then
    if (activate) then
      CS%rampStartTime = Time ! Record the current time
      CS%CFLrampingIsActivated = .true.
    endif
  endif
  if (.not.CS%CFLrampingIsActivated) return
  deltaTime = max(0., time_minus_signed(Time, CS%rampStartTime, scale=US%s_to_T))
  if (deltaTime >= CS%truncRampTime) then
    CS%CFL_trunc = CS%CFL_truncE
    CS%truncRampTime = 0. ! This turns off ramping after this call
  else
    wghtA = min( 1., deltaTime / CS%truncRampTime ) ! Linear profile in time
    !wghtA = wghtA*wghtA ! Convert linear profile to parabolic profile in time
    !wghtA = wghtA*wghtA*(3. - 2.*wghtA) ! Convert linear profile to cosine profile
    wghtA = 1. - ( (1. - wghtA)**2 ) ! Convert linear profile to inverted parabolic profile
    CS%CFL_trunc = CS%CFL_truncS + wghtA * ( CS%CFL_truncE - CS%CFL_truncS )
  endif
  !$omp target update to(CS%CFL_trunc)
  write(msg(1:12),'(es12.3)') CS%CFL_trunc
  call MOM_error(NOTE, "MOM_vert_friction: updateCFLtruncationValue set CFL limit to "//trim(msg))
end subroutine updateCFLtruncationValue

!> Clean up and deallocate the vertical friction module
subroutine vertvisc_end(CS)
  type(vertvisc_CS), intent(inout) :: CS  !< Vertical viscosity control structure that
                                          !! will be deallocated in this subroutine.

  if ((len_trim(CS%u_trunc_file) > 0) .or. (len_trim(CS%v_trunc_file) > 0)) &
    deallocate(CS%PointAccel_CSp)

  !$omp target exit data map(release: CS%a_u, CS%a_v, CS%h_u, CS%h_v)
  !$omp target exit data map(release: CS%a_u_gl90, CS%a_v_gl90)
  if (associated(CS%a1_shelf_u)) then
    !$omp target exit data map(release: CS%a1_shelf_u)
  endif
  if (associated(CS%a1_shelf_v)) then
    !$omp target exit data map(release: CS%a1_shelf_v)
  endif
  if (allocated(CS%kappa_gl90_2d)) then
    !$omp target exit data map(release: CS%kappa_gl90_2d)
  endif
  !$omp target exit data map(release: CS)

  DEALLOC_(CS%a_u) ; DEALLOC_(CS%h_u)
  DEALLOC_(CS%a_v) ; DEALLOC_(CS%h_v)
  if (associated(CS%a1_shelf_u)) deallocate(CS%a1_shelf_u)
  if (associated(CS%a1_shelf_v)) deallocate(CS%a1_shelf_v)
  if (allocated(CS%kappa_gl90_2d)) deallocate(CS%kappa_gl90_2d)
end subroutine vertvisc_end

!> \namespace mom_vert_friction
!! \author Robert Hallberg
!! \date April 1994 - October 2006
!!
!!  The vertical diffusion of momentum is fully implicit.  This is
!!  necessary to allow for vanishingly small layers.  The coupling
!!  is based on the distance between the centers of adjacent layers,
!!  except where a layer is close to the bottom compared with a
!!  bottom boundary layer thickness when a bottom drag law is used.
!!  A stress top b.c. and a no slip bottom  b.c. are used.  There
!!  is no limit on the time step for vertvisc.
!!
!!  Near the bottom, the horizontal thickness interpolation scheme
!!  changes to an upwind biased estimate to control the effect of
!!  spurious Montgomery potential gradients at the bottom where
!!  nearly massless layers layers ride over the topography.  Within a
!!  few boundary layer depths of the bottom, the harmonic mean
!!  thickness (i.e. (2 h+ h-) / (h+ + h-) ) is used if the velocity
!!  is from the thinner side and the arithmetic mean thickness
!!  (i.e. (h+ + h-)/2) is used if the velocity is from the thicker
!!  side.  Both of these thickness estimates are second order
!!  accurate.  Above this the arithmetic mean thickness is used.
!!
!!  In addition, vertvisc truncates any velocity component that exceeds a
!!  maximum CFL number to a fraction of this value.  This basically keeps
!!  instabilities spatially localized.  The number of times the velocity is
!!  truncated is reported each time the energies are saved, and if
!!  exceeds CS%Maxtrunc the model will stop itself and change the time
!!  to a large value.  This has proven very useful in (1) diagnosing
!!  model failures and (2) letting the model settle down to a
!!  meaningful integration from a poorly specified initial condition.
!!
!!  The same code is used for the two velocity components, by
!!  indirectly referencing the velocities and defining a handful of
!!  direction-specific defined variables.
!!
!!  Macros written all in capital letters are defined in MOM_memory.h.
!!
!!     A small fragment of the grid is shown below:
!! \verbatim
!!    j+1  x ^ x ^ x   At x:  q
!!    j+1  > o > o >   At ^:  v, frhatv, tauy
!!    j    x ^ x ^ x   At >:  u, frhatu, taux
!!    j    > o > o >   At o:  h
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!!  The boundaries always run through q grid points (x).
end module MOM_vert_friction
