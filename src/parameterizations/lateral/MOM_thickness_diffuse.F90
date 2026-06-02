! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Isopycnal height diffusion (or Gent McWilliams diffusion)
module MOM_thickness_diffuse

use MOM_debugging,             only : hchksum, uvchksum
use MOM_diag_mediator,         only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,         only : register_diag_field, safe_alloc_ptr, time_type
use MOM_diag_mediator,         only : diag_update_remap_grids
use MOM_domains,               only : pass_var, CORNER, pass_vector
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_EOS,                   only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_EOS,                   only : calculate_density_second_derivs
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_io,                    only : MOM_read_data, slasher
use MOM_interface_heights,     only : find_eta, thickness_to_dz
use MOM_isopycnal_slopes,      only : vert_fill_TS
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_stochastics,           only : stochastic_CS
use MOM_unit_scaling,          only : unit_scale_type
use MOM_variables,             only : thermo_var_ptrs, cont_diag_ptrs
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_meso_sfn_ANN,          only : meso_sfn_ANN_compute, MESO_SFN_ANN_CS
use MOM_meso_sfn_ANN,          only : meso_sfn_ANN_init, meso_sfn_ANN_end

implicit none ; private

#include <MOM_memory.h>

public thickness_diffuse, thickness_diffuse_init, thickness_diffuse_end
public thickness_diffuse_get_KH

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for thickness_diffuse
type, public :: thickness_diffuse_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: Khth                !< Background isopycnal depth diffusivity [L2 T-1 ~> m2 s-1]
  real    :: Khth_Slope_Cff      !< Slope dependence coefficient of Khth [nondim]
  real    :: max_Khth_CFL        !< Maximum value of the diffusive CFL for isopycnal height diffusion [nondim]
  real    :: Khth_Min            !< Minimum value of Khth [L2 T-1 ~> m2 s-1]
  real    :: Khth_Max            !< Maximum value of Khth [L2 T-1 ~> m2 s-1], or 0 for no max
  real    :: Kh_eta_bg           !< Background isopycnal height diffusivity [L2 T-1 ~> m2 s-1]
  real    :: Kh_eta_vel          !< Velocity scale that is multiplied by the grid spacing to give
                                 !! the isopycnal height diffusivity [L T-1 ~> m s-1]
  real    :: slope_max           !< Slopes steeper than slope_max are limited in some way [Z L-1 ~> nondim]
  real    :: kappa_smooth        !< Vertical diffusivity used to interpolate more sensible values
                                 !! of T & S into thin layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  logical :: thickness_diffuse   !< If true, interfaces heights are diffused.
  logical :: full_depth_khth_min !< If true, KHTH_MIN is enforced throughout the whole water column.
                                 !! Otherwise, KHTH_MIN is only enforced at the surface. This parameter
                                 !! is only available when KHTH_USE_EBT_STRUCT=True and KHTH_MIN>0.
  logical :: use_FGNV_streamfn   !< If true, use the streamfunction formulation of
                                 !! Ferrari et al., 2010, which effectively emphasizes
                                 !! graver vertical modes by smoothing in the vertical.
  real    :: FGNV_scale          !< A coefficient scaling the vertical smoothing term in the
                                 !! Ferrari et al., 2010, streamfunction formulation [nondim].
  real    :: FGNV_c_min          !< A minimum wave speed used in the Ferrari et al., 2010,
                                 !! streamfunction formulation [L T-1 ~> m s-1].
  real    :: N2_floor            !< A floor for squared buoyancy frequency in the Ferrari et al., 2010,
                                 !! streamfunction formulation divided by aspect ratio rescaling factors
                                 !! [L2 Z-2 T-2 ~> s-2].
  logical :: detangle_interfaces !< If true, add 3-d structured interface height
                                 !! diffusivities to horizontally smooth jagged layers.
  real    :: detangle_time       !< If detangle_interfaces is true, this is the
                                 !! timescale over which maximally jagged grid-scale
                                 !! thickness variations are suppressed [T ~> s].  This must be
                                 !! longer than DT, or 0 (the default) to use DT.
  integer :: nkml                !< number of layers within mixed layer
  logical :: debug               !< write verbose checksums for debugging purposes
  logical :: use_GME_thickness_diffuse !< If true, passes GM coefficients to MOM_hor_visc for use
                                 !! with GME closure.
  logical :: MEKE_GEOMETRIC      !< If true, uses the GM coefficient formulation from the GEOMETRIC
                                 !! framework (Marshall et al., 2012)
  real    :: MEKE_GEOMETRIC_alpha!< The nondimensional coefficient governing the efficiency of
                                 !! the GEOMETRIC isopycnal height diffusion [nondim]
  real    :: MEKE_GEOMETRIC_epsilon !< Minimum Eady growth rate for the GEOMETRIC thickness
                                 !! diffusivity [T-1 ~> s-1].
  integer :: MEKE_GEOM_answer_date  !< The vintage of the expressions in the MEKE_GEOMETRIC
                                 !! calculation.  Values below 20190101 recover the answers from the
                                 !! original implementation, while higher values use expressions that
                                 !! satisfy rotational symmetry.
  logical :: Use_KH_in_MEKE      !< If true, uses the isopycnal height diffusivity calculated here to diffuse MEKE.
  real    :: MEKE_min_depth_diff !< The minimum total depth over which to average the diffusivity
                                 !! used for MEKE [H ~> m or kg m-2].  When the total depth is less
                                 !! than this, the diffusivity is scaled away.
  logical :: GM_src_alt          !< If true, use the GM energy conversion form S^2*N^2*kappa rather
                                 !! than the streamfunction for the GM source term for MEKE.
  integer :: MEKE_src_answer_date  !< The vintage of the expressions in the GM energy conversion
                                 !! calculation when MEKE_GM_SRC_ALT is true.  Values below 20240601
                                 !! recover the answers from the original implementation, while higher
                                 !! values use expressions that satisfy rotational symmetry.
  logical :: MEKE_src_slope_bug  !< If true, use a bug that limits the positive values, but not the
                                 !! negative values, of the slopes used when MEKE_GM_SRC_ALT is true.
                                 !! When this is true, it breaks rotational symmetry.
  logical :: use_GM_work_bug     !< If true, use the incorrect sign for the
                                 !! top-level work tendency on the top layer.
  logical :: read_khth           !< If true, read a file containing the spatially varying horizontal
                                 !! isopycnal height diffusivity
  logical :: use_stanley_gm      !< If true, also use the Stanley parameterization in MOM_thickness_diffuse

  logical :: use_meso_sfn_ANN  !< If true, use the meso-scale streamfunction ANN parameterization
  type(MESO_SFN_ANN_CS) :: meso_sfn_ANN_CS !< Control structure for the meso-scale streamfunction ANN parameterization

  type(diag_ctrl), pointer :: diag => NULL() !< structure used to regulate timing of diagnostics
  real, allocatable :: GMwork(:,:)        !< Work by isopycnal height diffusion [R Z L2 T-3 ~> W m-2]
  real, allocatable :: diagSlopeX(:,:,:)  !< Diagnostic: zonal neutral slope [Z L-1 ~> nondim]
  real, allocatable :: diagSlopeY(:,:,:)  !< Diagnostic: zonal neutral slope [Z L-1 ~> nondim]

  real, allocatable :: Kh_eta_u(:,:)    !< Isopycnal height diffusivities at u points [L2 T-1 ~> m2 s-1]
  real, allocatable :: Kh_eta_v(:,:)    !< Isopycnal height diffusivities in v points [L2 T-1 ~> m2 s-1]

  real, allocatable :: KH_u_GME(:,:,:)  !< Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
  real, allocatable :: KH_v_GME(:,:,:)  !< Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
  real, allocatable :: khth2d(:,:)      !< 2D isopycnal height diffusivity at h-points [L2 T-1 ~> m2 s-1]

  !>@{
  !! Diagnostic identifier
  integer :: id_uhGM    = -1, id_vhGM    = -1, id_GMwork = -1
  integer :: id_KH_u    = -1, id_KH_v    = -1, id_KH_t   = -1
  integer :: id_KH_u1   = -1, id_KH_v1   = -1, id_KH_t1  = -1
  integer :: id_slope_x = -1, id_slope_y = -1
  integer :: id_sfn_unlim_x = -1, id_sfn_unlim_y = -1, id_sfn_x = -1, id_sfn_y = -1
  !>@}
end type thickness_diffuse_CS

interface

!> Calculates isopycnal height diffusion coefficients and applies isopycnal height diffusion
!! by modifying to the layer thicknesses, h. Diffusivities are limited to ensure stability.
!! Also returns along-layer mass fluxes used in the continuity equation.
module subroutine thickness_diffuse(h, uhtr, vhtr, tv, dt, G, GV, US, MEKE, VarMix, CDp, CS, STOCH, u, v)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !! [L2 H ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !! [L2 H ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamics structure
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  type(MEKE_type),                            intent(inout) :: MEKE   !< MEKE fields
  type(VarMix_CS), target,                    intent(in)    :: VarMix !< Variable mixing coefficients
  type(cont_diag_ptrs),                       intent(inout) :: CDp    !< Diagnostics for the continuity equation
  type(thickness_diffuse_CS),                 intent(inout) :: CS     !< Control structure for thickness_diffuse
  type(stochastic_CS),                        intent(inout) :: STOCH !< Stochastic control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v !< Meridional velocity [L T-1 ~> m s-1].

end subroutine thickness_diffuse

!> Calculates parameterized layer transports for use in the continuity equation.
!! Fluxes are limited to give positive definite thicknesses.
!! Called by thickness_diffuse().
module subroutine thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, &
                                  CS, int_slope_u, int_slope_v, slope_x, slope_y, STOCH, VarMix, &
                                  Sfn_unlim_u_3D, Sfn_unlim_v_3D)
  type(ocean_grid_type),                        intent(in)  :: G     !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)  :: GV    !< Vertical grid structure
  type(unit_scale_type),                        intent(in)  :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)  :: h     !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)  :: e     !< Interface positions [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: Kh_u  !< Isopycnal height diffusivity
                                                                     !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)  :: Kh_v  !< Isopycnal height diffusivity
                                                                     !! at v points [L2 T-1 ~> m2 s-1]
  type(thermo_var_ptrs),                        intent(in)  :: tv    !< Thermodynamics structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(out) :: uhD   !< Zonal mass fluxes
                                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(out) :: vhD   !< Meridional mass fluxes
                                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(:,:),                         pointer     :: cg1   !< Wave speed [L T-1 ~> m s-1]
  real,                                         intent(in)  :: dt    !< Time increment [T ~> s]
  type(MEKE_type),                              intent(inout) :: MEKE !< MEKE fields
  type(thickness_diffuse_CS),                   intent(inout) :: CS  !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(in)  :: int_slope_u !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from
                                                                     !! the interface slopes without consideration of
                                                                     !! density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(in)  :: int_slope_v !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from
                                                                     !! the interface slopes without consideration of
                                                                     !! density gradients [nondim].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), optional, intent(in)  :: slope_x !< Isopyc. slope at u [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), optional, intent(in)  :: slope_y !< Isopyc. slope at v [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), optional, intent(in) :: Sfn_unlim_u_3D !< ANN streamfunction
                                                                      !! at u [Z L2 T-1 ~> m3 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), optional, intent(in) :: Sfn_unlim_v_3D !< ANN streamfunction
                                                                      !! at v [Z L2 T-1 ~> m3 s-1]
  type(stochastic_CS),                       optional, intent(inout)  :: STOCH !< Stochastic control structure
  type(VarMix_CS), target,                      optional, intent(in)  :: VarMix !< Variable mixing coefficents

end subroutine thickness_diffuse_full
end interface

contains

!> Tridiagonal solver for streamfunction at interfaces
subroutine streamfn_solver(nk, c2_h, hN2, sfn)
  integer,               intent(in)    :: nk   !< Number of layers
  real, dimension(nk),   intent(in)    :: c2_h !< Wave speed squared over thickness in layers, rescaled to
                                               !! [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2]
  real, dimension(nk+1), intent(in)    :: hN2  !< Thickness times N2 at interfaces times rescaling factors
                                               !! [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2]
  real, dimension(nk+1), intent(inout) :: sfn  !< Streamfunction [H L2 T-1 ~> m3 s-1 or kg s-1] or arbitrary units
                                               !! On entry, equals diffusivity times slope.
                                               !! On exit, equals the streamfunction.
  ! Local variables
  real :: c1(nk)  ! The dependence of the final streamfunction on the values below [nondim]
  real :: d1      ! The complement of c1(k) (i.e., 1 - c1(k)) [nondim]
  real :: b_denom ! A term in the denominator of beta [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2]
  real :: beta    ! The normalization for the pivot [Z2 T2 H-1 L-2 ~> s2 m-1 or m2 s2 kg-1]
  integer :: k

  sfn(1) = 0.
  b_denom = hN2(2) + c2_h(1)
  beta = 1.0 / ( b_denom + c2_h(2) )
  d1 = beta * b_denom
  sfn(2) = ( beta * hN2(2) )*sfn(2)
  do K=3,nk
    c1(k-1) = beta * c2_h(k-1)
    b_denom = hN2(K) + d1*c2_h(k-1)
    beta = 1.0 / (b_denom + c2_h(k))
    d1 = beta * b_denom
    sfn(K) = beta * (hN2(K)*sfn(K) + c2_h(k-1)*sfn(K-1))
  enddo
  c1(nk) = beta * c2_h(nk)
  sfn(nk+1) = 0.
  do K=nk,2,-1
    sfn(K) = sfn(K) + c1(k)*sfn(K+1)
  enddo

end subroutine streamfn_solver

!> Add a diffusivity that acts on the isopycnal heights, regardless of the densities
subroutine add_interface_Kh(G, GV, US, CS, Kh_u, Kh_v, Kh_u_CFL, Kh_v_CFL, int_slope_u, int_slope_v)
  type(ocean_grid_type),                        intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US   !< A dimensional unit scaling type
  type(thickness_diffuse_CS),                   intent(in)    :: CS   !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kh_u !< Isopycnal height diffusivity
                                                                      !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: Kh_v !< Isopycnal height diffusivity
                                                                      !! at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),            intent(in)    :: Kh_u_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G)),            intent(in)    :: Kh_v_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: int_slope_u !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: int_slope_v !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].

  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do k=1,nz+1 ; do j=js,je ; do I=is-1,ie ; if (CS%Kh_eta_u(I,j) > 0.0) then
    int_slope_u(I,j,K) = (int_slope_u(I,j,K)*Kh_u(I,j,K) + CS%Kh_eta_u(I,j)) / &
                         (Kh_u(I,j,K) + CS%Kh_eta_u(I,j))
    Kh_u(I,j,K) = min(Kh_u(I,j,K) + CS%Kh_eta_u(I,j), Kh_u_CFL(I,j))
  endif ; enddo ; enddo ; enddo

  do k=1,nz+1 ; do J=js-1,je ; do i=is,ie ; if (CS%Kh_eta_v(i,J) > 0.0) then
    int_slope_v(i,J,K) = (int_slope_v(i,J,K)*Kh_v(i,J,K) + CS%Kh_eta_v(i,J)) / &
                         (Kh_v(i,J,K) + CS%Kh_eta_v(i,J))
    Kh_v(i,J,K) = min(Kh_v(i,J,K) + CS%Kh_eta_v(i,J), Kh_v_CFL(i,J))
  endif ; enddo ; enddo ; enddo

end subroutine add_interface_Kh

!> Modifies isopycnal height diffusivities to untangle layer structures
subroutine add_detangling_Kh(h, e, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, tv, dt, G, GV, US, CS, &
                             int_slope_u, int_slope_v)
  type(ocean_grid_type),                        intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV   !< Vertical grid structure
  type(unit_scale_type),                        intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h    !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1),  intent(in)    :: e    !< Interface positions [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: Kh_u !< Isopycnal height diffusivity
                                                                      !! at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: Kh_v !< Isopycnal height diffusivity
                                                                      !! at v points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G)),            intent(in)    :: Kh_u_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at u points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G)),            intent(in)    :: Kh_v_CFL !< Maximum stable isopycnal height
                                                                      !! diffusivity at v points [L2 T-1 ~> m2 s-1]
  type(thermo_var_ptrs),                        intent(in)    :: tv   !< Thermodynamics structure
  real,                                         intent(in)    :: dt   !< Time increment [T ~> s]
  type(thickness_diffuse_CS),                   intent(in)    :: CS   !< Control structure for thickness_diffuse
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: int_slope_u !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: int_slope_v !< Ratio that determine how much of
                                                                      !! the isopycnal slopes are taken directly from
                                                                      !! the interface slopes without consideration
                                                                      !! of density gradients [nondim].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    de_top     ! The distances between the top of a layer and the top of the
               ! region where the detangling is applied [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    Kh_lay_u   ! The tentative isopycnal height diffusivity for each layer at
               ! u points [L2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    Kh_lay_v   ! The tentative isopycnal height diffusivity for each layer at
               ! v points [L2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    de_bot     ! The distances from the bottom of the region where the
               ! detangling is applied [H ~> m or kg m-2].
  real :: h1, h2    ! The thinner and thicker surrounding thicknesses [H ~> m or kg m-2],
                    ! with the thinner modified near the boundaries to mask out
                    ! thickness variations due to topography, etc.
  real :: jag_Rat   ! The nondimensional jaggedness ratio for a layer, going
                    ! from 0 (smooth) to 1 (jagged) [nondim].  This is the difference
                    ! between the arithmetic and harmonic mean thicknesses
                    ! normalized by the arithmetic mean thickness.
  real :: Kh_scale  ! A ratio by which Kh_u_CFL is scaled for maximally jagged
                    ! layers [nondim].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: I_sl      ! The absolute value of the larger in magnitude of the slopes
                    ! above and below [L Z-1 ~> nondim].
  real :: Rsl       ! The ratio of the smaller magnitude slope to the larger
                    ! magnitude one [nondim]. 0 <= Rsl <1.
  real :: IRsl      ! The (limited) inverse of Rsl [nondim]. 1 < IRsl <= 1e9.
  real :: dH        ! The thickness gradient divided by the damping timescale
                    ! and the ratio of the face length to the adjacent cell
                    ! areas for comparability with the diffusivities [L Z T-1 ~> m2 s-1].
  real :: adH       ! The absolute value of dH [L Z T-1 ~> m2 s-1].
  real :: sign      ! 1 or -1, with the same sign as the layer thickness gradient [nondim].
  real :: sl_K      ! The sign-corrected slope of the interface above [Z L-1 ~> nondim].
  real :: sl_Kp1    ! The sign-corrected slope of the interface below [Z L-1 ~> nondim].
  real :: I_sl_K    ! The (limited) inverse of sl_K [L Z-1 ~> nondim].
  real :: I_sl_Kp1  ! The (limited) inverse of sl_Kp1 [L Z-1 ~> nondim].
  real :: I_4t      ! A quarter of a flux scaling factor divided by
                    ! the damping timescale [T-1 ~> s-1].
  real :: Fn_R      ! A function of Rsl, such that Rsl < Fn_R < 1 [nondim]
  real :: Idx_eff   ! The effective inverse x-grid spacing at a u-point [L-1 ~> m-1]
  real :: Idy_eff   ! The effective inverse y-grid spacing at a v-point [L-1 ~> m-1]
  real :: slope_sq  ! The sum of the squared slopes above and below a layer [Z2 L-2 ~> nondim]
  real :: Kh_max    ! A local ceiling on the diffusivity [L2 T-1 ~> m2 s-1].
  real :: wt1, wt2  ! Nondimensional weights [nondim].
  !   Variables used only in testing code.
  ! real, dimension(SZK_(GV)) :: uh_here ! The transport in a layer [Z L2 T-1 ~> m3 s-1]
  ! real, dimension(SZK_(GV)+1) :: Sfn ! The streamfunction at an interface [Z L T-1 ~> m2 s-1]
  real :: dKh       ! An increment in the diffusivity [L2 T-1 ~> m2 s-1].

  real, dimension(SZIB_(G),SZK_(GV)+1) :: &
    Kh_bg, &        ! The background (floor) value of Kh [L2 T-1 ~> m2 s-1].
    Kh, &           ! The tentative value of Kh [L2 T-1 ~> m2 s-1].
    Kh_detangle, &  ! The detangling diffusivity that could be used [L2 T-1 ~> m2 s-1].
    Kh_min_max_p, & ! The smallest ceiling that can be placed on Kh(I,K)
                    ! based on the value of Kh(I,K+1) [L2 T-1 ~> m2 s-1].
    Kh_min_max_m, & ! The smallest ceiling that can be placed on Kh(I,K)
                    ! based on the value of Kh(I,K-1) [L2 T-1 ~> m2 s-1].
    ! The following are variables that define the relationships between
    ! successive values of Kh.
    ! Search for Kh that satisfy...
    !    Kh(I,K) >= Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K)
    !    Kh(I,K) >= Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K)
    !    Kh(I,K) <= Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K)
    !    Kh(I,K) <= Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K)
    Kh_min_m , &   ! See above [nondim].
    Kh0_min_m , &  ! See above [L2 T-1 ~> m2 s-1].
    Kh_max_m , &   ! See above [nondim].
    Kh0_max_m, &   ! See above [L2 T-1 ~> m2 s-1].
    Kh_min_p , &   ! See above [nondim].
    Kh0_min_p , &  ! See above [L2 T-1 ~> m2 s-1].
    Kh_max_p , &   ! See above [nondim].
    Kh0_max_p      ! See above [L2 T-1 ~> m2 s-1].
  real, dimension(SZIB_(G)) :: &
    Kh_max_max  ! The maximum diffusivity permitted in a column [L2 T-1 ~> m2 s-1]
  logical, dimension(SZIB_(G)) :: &
    do_i        ! If true, work on a column.
  integer :: i, j, k, n, ish, jsh, is, ie, js, je, nz, k_top
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  k_top = GV%nk_rho_varies + 1
  h_neglect = GV%H_subroundoff
  !   The 0.5 is because we are not using uniform weightings, but are
  ! distributing the diffusivities more effectively (with wt1 & wt2), but this
  ! means that the additions to a single interface can be up to twice as large.
  Kh_scale = 0.5
  if (CS%detangle_time > dt) Kh_scale = 0.5 * dt / CS%detangle_time

  do j=js-1,je+1 ; do i=is-1,ie+1
    de_top(i,j,k_top) = 0.0 ; de_bot(i,j) = 0.0
  enddo ; enddo
  do k=k_top+1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
    de_top(i,j,k) = de_top(i,j,k-1) + h(i,j,k-1)
  enddo ; enddo ; enddo

  do j=js,je ; do I=is-1,ie
    Kh_lay_u(I,j,nz) = 0.0 ; Kh_lay_u(I,j,k_top) = 0.0
  enddo ; enddo
  do J=js-1,je ; do i=is,ie
    Kh_lay_v(i,J,nz) = 0.0 ; Kh_lay_v(i,J,k_top) = 0.0
  enddo ; enddo

  do k=nz-1,k_top+1,-1
    ! Find the diffusivities associated with each layer.
    do j=js-1,je+1 ; do i=is-1,ie+1
      de_bot(i,j) = de_bot(i,j) + h(i,j,k+1)
    enddo ; enddo

    do j=js,je ; do I=is-1,ie ; if (G%OBCmaskCu(I,j) > 0.0) then
      if (h(i,j,k) > h(i+1,j,k)) then
        h2 = h(i,j,k)
        h1 = max( h(i+1,j,k), h2 - min(de_bot(i+1,j), de_top(i+1,j,k)) )
      else
        h2 = h(i+1,j,k)
        h1 = max( h(i,j,k), h2 - min(de_bot(i,j), de_top(i,j,k)) )
      endif
      jag_Rat = (h2 - h1)**2 / (h2 + h1 + h_neglect)**2
      KH_lay_u(I,j,k) = (Kh_scale * KH_u_CFL(I,j)) * jag_Rat**2
    endif ; enddo ; enddo

    do J=js-1,je ; do i=is,ie ; if (G%OBCmaskCv(i,J) > 0.0) then
      if (h(i,j,k) > h(i,j+1,k)) then
        h2 = h(i,j,k)
        h1 = max( h(i,j+1,k), h2 - min(de_bot(i,j+1), de_top(i,j+1,k)) )
      else
        h2 = h(i,j+1,k)
        h1 = max( h(i,j,k), h2 - min(de_bot(i,j), de_top(i,j,k)) )
      endif
      jag_Rat = (h2 - h1)**2 / (h2 + h1 + h_neglect)**2
      KH_lay_v(i,J,k) = (Kh_scale * KH_v_CFL(i,J)) * jag_Rat**2
    endif ; enddo ; enddo
  enddo

  ! Limit the diffusivities

  I_4t = Kh_scale / (4.0 * dt)

  do n=1,2
    if (n==1) then ; jsh = js ; ish = is-1
    else ; jsh = js-1 ; ish = is ; endif

    do j=jsh,je

      ! First, populate the diffusivities
      if (n==1) then ! This is a u-column.
        do i=ish,ie
          do_i(I) = (G%OBCmaskCu(I,j) > 0.0)
          Kh_Max_max(I) = KH_u_CFL(I,j)
        enddo
        do K=1,nz+1 ; do i=ish,ie
          Kh_bg(I,K) = KH_u(I,j,K) ; Kh(I,K) = Kh_bg(I,K)
          Kh_min_max_p(I,K) = Kh_bg(I,K) ; Kh_min_max_m(I,K) = Kh_bg(I,K)
          Kh_detangle(I,K) = 0.0
        enddo ; enddo
      else ! This is a v-column.
        do i=ish,ie
          do_i(i) = (G%OBCmaskCv(i,J) > 0.0) ; Kh_Max_max(I) = KH_v_CFL(i,J)
        enddo
        do K=1,nz+1 ; do i=ish,ie
          Kh_bg(I,K) = KH_v(I,j,K) ; Kh(I,K) = Kh_bg(I,K)
          Kh_min_max_p(I,K) = Kh_bg(I,K) ; Kh_min_max_m(I,K) = Kh_bg(I,K)
          Kh_detangle(I,K) = 0.0
        enddo ; enddo
      endif

      ! Determine the limits on the diffusivities.
      do k=k_top,nz ; do i=ish,ie ; if (do_i(i)) then
        if (n==1) then ! This is a u-column.
          dH = 0.0
          Idx_eff = ((G%IareaT(i+1,j) + G%IareaT(i,j)) * G%dy_Cu(I,j))
          !   This expression uses differences in e in place of h for better
          ! consistency with the slopes.
          if (Idx_eff > 0.0) &
            dH = I_4t * ((e(i+1,j,K) - e(i+1,j,K+1)) - &
                         (e(i,j,K) - e(i,j,K+1))) / Idx_eff
           ! dH = I_4t * (h(i+1,j,k) - h(i,j,k)) / Idx_eff

          adH = abs(dH)
          sign = 1.0 ; if (dH < 0) sign = -1.0
          sl_K = sign * (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
          sl_Kp1 = sign * (e(i+1,j,K+1)-e(i,j,K+1)) * G%IdxCu(I,j)

          ! Add the incremental diffusivities to the surrounding interfaces.
          ! Adding more to the more steeply sloping layers (as below) makes
          ! the diffusivities more than twice as effective.
          slope_sq = (sl_K**2 + sl_Kp1**2)
          wt1 = 0.5 ; wt2 = 0.5
          if (slope_sq > 0.0) then
            wt1 = sl_K**2 / slope_sq ; wt2 = sl_Kp1**2 / slope_sq
          endif
          Kh_detangle(I,K) = Kh_detangle(I,K) + wt1*KH_lay_u(I,j,k)
          Kh_detangle(I,K+1) = Kh_detangle(I,K+1) + wt2*KH_lay_u(I,j,k)
        else ! This is a v-column.
          dH = 0.0
          Idy_eff = ((G%IareaT(i,j+1) + G%IareaT(i,j)) * G%dx_Cv(I,j))
          if (Idy_eff > 0.0) &
            dH = I_4t * ((e(i,j+1,K) - e(i,j+1,K+1)) - &
                         (e(i,j,K) - e(i,j,K+1))) / Idy_eff
           ! dH = I_4t * (h(i,j+1,k) - h(i,j,k)) / Idy_eff

          adH = abs(dH)
          sign = 1.0 ; if (dH < 0) sign = -1.0
          sl_K = sign * (e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J)
          sl_Kp1 = sign * (e(i,j+1,K+1)-e(i,j,K+1)) * G%IdyCv(i,J)

          ! Add the incremental diffusivities to the surrounding interfaces.
          ! Adding more to the more steeply sloping layers (as below) makes
          ! the diffusivities more than twice as effective.
          slope_sq = (sl_K**2 + sl_Kp1**2)
          wt1 = 0.5 ; wt2 = 0.5
          if (slope_sq > 0.0) then
            wt1 = sl_K**2 / slope_sq ; wt2 = sl_Kp1**2 / slope_sq
          endif
          Kh_detangle(I,K) = Kh_detangle(I,K) + wt1*KH_lay_v(i,J,k)
          Kh_detangle(I,K+1) = Kh_detangle(I,K+1) + wt2*KH_lay_v(i,J,k)
        endif

        if (adH == 0.0) then
          Kh_min_m(I,K+1) = 1.0 ; Kh0_min_m(I,K+1) = 0.0
          Kh_max_m(I,K+1) = 1.0 ; Kh0_max_m(I,K+1) = 0.0
          Kh_min_p(I,K) = 1.0 ; Kh0_min_p(I,K) = 0.0
          Kh_max_p(I,K) = 1.0 ; Kh0_max_p(I,K) = 0.0
        elseif (adH > 0.0) then
          if (sl_K <= sl_Kp1) then
            ! This case should only arise from nonlinearities in the equation of state.
            ! Treat it as though dedx(K) = dedx(K+1) & dH = 0.
            Kh_min_m(I,K+1) = 1.0 ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = 1.0 ; Kh0_max_m(I,K+1) = 0.0
            Kh_min_p(I,K) = 1.0 ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = 1.0 ; Kh0_max_p(I,K) = 0.0
          elseif (sl_K <= 0.0) then   ! Both slopes are opposite to dH
            I_sl = -1.0 / sl_Kp1
            Rsl = -sl_K * I_sl                            ! 0 <= Rsl < 1
            IRsl = 1e9 ; if (Rsl > 1e-9) IRsl = 1.0/Rsl   ! 1 < IRsl <= 1e9

            Fn_R = Rsl
            if (Kh_max_max(I) > 0) &
              Fn_R = min(sqrt(Rsl), Rsl + (adH * I_sl) / (Kh_Max_max(I)))

            Kh_min_m(I,K+1) = Fn_R ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = Rsl ; Kh0_max_m(I,K+1) = adH * I_sl
            Kh_min_p(I,K) = IRsl ; Kh0_min_p(I,K) = -adH * (I_sl*IRsl)
            Kh_max_p(I,K) = 1.0/(Fn_R + 1.0e-30) ; Kh0_max_p(I,K) = 0.0
          elseif (sl_Kp1 < 0.0) then  ! Opposite (nonzero) signs of slopes.
            I_sl_K = 1e18*US%Z_to_L ; if (sl_K > 1e-18*US%L_to_Z) I_sl_K = 1.0 / sl_K
            I_sl_Kp1 = 1e18*US%Z_to_L ; if (-sl_Kp1 > 1e-18*US%L_to_Z) I_sl_Kp1 = -1.0 / sl_Kp1

            Kh_min_m(I,K+1) = 0.0 ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = - sl_K*I_sl_Kp1 ; Kh0_max_m(I,K+1) = adH*I_sl_Kp1
            Kh_min_p(I,K) = 0.0 ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = sl_Kp1*I_sl_K ; Kh0_max_p(I,K) = adH*I_sl_K

            ! This limit does not use the slope weighting so that potentially
            ! sharp gradients in diffusivities are not forced to occur.
            Kh_Max = adH / (sl_K - sl_Kp1)
            Kh_min_max_p(I,K) = max(Kh_min_max_p(I,K), Kh_Max)
            Kh_min_max_m(I,K+1) = max(Kh_min_max_m(I,K+1), Kh_Max)
          else ! Both slopes are of the same sign as dH.
            I_sl = 1.0 / sl_K
            Rsl = sl_Kp1 * I_sl                           ! 0 <= Rsl < 1
            IRsl = 1e9 ; if (Rsl > 1e-9) IRsl = 1.0/Rsl   ! 1 < IRsl <= 1e9

            ! Rsl <= Fn_R <= 1
            Fn_R = Rsl
            if (Kh_max_max(I) > 0) &
              Fn_R = min(sqrt(Rsl), Rsl + (adH * I_sl) / Kh_Max_max(I))

            Kh_min_m(I,K+1) = IRsl ; Kh0_min_m(I,K+1) = -adH * (I_sl*IRsl)
            Kh_max_m(I,K+1) = 1.0/(Fn_R + 1.0e-30) ; Kh0_max_m(I,K+1) = 0.0
            Kh_min_p(I,K) = Fn_R ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = Rsl ; Kh0_max_p(I,K) = adH * I_sl
          endif
        endif
      endif ; enddo ; enddo ! I-loop & k-loop

      do k=k_top,nz+1,nz+1-k_top ; do i=ish,ie ; if (do_i(i)) then
        ! The diffusivities at k_top and nz+1 are both fixed.
        Kh_min_m(I,k) = 0.0 ; Kh0_min_m(I,k) = 0.0
        Kh_max_m(I,k) = 0.0 ; Kh0_max_m(I,k) = 0.0
        Kh_min_p(I,k) = 0.0 ; Kh0_min_p(I,k) = 0.0
        Kh_max_p(I,k) = 0.0 ; Kh0_max_p(I,k) = 0.0
        Kh_min_max_p(I,K) = Kh_bg(I,K)
        Kh_min_max_m(I,K) = Kh_bg(I,K)
      endif ; enddo ; enddo ! I-loop and k_top/nz+1 loop

      ! Search for Kh that satisfy...
      !    Kh(I,K) >= Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K)
      !    Kh(I,K) >= Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K)
      !    Kh(I,K) <= Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K)
      !    Kh(I,K) <= Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K)

      ! Increase the diffusivities to satisfy the min constraints.
      ! All non-zero min constraints on one diffusivity are max constraints on another.
      do K=k_top+1,nz ; do i=ish,ie ; if (do_i(i)) then
        Kh(I,K) = max(Kh_bg(I,K), Kh_detangle(I,K), &
                      min(Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K), Kh(I,K-1)))

        if (Kh0_max_m(I,K) > Kh_bg(I,K)) Kh(I,K) = min(Kh(I,K), Kh0_max_m(I,K))
        if (Kh0_max_p(I,K) > Kh_bg(I,K)) Kh(I,K) = min(Kh(I,K), Kh0_max_p(I,K))
      endif ; enddo ; enddo ! I-loop & k-loop
      ! This is still true... do i=ish,ie ; Kh(I,nz+1) = Kh_bg(I,nz+1) ; enddo
      do K=nz,k_top+1,-1 ; do i=ish,ie ; if (do_i(i)) then
        Kh(I,k) = max(Kh(I,K), min(Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K), Kh(I,K+1)))

        Kh_Max = max(Kh_min_max_p(I,K), Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K))
        Kh(I,k) = min(Kh(I,k), Kh_Max)
      endif ; enddo ; enddo ! I-loop & k-loop
      !  All non-zero min constraints on one diffusivity are max constraints on
      ! another layer, so the min constraints can now be discounted.

      ! Decrease the diffusivities to satisfy the max constraints.
        do K=k_top+1,nz ; do i=ish,ie ; if (do_i(i)) then
          Kh_Max = max(Kh_min_max_m(I,K), Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K))
          if (Kh(I,k) > Kh_Max) Kh(I,k) = Kh_Max
        endif ; enddo ; enddo  ! i- and K-loops

      ! This code tests the solutions...
!     do i=ish,ie
!       Sfn(:) = 0.0 ; uh_here(:) = 0.0
!       do K=k_top,nz
!         if ((Kh(i,K) > Kh_bg(i,K)) .or. (Kh(i,K+1) > Kh_bg(i,K+1))) then
!           if (n==1) then ! u-point.
!             if ((h(i+1,j,k) - h(i,j,k)) * &
!                 ((e(i+1,j,K)-e(i+1,j,K+1)) - (e(i,j,K)-e(i,j,K+1))) > 0.0) then
!               Sfn(K) = -Kh(i,K) * (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
!               Sfn(K+1) = -Kh(i,K+1) * (e(i+1,j,K+1)-e(i,j,K+1)) * G%IdxCu(I,j)
!               uh_here(k) = (Sfn(K) - Sfn(K+1))*G%dy_Cu(I,j)
!               if (abs(uh_here(k)) * min(G%IareaT(i,j), G%IareaT(i+1,j)) > &
!                   (1e-10*GV%m_to_H)) then
!                 if (uh_here(k) * (h(i+1,j,k) - h(i,j,k)) > 0.0) then
!                   call MOM_error(WARNING, "Corrective u-transport is up the thickness gradient.", .true.)
!                 endif
!                 if (((h(i,j,k) - 4.0*dt*G%IareaT(i,j)*uh_here(k)) - &
!                      (h(i+1,j,k) + 4.0*dt*G%IareaT(i+1,j)*uh_here(k))) * &
!                     (h(i,j,k) - h(i+1,j,k)) < 0.0) then
!                   call MOM_error(WARNING, "Corrective u-transport is too large.", .true.)
!                 endif
!               endif
!             endif
!           else ! v-point
!             if ((h(i,j+1,k) - h(i,j,k)) * &
!                 ((e(i,j+1,K)-e(i,j+1,K+1)) - (e(i,j,K)-e(i,j,K+1))) > 0.0) then
!               Sfn(K) = -Kh(i,K) * (e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J)
!               Sfn(K+1) = -Kh(i,K+1) * (e(i,j+1,K+1)-e(i,j,K+1)) * G%IdyCv(i,J)
!               uh_here(k) = (Sfn(K) - Sfn(K+1))*G%dx_Cv(i,J)
!               if (abs(uh_here(K)) * min(G%IareaT(i,j), G%IareaT(i,j+1)) > &
!                   (1e-10*GV%m_to_H)) then
!                 if (uh_here(K) * (h(i,j+1,k) - h(i,j,k)) > 0.0) then
!                   call MOM_error(WARNING, &
!                          "Corrective v-transport is up the thickness gradient.", .true.)
!                 endif
!                 if (((h(i,j,k) - 4.0*dt*G%IareaT(i,j)*uh_here(K)) - &
!                      (h(i,j+1,k) + 4.0*dt*G%IareaT(i,j+1)*uh_here(K))) * &
!                     (h(i,j,k) - h(i,j+1,k)) < 0.0) then
!                   call MOM_error(WARNING, &
!                          "Corrective v-transport is too large.", .true.)
!                 endif
!               endif
!             endif
!           endif ! u- or v- selection.
!          !  de_dx(I,K) = (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
!         endif
!       enddo
!     enddo

      if (n==1) then ! This is a u-column.
        do K=k_top+1,nz ; do i=ish,ie
          if (Kh(I,K) > KH_u(I,j,K)) then
            dKh = (Kh(I,K) - KH_u(I,j,K))
            int_slope_u(I,j,K) = dKh / Kh(I,K)
            KH_u(I,j,K) = Kh(I,K)
          endif
        enddo ; enddo
      else ! This is a v-column.
        do K=k_top+1,nz ; do i=ish,ie
          if (Kh(i,K) > KH_v(i,J,K)) then
            dKh = Kh(i,K) - KH_v(i,J,K)
            int_slope_v(i,J,K) = dKh / Kh(i,K)
            KH_v(i,J,K) = Kh(i,K)
          endif
        enddo ; enddo
      endif

    enddo ! j-loop
  enddo  ! n-loop over u- and v- directions.

end subroutine add_detangling_Kh

!> Initialize the isopycnal height diffusion module and its control structure
subroutine thickness_diffuse_init(Time, G, GV, US, param_file, diag, CDp, CS)
  type(time_type),         intent(in) :: Time    !< Current model time
  type(ocean_grid_type),   intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV      !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(cont_diag_ptrs),    intent(inout) :: CDp  !< Continuity equation diagnostics
  type(thickness_diffuse_CS), intent(inout) :: CS !< Control structure for thickness_diffuse

  ! Local variables
  character(len=40)  :: mdl = "MOM_thickness_diffuse" ! This module's name.
  character(len=200) :: khth_file, inputdir, khth_varname
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  real :: grid_sp      ! The local grid spacing [L ~> m]
  real :: omega        ! The Earth's rotation rate [T-1 ~> s-1]
  real :: strat_floor  ! A floor for buoyancy frequency in the Ferrari et al. 2010,
                       ! streamfunction formulation, expressed as a fraction of planetary
                       ! rotation divided by an aspect ratio rescaling factor [L Z-1 ~> nondim]
  real :: Stanley_coeff ! Coefficient relating the temperature gradient and sub-gridscale
                        ! temperature variance [nondim]
  logical :: khth_use_ebt_struct ! If true, uses the equivalent barotropic structure
                                 ! as the vertical structure of thickness diffusivity.
                                 ! Used to determine if FULL_DEPTH_KHTH_MIN should be
                                 ! available.
  logical :: use_meke = .false. ! If true, use the MEKE formulation for the thickness diffusivity.
  integer :: default_answer_date ! The default setting for the various ANSWER_DATE flags.
  logical :: stoch_eos           ! Can't use Stanley param here unless stoch_eos is true
  integer :: i, j

  CS%initialized = .true.
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "THICKNESSDIFFUSE", CS%thickness_diffuse, &
                 "If true, interface heights are diffused with a "//&
                 "coefficient of KHTH.", default=.false.)
  call get_param(param_file, mdl, "USE_THICKNESS_DIFFUSE_ANN", CS%use_meso_sfn_ANN, &
                 "If true, use the ANN to compute the mesoscale streamfunction "//&
                 "for thickness diffusivity.", default=.false.)
  if (CS%use_meso_sfn_ANN) then
    call meso_sfn_ANN_init(Time, G, GV, US, param_file, diag, CS%meso_sfn_ANN_CS)
  endif
  call get_param(param_file, mdl, "KHTH", CS%Khth, &
                 "The background horizontal thickness diffusivity.", &
                 default=0.0, units="m2 s-1", scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "READ_KHTH", CS%read_khth, &
                 "If true, read a file (given by KHTH_FILE) containing the "//&
                 "spatially varying horizontal isopycnal height diffusivity.", &
                 default=.false.)
  if (CS%read_khth) then
    if (CS%Khth > 0) then
        call MOM_error(FATAL, "thickness_diffuse_init: KHTH > 0 is not "// &
              "compatible with READ_KHTH = TRUE. ")
    endif
    call get_param(param_file, mdl, "INPUTDIR", inputdir, &
                 "The directory in which all input files are found.", &
                 default=".", do_not_log=.true.)
    inputdir = slasher(inputdir)
    call get_param(param_file, mdl, "KHTH_FILE", khth_file, &
                 "The file containing the spatially varying horizontal "//&
                 "isopycnal height diffusivity.", default="khth.nc")
    call get_param(param_file, mdl, "KHTH_VARIABLE", khth_varname, &
                 "The name of the isopycnal height diffusivity variable to read "//&
                 "from KHTH_FILE.", &
                 default="khth")
    khth_file = trim(inputdir) // trim(khth_file)

    allocate(CS%khth2d(G%isd:G%ied, G%jsd:G%jed), source=0.0)
    call MOM_read_data(khth_file, khth_varname, CS%khth2d(:,:), G%domain, scale=US%m_to_L**2*US%T_to_s)
    call pass_var(CS%khth2d, G%domain)
  endif
  call get_param(param_file, mdl, "KHTH_SLOPE_CFF", CS%KHTH_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula for "//&
                 "the interface depth diffusivity", units="nondim", default=0.0)
  call get_param(param_file, mdl, "KHTH_MIN", CS%KHTH_Min, &
                 "The minimum horizontal thickness diffusivity.", &
                 default=0.0, units="m2 s-1", scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "KHTH_USE_EBT_STRUCT", khth_use_ebt_struct, &
                 "If true, uses the equivalent barotropic structure "//&
                 "as the vertical structure of thickness diffusivity.",&
                 default=.false., do_not_log=.true.)
  if (khth_use_ebt_struct .and. CS%KHTH_Min>0.0) then
    call get_param(param_file, mdl, "FULL_DEPTH_KHTH_MIN", CS%full_depth_khth_min, &
                   "If true, KHTH_MIN is enforced throughout the whole water column. "//&
                   "Otherwise, KHTH_MIN is only enforced at the surface. This parameter "//&
                   "is only available when KHTH_USE_EBT_STRUCT=True and KHTH_MIN>0.",      &
                   default=.false.)
  endif
  call get_param(param_file, mdl, "KHTH_MAX", CS%KHTH_Max, &
                 "The maximum horizontal thickness diffusivity.", &
                 default=0.0, units="m2 s-1", scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "KHTH_MAX_CFL", CS%max_Khth_CFL, &
                 "The maximum value of the local diffusive CFL ratio that "//&
                 "is permitted for the thickness diffusivity. 1.0 is the "//&
                 "marginally unstable value in a pure layered model, but "//&
                 "much smaller numbers (e.g. 0.1) seem to work better for "//&
                 "ALE-based models.", units="nondimensional", default=0.8)

  call get_param(param_file, mdl, "KH_ETA_CONST", CS%Kh_eta_bg, &
                 "The background horizontal diffusivity of the interface heights (without "//&
                 "considering the layer density structure).  If diffusive CFL limits are "//&
                 "encountered, the diffusivities of the isopycnals and the interfaces heights "//&
                 "are scaled back proportionately.", &
                 default=0.0, units="m2 s-1", scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "KH_ETA_VEL_SCALE", CS%Kh_eta_vel, &
                 "A velocity scale that is multiplied by the grid spacing to give a contribution "//&
                 "to the horizontal diffusivity of the interface heights (without considering "//&
                 "the layer density structure).", &
                 default=0.0, units="m s-1", scale=US%m_to_L*US%T_to_s)

  if ((CS%Kh_eta_bg > 0.0) .or. (CS%Kh_eta_vel > 0.0)) then
    allocate(CS%Kh_eta_u(G%IsdB:G%IedB, G%jsd:G%jed), source=0.)
    allocate(CS%Kh_eta_v(G%isd:G%ied, G%JsdB:G%JedB), source=0.)
    do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
      grid_sp = sqrt((2.0*G%dxCu(I,j)**2 * G%dyCu(I,j)**2) / ((G%dxCu(I,j)**2) + (G%dyCu(I,j)**2)))
      CS%Kh_eta_u(I,j) = G%OBCmaskCu(I,j) * MAX(0.0, CS%Kh_eta_bg + CS%Kh_eta_vel * grid_sp)
    enddo ; enddo
    do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
      grid_sp = sqrt((2.0*G%dxCv(i,J)**2 * G%dyCv(i,J)**2) / ((G%dxCv(i,J)**2) + (G%dyCv(i,J)**2)))
      CS%Kh_eta_v(i,J) = G%OBCmaskCv(i,J) * MAX(0.0, CS%Kh_eta_bg + CS%Kh_eta_vel * grid_sp)
    enddo ; enddo
  endif

  if (CS%max_Khth_CFL < 0.0) CS%max_Khth_CFL = 0.0
  call get_param(param_file, mdl, "DETANGLE_INTERFACES", CS%detangle_interfaces, &
                 "If defined add 3-d structured enhanced interface height "//&
                 "diffusivities to horizontally smooth jagged layers.", &
                 default=.false.)
  CS%detangle_time = 0.0
  if (CS%detangle_interfaces) &
    call get_param(param_file, mdl, "DETANGLE_TIMESCALE", CS%detangle_time, &
                 "A timescale over which maximally jagged grid-scale "//&
                 "thickness variations are suppressed.  This must be "//&
                 "longer than DT, or 0 to use DT.", units="s", default=0.0, scale=US%s_to_T)
  call get_param(param_file, mdl, "KHTH_SLOPE_MAX", CS%slope_max, &
                 "A slope beyond which the calculated isopycnal slope is "//&
                 "not reliable and is scaled away.", units="nondim", default=0.01, scale=US%L_to_Z)
  call get_param(param_file, mdl, "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=GV%m2_s_to_HZ_T)
  call get_param(param_file, mdl, "KHTH_USE_FGNV_STREAMFUNCTION", CS%use_FGNV_streamfn, &
                 "If true, use the streamfunction formulation of "//&
                 "Ferrari et al., 2010, which effectively emphasizes "//&
                 "graver vertical modes by smoothing in the vertical.",  &
                 default=.false.)
  call get_param(param_file, mdl, "FGNV_FILTER_SCALE", CS%FGNV_scale, &
                 "A coefficient scaling the vertical smoothing term in the "//&
                 "Ferrari et al., 2010, streamfunction formulation.", &
                 units="nondim", default=1., do_not_log=.not.CS%use_FGNV_streamfn)
  call get_param(param_file, mdl, "FGNV_C_MIN", CS%FGNV_c_min, &
                 "A minium wave speed used in the Ferrari et al., 2010, "//&
                 "streamfunction formulation.", &
                 default=0., units="m s-1", scale=US%m_s_to_L_T, do_not_log=.not.CS%use_FGNV_streamfn)
  call get_param(param_file, mdl, "FGNV_STRAT_FLOOR", strat_floor, &
                 "A floor for Brunt-Vasaila frequency in the Ferrari et al., 2010, "//&
                 "streamfunction formulation, expressed as a fraction of planetary "//&
                 "rotation, OMEGA. This should be tiny but non-zero to avoid degeneracy.", &
                 default=1.e-15, units="nondim", scale=US%Z_to_L, do_not_log=.not.CS%use_FGNV_streamfn)
  call get_param(param_file, mdl, "STOCH_EOS", stoch_eos, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "USE_STANLEY_GM", CS%use_stanley_gm, &
                 "If true, turn on Stanley SGS T variance parameterization "// &
                 "in GM code.", default=.false.)
  if (CS%use_Stanley_GM .and. .not.stoch_eos) then
    call MOM_error(FATAL, "thickness_diffuse_init: USE_STANLEY_GM requires STOCH_EOS")
  endif
  call get_param(param_file, mdl, "OMEGA", omega, &
                 "The rotation rate of the earth.", &
                 default=7.2921e-5, units="s-1", scale=US%T_to_s, do_not_log=.not.CS%use_FGNV_streamfn)
  CS%N2_floor = 0.
  if (CS%use_FGNV_streamfn) CS%N2_floor = (strat_floor*omega)**2
  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231, do_not_log=.true.)

  call get_param(param_file, mdl, "MEKE_GM_SRC_ALT", CS%GM_src_alt, &
                 "If true, use the GM energy conversion form S^2*N^2*kappa rather "//&
                 "than the streamfunction for the GM source term.", default=.false.)
  call get_param(param_file, mdl, "MEKE_GM_SRC_ANSWER_DATE", CS%MEKE_src_answer_date, &
                 "The vintage of the expressions in the GM energy conversion calculation when "//&
                 "MEKE_GM_SRC_ALT is true.  Values below 20240601 recover the answers from the "//&
                 "original implementation, while higher values use expressions that satisfy "//&
                 "rotational symmetry.", &
                 default=default_answer_date, do_not_log=.not.CS%GM_src_alt)
  call get_param(param_file, mdl, "MEKE_GM_SRC_ALT_SLOPE_BUG", CS%MEKE_src_slope_bug, &
                 "If true, use a bug that limits the positive values, but not the negative values, "//&
                 "of the slopes used when MEKE_GM_SRC_ALT is true.  When this is true, it breaks "//&
                 "all of the symmetry rules that MOM6 is supposed to obey.", &
                 default=.false., do_not_log=.not.CS%GM_src_alt)

  call get_param(param_file, mdl, "MEKE_GEOMETRIC", CS%MEKE_GEOMETRIC, &
                 "If true, uses the GM coefficient formulation from the GEOMETRIC "//&
                 "framework (Marshall et al., 2012).", default=.false.)
  if (CS%MEKE_GEOMETRIC) then
    call get_param(param_file, mdl, "MEKE_GEOMETRIC_EPSILON", CS%MEKE_GEOMETRIC_epsilon, &
                 "Minimum Eady growth rate used in the calculation of GEOMETRIC "//&
                 "thickness diffusivity.", units="s-1", default=1.0e-7, scale=US%T_to_s)
    call get_param(param_file, mdl, "MEKE_GEOMETRIC_ALPHA", CS%MEKE_GEOMETRIC_alpha, &
                 "The nondimensional coefficient governing the efficiency of the GEOMETRIC "//&
                 "thickness diffusion.", units="nondim", default=0.05)

    call get_param(param_file, mdl, "MEKE_GEOMETRIC_ANSWER_DATE", CS%MEKE_GEOM_answer_date, &
                 "The vintage of the expressions in the MEKE_GEOMETRIC calculation.  "//&
                 "Values below 20190101 recover the answers from the original implementation, "//&
                 "while higher values use expressions that satisfy rotational symmetry.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
    if (.not.GV%Boussinesq) CS%MEKE_GEOM_answer_date = max(CS%MEKE_GEOM_answer_date, 20230701)
  endif

  call get_param(param_file, mdl, "USE_MEKE", use_meke, default=.false., do_not_log=.true.)
  if (use_meke) then
    call get_param(param_file, mdl, "USE_KH_IN_MEKE", CS%Use_KH_in_MEKE, &
                   "If true, uses the thickness diffusivity calculated here to diffuse MEKE.", &
                   default=.false.)
    call get_param(param_file, mdl, "MEKE_MIN_DEPTH_DIFF", CS%MEKE_min_depth_diff, &
                   "The minimum total depth over which to average the diffusivity used for MEKE.  "//&
                   "When the total depth is less than this, the diffusivity is scaled away.", &
                   units="m", default=1.0, scale=GV%m_to_H, do_not_log=.not.CS%Use_KH_in_MEKE)
  else
    CS%Use_KH_in_MEKE = .false.
  endif

  call get_param(param_file, mdl, "USE_GME", CS%use_GME_thickness_diffuse, &
                 "If true, use the GM+E backscatter scheme in association "//&
                 "with the Gent and McWilliams parameterization.", default=.false.)

  call get_param(param_file, mdl, "USE_GM_WORK_BUG", CS%use_GM_work_bug, &
                 "If true, compute the top-layer work tendency on the u-grid "//&
                 "with the incorrect sign, for legacy reproducibility.", &
                 default=.false.)

  if (CS%use_GME_thickness_diffuse) then
    allocate(CS%KH_u_GME(G%IsdB:G%IedB, G%jsd:G%jed, GV%ke+1), source=0.)
    allocate(CS%KH_v_GME(G%isd:G%ied, G%JsdB:G%JedB, GV%ke+1), source=0.)
  endif

  CS%id_uhGM = register_diag_field('ocean_model', 'uhGM', diag%axesCuL, Time, &
           'Time Mean Diffusive Zonal Thickness Flux', &
           'kg s-1', conversion=GV%H_to_kg_m2*US%L_to_m**2*US%s_to_T, &
           y_cell_method='sum', v_extensive=.true.)
  if (CS%id_uhGM > 0) call safe_alloc_ptr(CDp%uhGM,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  CS%id_vhGM = register_diag_field('ocean_model', 'vhGM', diag%axesCvL, Time, &
           'Time Mean Diffusive Meridional Thickness Flux', &
           'kg s-1', conversion=GV%H_to_kg_m2*US%L_to_m**2*US%s_to_T, &
           x_cell_method='sum', v_extensive=.true.)
  if (CS%id_vhGM > 0) call safe_alloc_ptr(CDp%vhGM,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)

  CS%id_GMwork = register_diag_field('ocean_model', 'GMwork', diag%axesT1, Time, &
          'Integrated Tendency of Ocean Mesoscale Eddy KE from Parameterized Eddy Advection', &
          'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2, cmor_field_name='tnkebto', &
          cmor_long_name='Integrated Tendency of Ocean Mesoscale Eddy KE from Parameterized Eddy Advection', &
          cmor_standard_name='tendency_of_ocean_eddy_kinetic_energy_content_due_to_parameterized_eddy_advection')
  if (CS%id_GMwork > 0) &
    allocate(CS%GMwork(G%isd:G%ied,G%jsd:G%jed), source=0.)

  CS%id_KH_u = register_diag_field('ocean_model', 'KHTH_u', diag%axesCui, Time, &
           'Parameterized mesoscale eddy advection diffusivity at U-point', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  CS%id_KH_v = register_diag_field('ocean_model', 'KHTH_v', diag%axesCvi, Time, &
           'Parameterized mesoscale eddy advection diffusivity at V-point', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  CS%id_KH_t = register_diag_field('ocean_model', 'KHTH_t', diag%axesTL, Time, &
          'Ocean Tracer Diffusivity due to Parameterized Mesoscale Advection', &
          'm2 s-1', conversion=US%L_to_m**2*US%s_to_T, &
          cmor_field_name='diftrblo', &
          cmor_long_name='Ocean Tracer Diffusivity due to Parameterized Mesoscale Advection', &
          cmor_standard_name='ocean_tracer_diffusivity_due_to_parameterized_mesoscale_advection')

  CS%id_KH_u1 = register_diag_field('ocean_model', 'KHTH_u1', diag%axesCu1, Time,         &
           'Parameterized mesoscale eddy advection diffusivity at U-points (2-D)', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  CS%id_KH_v1 = register_diag_field('ocean_model', 'KHTH_v1', diag%axesCv1, Time,         &
           'Parameterized mesoscale eddy advection diffusivity at V-points (2-D)', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  CS%id_KH_t1 = register_diag_field('ocean_model', 'KHTH_t1', diag%axesT1, Time, &
           'Parameterized mesoscale eddy advection diffusivity at T-points (2-D)', &
           'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)

  CS%id_slope_x =  register_diag_field('ocean_model', 'neutral_slope_x', diag%axesCui, Time, &
           'Zonal slope of neutral surface', 'nondim', conversion=US%Z_to_L)
  if (CS%id_slope_x > 0) &
    allocate(CS%diagSlopeX(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke+1), source=0.)

  CS%id_slope_y =  register_diag_field('ocean_model', 'neutral_slope_y', diag%axesCvi, Time, &
           'Meridional slope of neutral surface', 'nondim', conversion=US%Z_to_L)
  if (CS%id_slope_y > 0) &
    allocate(CS%diagSlopeY(G%isd:G%ied,G%JsdB:G%JedB,GV%ke+1), source=0.)

  CS%id_sfn_x =  register_diag_field('ocean_model', 'GM_sfn_x', diag%axesCui, Time, &
           'Parameterized Zonal Overturning Streamfunction', &
           'm3 s-1', conversion=GV%H_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_sfn_y =  register_diag_field('ocean_model', 'GM_sfn_y', diag%axesCvi, Time, &
           'Parameterized Meridional Overturning Streamfunction', &
           'm3 s-1', conversion=GV%H_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_sfn_unlim_x =  register_diag_field('ocean_model', 'GM_sfn_unlim_x', diag%axesCui, Time, &
           'Parameterized Zonal Overturning Streamfunction before limiting/smoothing', &
           'm3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_sfn_unlim_y =  register_diag_field('ocean_model', 'GM_sfn_unlim_y', diag%axesCvi, Time, &
           'Parameterized Meridional Overturning Streamfunction before limiting/smoothing', &
           'm3 s-1', conversion=US%Z_to_m*US%L_to_m**2*US%s_to_T)

end subroutine thickness_diffuse_init

!> Copies KH_u_GME and KH_v_GME from private type into arrays provided as arguments
subroutine thickness_diffuse_get_KH(CS, KH_u_GME, KH_v_GME, G, GV)
  type(thickness_diffuse_CS),          intent(in)  :: CS   !< Control structure for this module
  type(ocean_grid_type),               intent(in)  :: G    !< Grid structure
  type(verticalGrid_type),             intent(in)  :: GV   !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: KH_u_GME !< Isopycnal height
                                                   !! diffusivities at u-faces [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), intent(inout) :: KH_v_GME !< Isopycnal height
                                                   !! diffusivities at v-faces [L2 T-1 ~> m2 s-1]
  ! Local variables
  integer :: i,j,k

  do k=1,GV%ke+1 ; do j = G%jsc, G%jec ; do I = G%isc-1, G%iec
    KH_u_GME(I,j,k) = CS%KH_u_GME(I,j,k)
  enddo ; enddo ; enddo

  do k=1,GV%ke+1 ; do J = G%jsc-1, G%jec ; do i = G%isc, G%iec
    KH_v_GME(i,J,k) = CS%KH_v_GME(i,J,k)
  enddo ; enddo ; enddo

end subroutine thickness_diffuse_get_KH

!> Deallocate the thickness_diffus3 control structure
subroutine thickness_diffuse_end(CS, CDp)
  type(thickness_diffuse_CS), intent(inout) :: CS !< Control structure for thickness_diffuse
  type(cont_diag_ptrs), intent(inout) :: CDp      !< Continuity diagnostic control structure

  if (CS%id_slope_x > 0) deallocate(CS%diagSlopeX)
  if (CS%id_slope_y > 0) deallocate(CS%diagSlopeY)

  if (CS%id_GMwork > 0) deallocate(CS%GMwork)

  ! NOTE: [uv]hGM may be allocated either here or the diagnostic module
  if (associated(CDp%uhGM)) deallocate(CDp%uhGM)
  if (associated(CDp%vhGM)) deallocate(CDp%vhGM)

  if (CS%use_GME_thickness_diffuse) then
    deallocate(CS%KH_u_GME)
    deallocate(CS%KH_v_GME)
  endif

  if (allocated(CS%khth2d)) deallocate(CS%khth2d)
end subroutine thickness_diffuse_end

!> \namespace mom_thickness_diffuse
!!
!! \section section_gm Isopycnal height diffusion (aka Gent-McWilliams)
!!
!! Isopycnal height diffusion is implemented via along-layer mass fluxes
!! \f[
!! h^\dagger \leftarrow h^n - \Delta t \nabla \cdot ( \vec{uh}^* )
!! \f]
!! where the mass fluxes are cast as the difference in vector streamfunction
!!
!! \f[
!! \vec{uh}^* = \delta_k \vec{\psi} .
!! \f]
!!
!! The GM implementation of isopycnal height diffusion made the streamfunction proportional
!! to the potential density slope
!! \f[
!! \vec{\psi} = - \kappa_h \frac{\nabla_z \rho}{\partial_z \rho}
!! = \frac{g\kappa_h}{\rho_o} \frac{\nabla \rho}{N^2} = -\kappa_h \frac{M^2}{N^2}
!! \f]
!! but for robustness the scheme is implemented as
!! \f[
!! \vec{\psi} = -\kappa_h \frac{M^2}{\sqrt{N^4 + M^4}}
!! \f]
!! since the quantity \f$\frac{M^2}{\sqrt{N^4 + M^4}}\f$ is bounded between $-1$ and $1$ and does not change sign
!! if \f$N^2<0\f$.
!!
!! Optionally, the method of \cite ferrari2010, can be used to obtain the streamfunction which solves the
!! vertically elliptic equation:
!! \f[
!! \gamma_F \partial_z c^2 \partial_z \psi - N_*^2 \psi  = -( 1 + \gamma_F ) \kappa_h N_*^2 \frac{M^2}{\sqrt{N^4+M^4}}
!! \f]
!! which recovers the previous streamfunction relation in the limit that \f$ c \rightarrow 0 \f$.
!! Here, \f$c=\max(c_{min},c_g)\f$ is the maximum of either \f$c_{min}\f$ and either the first baroclinic mode
!! wave-speed or the equivalent barotropic mode wave-speed.
!! \f$N_*^2 = \max(N^2,0)\f$ is a non-negative form of the square of the buoyancy frequency.
!! The parameter \f$\gamma_F\f$ is used to reduce the vertical smoothing length scale.
!! \f[
!! \kappa_h = \left( \kappa_o + \alpha_{s} L_{s}^2 < S N > + \alpha_{M} \kappa_{M} \right) r(\Delta x,L_d)
!! \f]
!! where \f$ S \f$ is the isoneutral slope magnitude, \f$ N \f$ is the buoyancy frequency,
!! \f$\kappa_{M}\f$ is the diffusivity calculated by the MEKE parameterization (mom_meke module) and
!! \f$ r(\Delta x,L_d) \f$ is a function of the local resolution (ratio of grid-spacing, \f$\Delta x\f$,
!! to deformation radius, \f$L_d\f$). The length \f$L_s\f$ is provided by the mom_lateral_mixing_coeffs module
!! (enabled with <code>USE_VARIABLE_MIXING=True</code> and the term \f$<SN>\f$ is the vertical average slope
!! times the buoyancy frequency prescribed by \cite visbeck1996.
!!
!! The result of the above expression is subsequently bounded by minimum and maximum values, including an upper
!! diffusivity consistent with numerical stability (\f$ \kappa_{cfl} \f$ is calculated internally).
!! \f[
!! \kappa_h \leftarrow \min{\left( \kappa_{max}, \kappa_{cfl}, \max{\left( \kappa_{min}, \kappa_h \right)} \right)}
!!                      f(c_g,z)
!! \f]
!!
!! where \f$f(c_g,z)\f$ is a vertical structure function.
!! \f$f(c_g,z)\f$ is calculated in module mom_lateral_mixing_coeffs.
!! If <code>KHTH_USE_EBT_STRUCT=True</code> then \f$f(c_g,z)\f$ is set to look like the equivalent barotropic
!! modal velocity structure. Otherwise \f$f(c_g,z)=1\f$ and the diffusivity is independent of depth.
!!
!! In order to calculate meaningful slopes in vanished layers, temporary copies of the thermodynamic variables
!! are passed through a vertical smoother, function vert_fill_ts():
!! \f{eqnarray*}{
!! \left[ 1 + \Delta t \kappa_{smth} \frac{\partial^2}{\partial_z^2} \right] \theta & \leftarrow & \theta \\
!! \left[ 1 + \Delta t \kappa_{smth} \frac{\partial^2}{\partial_z^2} \right] s & \leftarrow & s
!! \f}
!!
!! \subsection section_khth_module_parameters Module mom_thickness_diffuse parameters
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>THICKNESSDIFFUSE</code> |
!! | \f$ \kappa_o \f$      | <code>KHTH</code> |
!! | \f$ \alpha_{s} \f$    | <code>KHTH_SLOPE_CFF</code> |
!! | \f$ \kappa_{min} \f$  | <code>KHTH_MIN</code> |
!! | \f$ \kappa_{max} \f$  | <code>KHTH_MAX</code> |
!! | -                     | <code>KHTH_MAX_CFL</code> |
!! | \f$ \kappa_{smth} \f$ | <code>KD_SMOOTH</code> |
!! | \f$ \alpha_{M} \f$    | <code>MEKE_KHTH_FAC</code> (from mom_meke module) |
!! | -                     | <code>KHTH_USE_EBT_STRUCT</code> (from mom_lateral_mixing_coeffs module) |
!! | -                     | <code>KHTH_USE_FGNV_STREAMFUNCTION</code> |
!! | \f$ \gamma_F \f$      | <code>FGNV_FILTER_SCALE</code> |
!! | \f$ c_{min} \f$       | <code>FGNV_C_MIN</code> |
!!
!! \subsection section_khth_module_reference References
!!
!! Ferrari, R., S.M. Griffies, A.J.G. Nurser and G.K. Vallis, 2010:
!! A boundary-value problem for the parameterized mesoscale eddy transport.
!! Ocean Modelling, 32, 143-156. http://doi.org/10.1016/j.ocemod.2010.01.004
!!
!! Visbeck, M., J.C. Marshall, H. Jones, 1996:
!! Dynamics of isolated convective regions in the ocean. J. Phys. Oceangr., 26, 1721-1734.
!! http://dx.doi.org/10.1175/1520-0485(1996)026%3C1721:DOICRI%3E2.0.CO;2

end module MOM_thickness_diffuse
