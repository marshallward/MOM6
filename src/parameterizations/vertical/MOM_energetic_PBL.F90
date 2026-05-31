! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

#include <MOM_memory.h>

!> Energetically consistent planetary boundary layer parameterization
module MOM_energetic_PBL

use MOM_coms, only : EFP_type, operator(+), assignment(=)
use MOM_coms, only : EFP_to_real
use MOM_coms, only : real_to_EFP
use MOM_coms, only : EFP_sum_across_PEs
use MOM_diag_mediator, only : diag_ctrl
use MOM_diag_mediator, only : register_diag_field
use MOM_diag_mediator, only : safe_alloc_alloc
use MOM_diag_mediator, only : time_type
use MOM_error_handler, only : MOM_error
use MOM_error_handler, only : MOM_mesg
use MOM_error_handler, only : FATAL, WARNING
use MOM_file_parser, only : get_param
use MOM_file_parser, only : log_param
use MOM_file_parser, only : log_version
use MOM_file_parser, only : param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_stochastics, only : stochastic_CS
use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_variables, only : vertvisc_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS

implicit none
private

public :: energetic_PBL
public :: energetic_PBL_init
public :: energetic_PBL_end
public :: energetic_PBL_get_MLD

!> This control structure holds parameters for the MOM_energetic_PBL module
type, public :: energetic_PBL_CS
  private
  logical :: initialized = .false.
    !< True if this control structure has been initialized.

  !/ Constants
  real :: VonKar
    !< The von Karman coefficient as used in the ePBL module [nondim]
  real :: omega
    !< The Earth's rotation rate [T-1 ~> s-1].
  real :: omega_frac
    !< When setting the decay scale for turbulence, use this fraction of the
    !! absolute rotation rate blended with the local value of f, as
    !! sqrt((1-omega_frac)*f^2 + omega_frac*4*omega^2) [nondim].

  !/ Convection related terms
  real :: nstar
    !< The fraction of the TKE input to the mixed layer available to drive
    !! entrainment [nondim]. This quantity is the vertically integrated
    !! buoyancy production minus the vertically integrated dissipation of
    !! TKE produced by buoyancy.

  !/ Mixing Length terms
  logical :: Use_MLD_iteration
    !< If true, use the proximity to the bottom of the actively turbulent
    !! surface boundary layer to constrain the mixing lengths.
  logical :: MLD_iteration_guess
    !< False to default to guessing half the
    !! ocean depth for the first iteration.
  logical :: MLD_bisection
    !< If true, use bisection with the iterative determination of the
    !! self-consistent mixed layer depth.  Otherwise use the false position
    !! after a maximum and minimum bound have been evaluated and the
    !! returned value from the previous guess or bisection before this.
  logical :: MLD_iter_bug
    !< If true use buggy logic that gives the wrong bounds for the next
    !! iteration when successive guesses increase by exactly EPBL_MLD_TOLERANCE.
  integer :: max_MLD_its
    !< The maximum number of iterations that can be used to find a
    !! self-consistent mixed layer depth with Use_MLD_iteration.
  real :: MixLenExponent
    !< Exponent in the mixing length shape-function [nondim].
    !! 1 is law-of-the-wall at top and bottom,
    !! 2 is more KPP like.
  real :: MKE_to_TKE_effic
    !< The efficiency with which mean kinetic energy released by
    !!  mechanically forced entrainment of the mixed layer is converted to
    !!  TKE, times conversion factors between the natural units of mean
    !!  kinetic energy and those used for TKE [Z2 L-2 ~> nondim].
  logical :: direct_calc
    !< If true and there is no conversion from mean kinetic energy to ePBL
    !! turbulent kinetic energy, use a direct calculation of the
    !! diffusivity that is supported by a given energy input instead of the
    !! more general but slower iterative solver.
  real :: ustar_min
    !< A minimum value of ustar to avoid numerical problems [Z T-1 ~> m s-1].
    !! If the value is small enough, this should not affect the solution.
  real :: Ekman_scale_coef
    !< A nondimensional scaling factor controlling the inhibition of the
    !! diffusive length scale by rotation [nondim].  Making this larger decreases
    !! the diffusivity in the planetary boundary layer.
  real :: transLay_scale
    !< A scale for the mixing length in the transition layer
    !! at the edge of the boundary layer as a fraction of the
    !! boundary layer thickness [nondim].  The default is 0, but a
    !! value of 0.1 might be better justified by observations.
  real :: MLD_tol
    !< A tolerance for determining the boundary layer thickness when
    !! Use_MLD_iteration is true [Z ~> m].
  real :: min_mix_len
    !< The minimum mixing length scale that will be used by ePBL [Z ~> m].
    !! The default (0) does not set a minimum.

  !/ Velocity scale terms
  integer :: wT_scheme
    !< An enumerated value indicating the method for finding the turbulent
    !! velocity scale.  There are currently two options:
    !! wT_mwT_from_cRoot_TKE is the original (TKE_remaining)^1/3
    !! wT_from_RH18 is the version described by Reichl and Hallberg, 2018
  real :: wstar_ustar_coef
    !< A ratio relating the efficiency with which convectively released
    !! energy is converted to a turbulent velocity, relative to
    !! mechanically forced turbulent kinetic energy [nondim].
    !! Making this larger increases the diffusivity.
  real :: vstar_surf_fac
    !< If (wT_scheme == wT_from_RH18) this is the proportionality coefficient between
    !! ustar and the surface mechanical contribution to vstar [nondim]
  real :: vstar_scale_fac
    !< An overall nondimensional scaling factor for vstar [nondim].  Making
    !! this larger increases the diffusivity.

  !mstar related options
  integer :: mstar_scheme
    !< An encoded integer to determine which formula is used to set mstar
  integer :: BBL_mstar_scheme
    !< An encoded integer to determine which formula is used to set mstar
  real :: mstar_cap
    !< Since mstar is restoring undissipated energy to mixing,
    !! there must be a cap on how large it can be [nondim].  This
    !! is definitely a function of latitude (Ekman limit),
    !! but will be taken as constant for now.

  !/ vertical decay related options
  real :: TKE_decay
    !< The ratio of the natural Ekman depth to the TKE decay scale [nondim].

  !/ mstar_scheme == 0
  real :: fixed_mstar
    !< mstar is the ratio of the friction velocity cubed to the TKE available to
    !! drive entrainment [nondim]. This quantity is the vertically
    !! integrated shear production minus the vertically integrated
    !! dissipation of TKE produced by shear.  This value is used if the option
    !! for using a fixed mstar is used.
  real :: BBL_fixed_mstar
    !< Similar to fixed_mstar, but for the bottom boundary layer

  !/ mstar_scheme == 2
  real :: C_Ek = 0.17
    !< mstar Coefficient in rotation limit for EPBL_MSTAR_SCHEME=OM4 [nondim]
  real :: mstar_coef = 0.3
    !< mstar coefficient in rotation/stabilizing balance for
    !! EPBL_MSTAR_SCHEME=OM4 [nondim]

  !/ mstar_scheme == 3
  real :: RH18_mstar_cN1
    !< mstar_N coefficient 1 (outer-most coefficient for fit) [nondim].
    !! Value of 0.275 in RH18.  Increasing this
    !! coefficient increases mechanical mixing for all values of Hf/ust,
    !! but is most effective at low values (weakly developed OSBLs).
  real :: RH18_mstar_cN2
    !< mstar_N coefficient 2 (coefficient outside of exponential decay) [nondim].
    !! Value of 8.0 in RH18.  Increasing this coefficient increases mstar
    !! for all values of HF/ust, with a consistent affect across
    !! a wide range of Hf/ust.
  real :: RH18_mstar_cN3
    !< mstar_N coefficient 3 (exponential decay coefficient) [nondim]. Value of
    !! -5.0 in RH18.  Increasing this increases how quickly the value
    !! of mstar decreases as Hf/ust increases.
  real :: RH18_mstar_cS1
    !< mstar_S coefficient for RH18 in stabilizing limit [nondim].
    !! Value of 0.2 in RH18.
  real :: RH18_mstar_cS2
    !< mstar_S exponent for RH18 in stabilizing limit [nondim].
    !! Value of 0.4 in RH18.

  !/ Coefficient for shear/convective turbulence interaction
  real :: mstar_convect_coef
    !< Factor to reduce mstar when statically unstable [nondim].

  !/ Langmuir turbulence related parameters
  logical :: Use_LT = .false.
    !< Flag for using LT in Energy calculation
  integer :: LT_enhance_form
    !< Integer for Enhancement functional form (various options)
  real :: LT_enhance_coef
    !< Coefficient in fit for Langmuir Enhancement [nondim]
  real :: LT_enhance_exp
    !< Exponent in fit for Langmuir Enhancement [nondim]
  real :: LaC_MLD_Ek
    !< Coefficient for Langmuir number modification based on the ratio of
    !! the mixed layer depth over the Ekman depth [nondim].
  real :: LaC_MLD_Ob_stab
    !< Coefficient for Langmuir number modification based on the ratio of
    !! the mixed layer depth over the Obukhov depth with stabilizing forcing [nondim].
  real :: LaC_Ek_Ob_stab
    !< Coefficient for Langmuir number modification based on the ratio of
    !! the Ekman depth over the Obukhov depth with stabilizing forcing [nondim].
  real :: LaC_MLD_Ob_un
    !< Coefficient for Langmuir number modification based on the ratio of
    !! the mixed layer depth over the Obukhov depth with destabilizing forcing [nondim].
  real :: LaC_Ek_Ob_un
    !< Coefficient for Langmuir number modification based on the ratio of
    !! the Ekman depth over the Obukhov depth with destabilizing forcing [nondim].
  real :: Max_Enhance_M = 5.
    !< The maximum allowed LT enhancement to the mixing [nondim].

  !/ Machine learned equation discovery model paramters
  logical :: eqdisc
    !< Uses machine learned shape function
  logical :: eqdisc_v0
    !< Uses machine learned velocity scale
  logical :: eqdisc_v0h
    !< Uses machine learned velocity scale that uses boundary layer depth as input
  real :: v0_lower_cap
    !< Lower cap to prevent v0 from attaining anomlously low values [Z T-1 ~> m s-1]
  real :: v0_upper_cap
    !< Upper cap to prevent v0 from attaining anomlously high values [Z T-1 ~> m s-1]
  real :: f_lower
    !< Lower cap of |f| i.e. absolute of Coriolis parameter [T-1 ~> s-1]
    !! Used only in get_eqdisc_v0 subroutine. Default is 0.1deg Lat
  real :: bflux_lower_cap
    !< Lower cap for capping blfux [Z2 T-3 ~> m2 s-3]
  real :: bflux_upper_cap
    !< Upper cap for capping blfux [Z2 T-3 ~> m2 s-3]
  real :: sigma_max_lower_cap
    !< Lower cap to prevent sigma_max from attaining low values [nondim]
  real :: sigma_max_upper_cap
    !< Upper cap to prevent sigma_max from attaining high values [nondim]
  real :: Eh_upper_cap
    !< Upper cap to prevent Eh = hf/(u__*) from attaining high values [nondim]
  real :: Lh_cap
    !< Cap to prevent Lh = h/Monin_Obukhov_depth from attaining beyond extreme values [nondim]
  real, allocatable, dimension(:) :: shape_function
    !< shape function used in machine learned diffusivity [nondim]

  !/ Coefficients used for Machine learned diffusivity
  real :: ML_c(18)
    !< Array of non-dimensional constants used in machine learned (ML) diffusivity [nondim]
  real :: shape_function_epsilon
    !< An small value of shape_function below the boundary layer depth [nondim]

  !/ Bottom boundary layer mixing related options
  real :: ePBL_BBL_effic
    !< The efficiency of bottom boundary layer mixing via ePBL driven by
    !! the bottom drag dissipation of mean kinetic energy, times
    !! conversion factors between the natural units of mean kinetic energy
    !! and those used for TKE [Z2 L-2 ~> nondim].
  real :: ePBL_tidal_effic
    !< The efficiency of bottom boundary layer mixing via ePBL driven by
    !! the bottom drag dissipation of tides, times conversion factors
    !! between the natural units of mean kinetic energy and those used for
    !! TKE [Z2 L-2 ~> nondim].
  logical :: Use_BBLD_iteration
    !< If true, use the proximity to the top of the actively turbulent
    !! bottom boundary layer to constrain the mixing lengths.
  real :: TKE_decay_BBL
    !< The ratio of the natural Ekman depth to the TKE decay scale for
    !! bottom boundary layer mixing [nondim]
  real :: min_BBL_mix_len
    !< The minimum mixing length scale that will be used by ePBL in the bottom
    !! boundary layer mixing [Z ~> m].  The default (0) does not set a minimum.
  real :: MixLenExponent_BBL
    !< Exponent in the bottom boundary layer mixing length shape-function [nondim].
    !! 1 is law-of-the-wall at top and bottom,
    !! 2 is more KPP like.
  real :: BBLD_tol
    !< The tolerance for the iteratively determined bottom boundary layer depth [Z ~> m].
    !! This is only used with USE_MLD_ITERATION.
  integer :: max_BBLD_its
    !< The maximum number of iterations that can be used to find a self-consistent
    !! bottom boundary layer depth.
  integer :: wT_scheme_BBL
    !< An enumerated value indicating the method for finding the bottom boundary
    !! layer turbulent velocity scale.  There are currently two options:
    !! wT_mwT_from_cRoot_TKE is the original (TKE_remaining)^1/3
    !! wT_from_RH18 is the version described by Reichl and Hallberg, 2018
  real :: vstar_scale_fac_BBL
    !< An overall nondimensional scaling factor for wT in the bottom boundary layer [nondim].
    !! Making this larger increases the bottom boundary layer diffusivity.", &
  real :: vstar_surf_fac_BBL
    !< If (wT_scheme_BBL == wT_from_RH18) this is the proportionality coefficient between
    !! ustar and the bottom boundayer layer mechanical contribution to vstar [nondim]
  real :: Ekman_scale_coef_BBL
    !< A nondimensional scaling factor controlling the inhibition of the
    !! diffusive length scale by rotation in the bottom boundary layer [nondim].
    !! Making this larger decreases the bottom boundary layer diffusivity.
  logical :: decay_adjusted_BBL_TKE
    !< If true, include an adjustment factor in the bottom boundary layer
    !! energetics that accounts for an exponential decay of TKE from a
    !! near-bottom source and an assumed piecewise linear linear profile
    !! of the buoyancy flux response to a change in a diffusivity.
  logical :: BBL_effic_bug
    !< If true, overestimate the efficiency of the non-tidal ePBL bottom boundary
    !! layer diffusivity by a factor of 1/sqrt(CDRAG), which is often a factor of
    !! about 18.3.
  logical :: ePBL_BBL_use_mstar
    !< If true, use an mstar*ustar^3 paramaterization to get the TKE available
    !! to drive mixing in the bottom boundary layer version of ePBL.  Otherwise,
    !! use the meanflow energy loss to bottom drag scaled by a constant efficiency.

  !/ Options for documenting differences from parameter choices
  integer :: options_diff
    !< If positive, this is a coded integer indicating a pair of
    !! settings whose differences are diagnosed in a passive diagnostic mode
    !! via extra calls to ePBL_column.  If this is 0 or negative no extra
    !! calls occur.

  !/ Others
  type(time_type), pointer :: Time=>NULL()
    !< A pointer to the ocean model's clock.

  logical :: TKE_diagnostics = .true.
    !< If true, diagnostics of the TKE budget are being calculated.
  integer :: answer_date
    !< The vintage of the order of arithmetic and expressions in the ePBL
    !! calculations.  Values below 20190101 recover the answers from the
    !! end of 2018, while higher values use updated and more robust forms
    !! of the same expressions.  Values below 20240101 use A**(1./3.) to
    !! estimate the cube root of A in several expressions, while higher
    !! values use the integer root function cuberoot(A) and therefore
    !! can work with scaled variables.
  logical :: orig_PE_calc
    !< If true, the ePBL code uses the original form of the
    !! potential energy change code.  Otherwise, it uses a newer version
    !! that can work with successive increments to the diffusivity in
    !! upward or downward passes.
  logical :: debug
    !< If true, write verbose checksums for debugging purposes.
  type(diag_ctrl), pointer :: diag=>NULL()
    !< A structure that is used to regulate the
    !! timing of diagnostic output.

  real, allocatable, dimension(:,:) :: ML_depth
    !< The mixed layer depth determined by active mixing in ePBL, which may
    !! be used for the first guess in the next time step [H ~> m or kg m-2]
  real, allocatable, dimension(:,:) :: BBL_depth
    !< The bottom boundary layer depth determined by active mixing in ePBL [H ~> m or kg m-2]

  type(EFP_type), dimension(2) :: sum_its
    !< The total number of iterations and columns worked on
  type(EFP_type), dimension(2) :: sum_its_BBL
    !< The total number of iterations and columns worked on

  !>@{ Diagnostic IDs
  integer :: id_Kd_ePBL_col_by_col = -1
  integer :: id_ML_depth = -1, id_hML_depth = -1, id_TKE_wind = -1, id_TKE_mixing = -1
  integer :: id_ustar_ePBL = -1, id_bflx_ePBL = -1
  integer :: id_TKE_MKE = -1, id_TKE_conv = -1, id_TKE_forcing = -1
  integer :: id_TKE_mech_decay = -1, id_TKE_conv_decay = -1
  integer :: id_Mixing_Length = -1, id_Velocity_Scale = -1
  integer :: id_Kd_BBL = -1, id_BBL_Mix_Length = -1, id_BBL_Vel_Scale = -1
  integer :: id_TKE_BBL = -1, id_TKE_BBL_mixing = -1, id_TKE_BBL_decay = -1
  integer :: id_ustar_BBL = -1, id_bflx_BBL = -1, id_BBL_decay_scale = -1, id_BBL_depth = -1
  integer :: id_mstar_sfc = -1, id_mstar_BBL = -1, id_LA_mod = -1, id_LA = -1, id_mstar_LT = -1

  ! The next options are used when passively diagnosing sensitivities from parameter choices
  integer :: id_opt_diff_Kd_ePBL = -1, id_opt_maxdiff_Kd_ePBL = -1, id_opt_diff_hML_depth = -1
  !>@}
end type energetic_PBL_CS

! Enumeration values for mstar_scheme
integer, parameter :: Use_Fixed_mstar = 0
  !< The value of mstar_scheme to use a constant mstar
integer, parameter :: mstar_from_Ekman = 2
  !< The value of mstar_scheme to base mstar on the ratio
  !! of the Ekman layer depth to the Obukhov depth
integer, parameter :: mstar_from_RH18 = 3
  !< The value of mstar_scheme to base mstar of of RH18
integer, parameter :: No_Langmuir = 0
  !< The value of LT_enhance_form not use Langmuir turbulence.
integer, parameter :: Langmuir_rescale = 2
  !< The value of LT_enhance_form to use a multiplicative
  !! rescaling of mstar to account for Langmuir turbulence.
integer, parameter :: Langmuir_add = 3
  !< The value of LT_enhance_form to add a contribution to
  !! mstar from Langmuir turbulence to other contributions.
integer, parameter :: wT_from_cRoot_TKE = 0
  !< Use a constant times the cube root of remaining TKE
  !! to calculate the turbulent velocity.
integer, parameter :: wT_from_RH18 = 1
  !< Use a scheme based on a combination of w* and v* as
  !! documented in Reichl & Hallberg (2018) to calculate
  !! the turbulent velocity.

character*(20), parameter :: CONSTANT_STRING = "CONSTANT"
character*(20), parameter :: OM4_STRING = "OM4"
character*(20), parameter :: RH18_STRING = "REICHL_H18"
character*(20), parameter :: ROOT_TKE_STRING = "CUBE_ROOT_TKE"
character*(20), parameter :: NONE_STRING = "NONE"
character*(20), parameter :: RESCALED_STRING = "RESCALE"
character*(20), parameter :: ADDITIVE_STRING = "ADDITIVE"


logical :: report_avg_its = .false.
  !< Report the average number of ePBL iterations for debugging.


interface
  module subroutine energetic_PBL(h_3d, u_3d, v_3d, tv, fluxes, visc, dt, &
      Kd_int, G, GV, US, CS, stoch_CS, dSV_dT, dSV_dS, TKE_forced, &
      buoy_flux, BBL_buoy_flux, Waves)

    type(ocean_grid_type), intent(inout) :: G
      !< The ocean's grid structure.
    type(verticalGrid_type), intent(in) :: GV
      !< The ocean's vertical grid structure.
    type(unit_scale_type), intent(in) :: US
      !< A dimensional unit scaling type
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h_3d
      !< Layer thicknesses [H ~> m or kg m-2]
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: u_3d
      !< Zonal velocities interpolated to h points [L T-1 ~> m s-1]
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: v_3d
      !< Zonal velocities interpolated to h points [L T-1 ~> m s-1]
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: dSV_dT
      !< The partial derivative of in-situ specific
      !! volume with potential temperature
      !! [R-1 C-1 ~> m3 kg-1 degC-1].
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: dSV_dS
      !< The partial derivative of in-situ specific
      !! volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1].
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: TKE_forced
      !< The forcing requirements to homogenize the
      !! forcing that has been applied to each layer
      !! [R Z3 T-2 ~> J m-2].
    type(thermo_var_ptrs), intent(inout) :: tv
      !< A structure containing pointers to any
      !! available thermodynamic fields. Absent fields
      !! have NULL ptrs.
    type(forcing), intent(inout) :: fluxes
      !< A structure containing pointers to any
      !! possible forcing fields. Unused fields have
      !! NULL ptrs.
    type(vertvisc_type), intent(in) :: visc
      !< Structure with vertical viscosities, BBL properties and related fields
    real, intent(in) :: dt
      !< Time increment [T ~> s].
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: Kd_int
      !< The diagnosed diffusivities at interfaces
      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
    type(energetic_PBL_CS),  intent(inout) :: CS
      !< Energetic PBL control structure
    real, dimension(SZI_(G),SZJ_(G)), intent(in) :: buoy_flux
      !< The surface buoyancy flux [Z2 T-3 ~> m2 s-3].
    real, dimension(SZI_(G),SZJ_(G)), intent(in) :: BBL_buoy_flux
      !< The bottom buoyancy flux [Z2 T-3 ~> m2 s-3].
    type(wave_parameters_CS), pointer :: Waves
      !< Waves control structure for Langmuir turbulence
    type(stochastic_CS), pointer :: stoch_CS
      !< The control structure returned by a previous
  end subroutine
end interface

contains

!> Copies the ePBL active mixed layer depth into MLD, in units of [Z ~> m]
!! unless other units are specified.
subroutine energetic_PBL_get_MLD(CS, MLD, G, US, m_to_MLD_units)
  type(energetic_PBL_CS), intent(in) :: CS
    !< Energetic PBL control structure
  type(ocean_grid_type), intent(in) :: G
    !< Grid structure
  type(unit_scale_type), intent(in) :: US
    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: MLD
    !< Depth of ePBL active mixing layer [Z ~> m]
    !! or other units
  real, optional, intent(in) :: m_to_MLD_units
    !< A conversion factor from meters
    !! to the desired units for MLD, sometimes [Z m-1 ~> 1]

  real :: scale
    ! A dimensional rescaling factor, often [nondim] or [m Z-1 ~> 1]
  integer :: i, j

  scale = 1.0
  if (present(m_to_MLD_units)) scale = US%Z_to_m * m_to_MLD_units

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    MLD(i,j) = scale*CS%ML_depth(i,j)
  enddo ; enddo
end subroutine energetic_PBL_get_MLD


!> This subroutine initializes the energetic_PBL module
subroutine energetic_PBL_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(energetic_PBL_CS),  intent(inout) :: CS   !< Energetic PBL control structure

  ! Local variables
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_energetic_PBL"  ! This module's name.
  character(len=20)  :: tmpstr  ! A string that is parsed for parameter settings
  character(len=20)  :: mstar_scheme ! A string that is parsed for mstar parameter settings
  character(len=20)  :: vel_scale_str ! A string that is parsed for velocity scale parameter settings
  character(len=120) :: diff_text ! A clause describing parameter setting that differ.
  real :: omega_frac_dflt  ! The default for omega_frac [nondim]
  integer :: isd, ied, jsd, jed
  integer :: mstar_mode, LT_enhance, wT_mode
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: use_omega
  logical :: no_BBL  ! If true, EPBL_BBL_EFFIC < 0 and EPBL_BBL_TIDAL_EFFIC < 0, so
                     ! bottom boundary layer mixing is not enabled.
  logical :: use_la_windsea
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%initialized = .true.
  CS%diag => diag
  CS%Time => Time

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")


!/1. General ePBL settings
  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "OMEGA", CS%omega, &
                 "The rotation rate of the earth.", &
                 units="s-1", default=7.2921e-5, scale=US%T_to_S)
  call get_param(param_file, mdl, "ML_USE_OMEGA", use_omega, &
                 "If true, use the absolute rotation rate instead of the "//&
                 "vertical component of rotation when setting the decay "//&
                 "scale for turbulence.", default=.false., do_not_log=.true.)
  omega_frac_dflt = 0.0
  if (use_omega) then
    call MOM_error(WARNING, "ML_USE_OMEGA is deprecated; use ML_OMEGA_FRAC=1.0 instead.")
    omega_frac_dflt = 1.0
  endif
  call get_param(param_file, mdl, "ML_OMEGA_FRAC", CS%omega_frac, &
                 "When setting the decay scale for turbulence, use this "//&
                 "fraction of the absolute rotation rate blended with the "//&
                 "local value of f, as sqrt((1-of)*f^2 + of*4*omega^2).", &
                 units="nondim", default=omega_frac_dflt)
  call get_param(param_file, mdl, "EKMAN_SCALE_COEF", CS%Ekman_scale_coef, &
                 "A nondimensional scaling factor controlling the inhibition "//&
                 "of the diffusive length scale by rotation. Making this larger "//&
                 "decreases the PBL diffusivity.", units="nondim", default=1.0)
  call get_param(param_file, mdl, 'VON_KARMAN_CONST', CS%vonKar, &
                 'The value the von Karman constant as used for mixed layer viscosity.', &
                 units='nondim', default=0.41)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "EPBL_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the energetic "//&
                 "PBL calculations.  Values below 20190101 recover the answers from the "//&
                 "end of 2018, while higher values use updated and more robust forms of the "//&
                 "same expressions.  Values below 20240101 use A**(1./3.) to estimate the cube "//&
                 "root of A in several expressions, while higher values use the integer root "//&
                 "function cuberoot(A) and therefore can work with scaled variables.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  call get_param(param_file, mdl, "EPBL_ORIGINAL_PE_CALC", CS%orig_PE_calc, &
                 "If true, the ePBL code uses the original form of the potential energy change "//&
                 "code.  Otherwise, the newer version that can work with successive increments "//&
                 "to the diffusivity in upward or downward passes is used.", &
                 default=.true.) ! Change the default to .false.?

  call get_param(param_file, mdl, "MKE_TO_TKE_EFFIC", CS%MKE_to_TKE_effic, &
                 "The efficiency with which mean kinetic energy released "//&
                 "by mechanically forced entrainment of the mixed layer "//&
                 "is converted to turbulent kinetic energy.", &
                 units="nondim", default=0.0, scale=US%L_to_Z**2)
  call get_param(param_file, mdl, "TKE_DECAY", CS%TKE_decay, &
                 "TKE_DECAY relates the vertical rate of decay of the TKE available "//&
                 "for mechanical entrainment to the natural Ekman depth.", &
                 units="nondim", default=2.5)
  call get_param(param_file, mdl, "DIRECT_EPBL_MIXING_CALC", CS%direct_calc, &
                 "If true and there is no conversion from mean kinetic energy to ePBL turbulent "//&
                 "kinetic energy, use a direct calculation of the diffusivity that is supported "//&
                 "by a given energy input instead of the more general but slower iterative solver.", &
                 default=.false., do_not_log=(CS%MKE_to_TKE_effic>0.0))


!/2. Options related to setting mstar

  call get_param(param_file, mdl, "EPBL_MSTAR_SCHEME", mstar_scheme, &
                 "EPBL_MSTAR_SCHEME selects the method for setting mstar.  Valid values are: \n"//&
                 "\t CONSTANT   - Use a fixed mstar given by MSTAR \n"//&
                 "\t OM4        - Use L_Ekman/L_Obukhov in the stabilizing limit, as in OM4 \n"//&
                 "\t REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.", &
                 default=CONSTANT_STRING, do_not_log=.true.)
  call get_param(param_file, mdl, "MSTAR_MODE", mstar_mode, default=-1)
  if (mstar_mode == 0) then
    mstar_scheme = CONSTANT_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = CONSTANT instead of the archaic MSTAR_MODE = 0.")
  elseif (mstar_mode == 1) then
    call MOM_error(FATAL, "You are using a legacy mstar mode in ePBL that has been phased out. "//&
                          "If you need to use this setting please report this error.  Also use "//&
                          "EPBL_MSTAR_SCHEME to specify the scheme for mstar.")
  elseif (mstar_mode == 2) then
    mstar_scheme = OM4_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = OM4 instead of the archaic MSTAR_MODE = 2.")
  elseif (mstar_mode == 3) then
    mstar_scheme = RH18_STRING
    call MOM_error(WARNING, "Use EPBL_MSTAR_SCHEME = REICHL_H18 instead of the archaic MSTAR_MODE = 3.")
  elseif (mstar_mode > 3) then
    call MOM_error(FATAL, "An unrecognized value of the obsolete parameter MSTAR_MODE was specified.")
  endif
  call log_param(param_file, mdl, "EPBL_MSTAR_SCHEME", mstar_scheme, &
                 "EPBL_MSTAR_SCHEME selects the method for setting mstar.  Valid values are: \n"//&
                 "\t CONSTANT   - Use a fixed mstar given by MSTAR \n"//&
                 "\t OM4        - Use L_Ekman/L_Obukhov in the stabilizing limit, as in OM4 \n"//&
                 "\t REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.", &
                 default=CONSTANT_STRING)
  mstar_scheme = uppercase(mstar_scheme)
  select case (mstar_scheme)
    case (CONSTANT_STRING)
      CS%mstar_scheme = Use_Fixed_mstar
    case (OM4_STRING)
      CS%mstar_scheme = mstar_from_Ekman
    case (RH18_STRING)
      CS%mstar_scheme = mstar_from_RH18
    case default
      call MOM_mesg('energetic_PBL_init: EPBL_MSTAR_SCHEME ="'//trim(mstar_scheme)//'"', 0)
      call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
            "EPBL_MSTAR_SCHEME = "//trim(mstar_scheme)//" found in input file.")
  end select
  call get_param(param_file, mdl, "MSTAR", CS%fixed_mstar, &
                 "The ratio of the friction velocity cubed to the TKE input to the "//&
                 "surface boundary layer.  This option is used if EPBL_MSTAR_SCHEME = CONSTANT.", &
                 units="nondim", default=1.2, do_not_log=(CS%mstar_scheme/=Use_Fixed_mstar))

  call get_param(param_file, mdl, "MSTAR_CAP", CS%mstar_cap, &
                 "If this value is positive, it sets the maximum value of mstar "//&
                 "allowed in ePBL.  (This is not used if EPBL_mstar_scheme = CONSTANT).", &
                 units="nondim", default=-1.0, do_not_log=(CS%mstar_scheme==Use_Fixed_mstar))
  ! mstar_scheme==mstar_from_Ekman options
  call get_param(param_file, mdl, "MSTAR2_COEF1", CS%mstar_coef, &
                 "Coefficient in computing mstar when rotation and stabilizing "//&
                 "effects are both important (used if EPBL_mstar_scheme = OM4).", &
                 units="nondim", default=0.3, do_not_log=(CS%mstar_scheme/=mstar_from_Ekman))
  call get_param(param_file, mdl, "MSTAR2_COEF2", CS%C_Ek, &
                 "Coefficient in computing mstar when only rotation limits "// &
                 "the total mixing (used if EPBL_MSTAR_SCHEME = OM4)", &
                 units="nondim", default=0.085, do_not_log=(CS%mstar_scheme/=mstar_from_Ekman))
  ! mstar_scheme==mstar_from_RH18 options
  call get_param(param_file, mdl, "RH18_MSTAR_CN1", CS%RH18_mstar_cn1,&
                 "MSTAR_N coefficient 1 (outer-most coefficient for fit). "//&
                 "The value of 0.275 is given in RH18.  Increasing this "//&
                 "coefficient increases mstar for all values of Hf/ust, but more "//&
                 "effectively at low values (weakly developed OSBLs).", &
                 units="nondim", default=0.275, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CN2", CS%RH18_mstar_cn2,&
                 "MSTAR_N coefficient 2 (coefficient outside of exponential decay). "//&
                 "The value of 8.0 is given in RH18.  Increasing this coefficient "//&
                 "increases mstar for all values of HF/ust, with a much more even "//&
                 "effect across a wide range of Hf/ust than CN1.", &
                 units="nondim", default=8.0, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CN3", CS%RH18_mstar_CN3,&
                 "MSTAR_N coefficient 3 (exponential decay coefficient). "//&
                 "The value of -5.0 is given in RH18.  Increasing this increases how "//&
                 "quickly the value of mstar decreases as Hf/ust increases.", &
                  units="nondim", default=-5.0, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CS1", CS%RH18_mstar_cs1,&
                 "MSTAR_S coefficient for RH18 in stabilizing limit. "//&
                 "The value of 0.2 is given in RH18 and increasing it increases "//&
                 "mstar in the presence of a stabilizing surface buoyancy flux.", &
                 units="nondim", default=0.2, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
  call get_param(param_file, mdl, "RH18_MSTAR_CS2", CS%RH18_mstar_cs2,&
                 "MSTAR_S exponent for RH18 in stabilizing limit. "//&
                 "The value of 0.4 is given in RH18 and increasing it increases mstar "//&
                 "exponentially in the presence of a stabilizing surface buoyancy flux.", &
                 Units="nondim", default=0.4, do_not_log=(CS%mstar_scheme/=mstar_from_RH18))
!/ BBL mstar related options
  call get_param(param_file, mdl, "EPBL_BBL_USE_MSTAR", CS%ePBL_BBL_use_mstar, &
                 "A logical to use mstar in the calculation of TKE in the ePBL BBL scheme", &
                 units="nondim", default=.false.)
  if (CS%ePBL_BBL_use_mstar) then
    call get_param(param_file, mdl, "EPBL_BBL_MSTAR_SCHEME", tmpstr, &
                   "EPBL_BBL_MSTAR_SCHEME selects the method for setting mstar in the BBL.  Valid values are: \n"//&
                   "\t CONSTANT   - Use a fixed mstar given by MSTAR_BBL \n"//&
                   "\t OM4        - Use L_Ekman/L_Obukhov in the stabilizing limit, as in OM4 \n"//&
                   "\t REICHL_H18 - Use the scheme documented in Reichl & Hallberg, 2018.", &
                   default=mstar_scheme)
    tmpstr = uppercase(tmpstr)
    select case (tmpstr)
      case (CONSTANT_STRING)
        CS%BBL_mstar_scheme = Use_Fixed_mstar
      case (OM4_STRING)
        CS%BBL_mstar_scheme = mstar_from_Ekman
      case (RH18_STRING)
        CS%BBL_mstar_scheme = mstar_from_RH18
      case default
        call MOM_mesg('energetic_PBL_init: EPBL_BBL_MSTAR_SCHEME ="'//trim(tmpstr)//'"', 0)
        call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
              "EPBL_BBL_MSTAR_SCHEME = "//trim(tmpstr)//" found in input file.")
    end select
    call get_param(param_file, mdl, "MSTAR_BBL", CS%BBL_fixed_mstar, &
                   "The ratio of the friction velocity cubed to the TKE input to the "//&
                   "bottom boundary layer.  This option is used if EPBL_BBL_MSTAR_SCHEME = CONSTANT.", &
                   units="nondim", default=1.2, do_not_log=(CS%BBL_mstar_scheme/=Use_Fixed_mstar))
  endif

!/ Convective turbulence related options
  call get_param(param_file, mdl, "NSTAR", CS%nstar, &
                 "The portion of the buoyant potential energy imparted by "//&
                 "surface fluxes that is available to drive entrainment "//&
                 "at the base of mixed layer when that energy is positive.", &
                 units="nondim", default=0.2)
  call get_param(param_file, mdl, "MSTAR_CONV_ADJ", CS%mstar_convect_coef, &
                 "Coefficient used for reducing mstar during convection "//&
                 "due to reduction of stable density gradient.", &
                 units="nondim", default=0.0)

!/ Mixing Length Options
  call get_param(param_file, mdl, "USE_MLD_ITERATION", CS%Use_MLD_iteration, &
                 "A logical that specifies whether or not to use the "//&
                 "distance to the bottom of the actively turbulent boundary "//&
                 "layer to help set the EPBL length scale.", default=.true.)
  call get_param(param_file, mdl, "EPBL_TRANSITION_SCALE", CS%transLay_scale, &
                 "A scale for the mixing length in the transition layer "//&
                 "at the edge of the boundary layer as a fraction of the "//&
                 "boundary layer thickness.", units="nondim", default=0.1)
  if ( CS%Use_MLD_iteration .and. abs(CS%transLay_scale-0.5) >= 0.5) then
    call MOM_error(FATAL, "If flag USE_MLD_ITERATION is true, then "//&
                 "EPBL_TRANSITION should be greater than 0 and less than 1.")
  endif

  call get_param(param_file, mdl, "MLD_ITERATION_GUESS", CS%MLD_ITERATION_GUESS, &
                 "If true, use the previous timestep MLD as a first guess in the MLD iteration, "//&
                 "otherwise use half the ocean depth as the first guess of the boundary layer "//&
                 "depth.  The default is false to facilitate reproducibility.", &
                 default=.false., do_not_log=.not.CS%Use_MLD_iteration)
  call get_param(param_file, mdl, "EPBL_MLD_TOLERANCE", CS%MLD_tol, &
                 "The tolerance for the iteratively determined mixed "//&
                 "layer depth.  This is only used with USE_MLD_ITERATION.", &
                 units="meter", default=1.0, scale=US%m_to_Z, do_not_log=.not.CS%Use_MLD_iteration)
  call get_param(param_file, mdl, "EPBL_MLD_BISECTION", CS%MLD_bisection, &
                 "If true, use bisection with the iterative determination of the self-consistent "//&
                 "mixed layer depth.  Otherwise use the false position after a maximum and minimum "//&
                 "bound have been evaluated and the returned value or bisection before this.", &
                 default=.false., do_not_log=.not.CS%Use_MLD_iteration)
   call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
   call get_param(param_file, mdl, "EPBL_MLD_ITER_BUG", CS%MLD_iter_bug, &
                 "If true, use buggy logic that gives the wrong bounds for the next iteration "//&
                 "when successive guesses increase by exactly EPBL_MLD_TOLERANCE.", &
                 default=enable_bugs, do_not_log=.not.CS%Use_MLD_iteration)
  call get_param(param_file, mdl, "EPBL_MLD_MAX_ITS", CS%max_MLD_its, &
                 "The maximum number of iterations that can be used to find a self-consistent "//&
                 "mixed layer depth.  If EPBL_MLD_BISECTION is true, the maximum number "//&
                 "of iterations needed is set by Depth/2^MAX_ITS < EPBL_MLD_TOLERANCE.", &
                 default=20, do_not_log=.not.CS%Use_MLD_iteration)
  if (.not.CS%Use_MLD_iteration) CS%Max_MLD_Its = 1
  call get_param(param_file, mdl, "EPBL_MIN_MIX_LEN", CS%min_mix_len, &
                 "The minimum mixing length scale that will be used "//&
                 "by ePBL.  The default (0) does not set a minimum.", &
                 units="meter", default=0.0, scale=US%m_to_Z)

  call get_param(param_file, mdl, "MIX_LEN_EXPONENT", CS%MixLenExponent, &
                 "The exponent applied to the ratio of the distance to the MLD "//&
                 "and the MLD depth which determines the shape of the mixing length. "//&
                 "This is only used if USE_MLD_ITERATION is True.", &
                 units="nondim", default=2.0)

!/ Turbulent velocity scale in mixing coefficient
  call get_param(param_file, mdl, "EPBL_VEL_SCALE_SCHEME", vel_scale_str, &
                 "Selects the method for translating TKE into turbulent velocities. "//&
                 "Valid values are: \n"//&
                 "\t CUBE_ROOT_TKE  - A constant times the cube root of remaining TKE. \n"//&
                 "\t REICHL_H18 - Use the scheme based on a combination of w* and v* as \n"//&
                 "\t              documented in Reichl & Hallberg, 2018.", &
                 default=ROOT_TKE_STRING, do_not_log=.true.)
  call get_param(param_file, mdl, "EPBL_VEL_SCALE_MODE", wT_mode, default=-1)
  if (wT_mode == 0) then
    vel_scale_str = ROOT_TKE_STRING
    call MOM_error(WARNING, "Use EPBL_VEL_SCALE_SCHEME = CUBE_ROOT_TKE instead of the archaic EPBL_VEL_SCALE_MODE = 0.")
  elseif (wT_mode == 1) then
    vel_scale_str = RH18_STRING
    call MOM_error(WARNING, "Use EPBL_VEL_SCALE_SCHEME = REICHL_H18 instead of the archaic EPBL_VEL_SCALE_MODE = 1.")
  elseif (wT_mode >= 2) then
    call MOM_error(FATAL, "An unrecognized value of the obsolete parameter EPBL_VEL_SCALE_MODE was specified.")
  endif
  call log_param(param_file, mdl, "EPBL_VEL_SCALE_SCHEME", vel_scale_str, &
                 "Selects the method for translating TKE into turbulent velocities. "//&
                 "Valid values are: \n"//&
                 "\t CUBE_ROOT_TKE  - A constant times the cube root of remaining TKE. \n"//&
                 "\t REICHL_H18 - Use the scheme based on a combination of w* and v* as \n"//&
                 "\t              documented in Reichl & Hallberg, 2018.", &
                 default=ROOT_TKE_STRING)
  vel_scale_str = uppercase(vel_scale_str)
  select case (vel_scale_str)
    case (ROOT_TKE_STRING)
      CS%wT_scheme = wT_from_cRoot_TKE
    case (RH18_STRING)
      CS%wT_scheme = wT_from_RH18
    case default
      call MOM_mesg('energetic_PBL_init: EPBL_VEL_SCALE_SCHEME ="'//trim(vel_scale_str)//'"', 0)
      call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
            "EPBL_VEL_SCALE_SCHEME = "//trim(vel_scale_str)//" found in input file.")
  end select

  call get_param(param_file, mdl, "WSTAR_USTAR_COEF", CS%wstar_ustar_coef, &
                 "A ratio relating the efficiency with which convectively "//&
                 "released energy is converted to a turbulent velocity, "//&
                 "relative to mechanically forced TKE. Making this larger "//&
                 "increases the BL diffusivity", units="nondim", default=1.0)
  call get_param(param_file, mdl, "EPBL_VEL_SCALE_FACTOR", CS%vstar_scale_fac, &
                 "An overall nondimensional scaling factor for wT. "//&
                 "Making this larger increases the PBL diffusivity.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "VSTAR_SURF_FAC", CS%vstar_surf_fac,&
                 "The proportionality times ustar to set vstar at the surface.", &
                 units="nondim", default=1.2)

  !/ Bottom boundary layer mixing related options
  call get_param(param_file, mdl, "EPBL_BBL_EFFIC", CS%ePBL_BBL_effic, &
                 "The efficiency of bottom boundary layer mixing via ePBL.  Setting this to a "//&
                 "value that is greater than 0 to enable bottom boundary layer mixing from EPBL.", &
                 units="nondim", default=0.0, scale=US%L_to_Z**2)
  call get_param(param_file, mdl, "EPBL_BBL_TIDAL_EFFIC", CS%ePBL_tidal_effic, &
                 "The efficiency of bottom boundary layer mixing via ePBL driven by the "//&
                 "bottom drag dissipation of tides, as provided in fluxes%BBL_tidal_dis.", &
                 units="nondim", default=0.0, scale=US%L_to_Z**2) !### Change the default to follow EPBL_BBL_EFFIC?
  no_BBL = ((CS%ePBL_BBL_effic <= 0.0) .and. (CS%ePBL_tidal_effic <= 0.0))

  call get_param(param_file, mdl, "USE_BBLD_ITERATION", CS%Use_BBLD_iteration, &
                 "A logical that specifies whether or not to use the distance to the top of the "//&
                 "actively turbulent bottom boundary layer to help set the EPBL length scale.", &
                 default=.true., do_not_log=no_BBL)
  call get_param(param_file, mdl, "TKE_DECAY_BBL", CS%TKE_decay_BBL, &
                 "TKE_DECAY_BBL relates the vertical rate of decay of the TKE available for "//&
                 "mechanical entrainment in the bottom boundary layer to the natural Ekman depth.", &
                 units="nondim", default=CS%TKE_decay, do_not_log=no_BBL)
  call get_param(param_file, mdl, "MIX_LEN_EXPONENT_BBL", CS%MixLenExponent_BBL, &
                 "The exponent applied to the ratio of the distance to the top of the BBL "//&
                 "and the total BBL depth which determines the shape of the mixing length. "//&
                 "This is only used if USE_MLD_ITERATION is True.", &
                 units="nondim", default=2.0, do_not_log=(no_BBL.or.(.not.CS%Use_BBLD_iteration)))
  call get_param(param_file, mdl, "EPBL_MIN_BBL_MIX_LEN", CS%min_BBL_mix_len, &
                 "The minimum mixing length scale that will be used by ePBL for bottom boundary "//&
                 "layer mixing.  Choosing (0) does not set a minimum.", &
                 units="meter", default=CS%min_mix_len, scale=US%m_to_Z, do_not_log=no_BBL)
  call get_param(param_file, mdl, "EPBL_BBLD_TOLERANCE", CS%BBLD_tol, &
                 "The tolerance for the iteratively determined bottom boundary layer depth.  "//&
                 "This is only used with USE_MLD_ITERATION.", &
                 units="meter", default=US%Z_to_m*CS%MLD_tol, scale=US%m_to_Z, &
                 do_not_log=(no_BBL.or.(.not.CS%Use_MLD_iteration)))
  call get_param(param_file, mdl, "EPBL_BBLD_MAX_ITS", CS%max_BBLD_its, &
                 "The maximum number of iterations that can be used to find a self-consistent "//&
                 "bottom boundary layer depth.", &
                 default=CS%max_MLD_its, do_not_log=(no_BBL.or.(.not.CS%Use_MLD_iteration)))
  if (.not.CS%Use_MLD_iteration) CS%max_BBLD_its = 1

  call get_param(param_file, mdl, "EPBL_BBL_VEL_SCALE_SCHEME", tmpstr, &
                 "Selects the method for translating bottom boundary layer TKE into turbulent velocities. "//&
                 "Valid values are: \n"//&
                 "\t CUBE_ROOT_TKE  - A constant times the cube root of remaining BBL TKE. \n"//&
                 "\t REICHL_H18 - Use the scheme based on a combination of w* and v* as \n"//&
                 "\t              documented in Reichl & Hallberg, 2018.", &
                 default=vel_scale_str, do_not_log=no_BBL)
  select case (tmpstr)
    case (ROOT_TKE_STRING)
      CS%wT_scheme_BBL = wT_from_cRoot_TKE
    case (RH18_STRING)
      CS%wT_scheme_BBL = wT_from_RH18
    case default
      call MOM_mesg('energetic_PBL_init: EPBL_BBL_VEL_SCALE_SCHEME ="'//trim(tmpstr)//'"', 0)
      call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
            "EPBL_BBL_VEL_SCALE_SCHEME = "//trim(tmpstr)//" found in input file.")
  end select
  call get_param(param_file, mdl, "EPBL_BBL_VEL_SCALE_FACTOR", CS%vstar_scale_fac_BBL, &
                 "An overall nondimensional scaling factor for wT in the bottom boundary layer. "//&
                 "Making this larger increases the bottom boundary layer diffusivity.", &
                 units="nondim", default=CS%vstar_scale_fac, do_not_log=no_BBL)
  call get_param(param_file, mdl, "VSTAR_BBL_SURF_FAC", CS%vstar_surf_fac_BBL,&
                 "The proportionality times ustar to set vstar in the bottom boundary layer.", &
                 units="nondim", default=CS%vstar_surf_fac, do_not_log=(no_BBL.or.(CS%wT_scheme_BBL/=wT_from_RH18)))
  call get_param(param_file, mdl, "EKMAN_SCALE_COEF_BBL", CS%Ekman_scale_coef_BBL, &
                 "A nondimensional scaling factor controlling the inhibition of the diffusive "//&
                 "length scale by rotation in the bottom boundary layer.  Making this larger "//&
                 "decreases the bottom boundary layer diffusivity.", &
                 units="nondim", default=CS%Ekman_scale_coef, do_not_log=no_BBL)
  call get_param(param_file, mdl, "EPBL_BBL_EFFIC_BUG", CS%BBL_effic_bug, &
                 "If true, overestimate the efficiency of the non-tidal ePBL bottom boundary "//&
                 "layer diffusivity by a factor of 1/sqrt(CDRAG), which is often a factor of "//&
                 "about 18.3.", default=.false., do_not_log=(CS%ePBL_BBL_effic<=0.0))

  call get_param(param_file, mdl, "DECAY_ADJUSTED_BBL_TKE", CS%decay_adjusted_BBL_TKE, &
                 "If true, include an adjustment factor in the bottom boundary layer energetics "//&
                 "that accounts for an exponential decay of TKE from a near-bottom source and "//&
                 "an assumed piecewise linear profile of the buoyancy flux response to a change "//&
                 "in a diffusivity.", &
                 default=.false., do_not_log=no_BBL)

  !/ Options related to Langmuir turbulence
  call get_param(param_file, mdl, "USE_LA_LI2016", use_LA_Windsea, &
       "A logical to use the Li et al. 2016 (submitted) formula to "//&
       "determine the Langmuir number.", default=.false.)
  ! Note this can be activated in other ways, but this preserves the old method.
  if (use_LA_windsea) then
    CS%use_LT = .true.
  else
    call get_param(param_file, mdl, "EPBL_LT", CS%use_LT, &
                 "A logical to use a LT parameterization.", default=.false.)
  endif
  if (CS%use_LT) then
    call get_param(param_file, mdl, "EPBL_LANGMUIR_SCHEME", tmpstr, &
                 "EPBL_LANGMUIR_SCHEME selects the method for including Langmuir turbulence. "//&
                 "Valid values are: \n"//&
                 "\t NONE     - Do not do any extra mixing due to Langmuir turbulence \n"//&
                 "\t RESCALE  - Use a multiplicative rescaling of mstar to account for Langmuir turbulence \n"//&
                 "\t ADDITIVE - Add a Langmuir turbulence contribution to mstar to other contributions", &
                 default=NONE_STRING, do_not_log=.true.)
    call get_param(param_file, mdl, "LT_ENHANCE", LT_enhance, default=-1)
    if (LT_enhance == 0) then
      tmpstr = NONE_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = NONE instead of the archaic LT_ENHANCE = 0.")
    elseif (LT_enhance == 1) then
      call MOM_error(FATAL, "You are using a legacy LT_ENHANCE mode in ePBL that has been phased out. "//&
                            "If you need to use this setting please report this error.  Also use "//&
                            "EPBL_LANGMUIR_SCHEME to specify the scheme for mstar.")
    elseif (LT_enhance == 2) then
      tmpstr = RESCALED_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = RESCALE instead of the archaic LT_ENHANCE = 2.")
    elseif (LT_enhance == 3) then
      tmpstr = ADDITIVE_STRING
      call MOM_error(WARNING, "Use EPBL_LANGMUIR_SCHEME = ADDITIVE instead of the archaic LT_ENHANCE = 3.")
    elseif (LT_enhance > 3) then
      call MOM_error(FATAL, "An unrecognized value of the obsolete parameter LT_ENHANCE was specified.")
    endif
    call log_param(param_file, mdl, "EPBL_LANGMUIR_SCHEME", tmpstr, &
                 "EPBL_LANGMUIR_SCHEME selects the method for including Langmuir turbulence. "//&
                 "Valid values are: \n"//&
                 "\t NONE     - Do not do any extra mixing due to Langmuir turbulence \n"//&
                 "\t RESCALE  - Use a multiplicative rescaling of mstar to account for Langmuir turbulence \n"//&
                 "\t ADDITIVE - Add a Langmuir turbulence contribution to mstar to other contributions", &
                 default=NONE_STRING)
    tmpstr = uppercase(tmpstr)
    select case (tmpstr)
      case (NONE_STRING)
        CS%LT_enhance_form = No_Langmuir
      case (RESCALED_STRING)
        CS%LT_enhance_form = Langmuir_rescale
      case (ADDITIVE_STRING)
        CS%LT_enhance_form = Langmuir_add
      case default
        call MOM_mesg('energetic_PBL_init: EPBL_LANGMUIR_SCHEME ="'//trim(tmpstr)//'"', 0)
        call MOM_error(FATAL, "energetic_PBL_init: Unrecognized setting "// &
              "EPBL_LANGMUIR_SCHEME = "//trim(tmpstr)//" found in input file.")
    end select

    call get_param(param_file, mdl, "LT_ENHANCE_COEF", CS%LT_enhance_coef, &
                 "Coefficient for Langmuir enhancement of mstar", &
                 units="nondim", default=0.447, do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_ENHANCE_EXP", CS%LT_enhance_exp, &
                 "Exponent for Langmuir enhancement of mstar", &
                 units="nondim", default=-1.33,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC1", CS%LaC_MLD_Ek, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching Ekman depth.", &
                 units="nondim", default=-0.87,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC2", CS%LaC_MLD_Ob_stab, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching stable Obukhov depth.", &
                 units="nondim", default=0.0,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC3", CS%LaC_MLD_Ob_un, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "MLD approaching unstable Obukhov depth.", &
                 units="nondim", default=0.0,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC4", CS%Lac_Ek_Ob_stab, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "ratio of Ekman to stable Obukhov depth.", &
                 units="nondim", default=0.95,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
    call get_param(param_file, mdl, "LT_MOD_LAC5", CS%Lac_Ek_Ob_un, &
                 "Coefficient for modification of Langmuir number due to "//&
                 "ratio of Ekman to unstable Obukhov depth.", &
                 units="nondim", default=0.95,  do_not_log=(CS%LT_enhance_form==No_Langmuir))
  endif

  !/Options related to Machine Learning Equation Discovery
  ! Logial flags for using shape function from equation discovery - machine learning
  ! EPBL_EQD_DIFFUSIVITY : EPBL + Equation Discovery Diffusivity parameters

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_SHAPE", CS%eqdisc, &
                 "Logical flag for activating ML equation for shape function "// &
                 "that uses forcing to change its structure.", &
                 units="nondim", default=.false.)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_VELOCITY", CS%eqdisc_v0, &
                   "Logical flag for activating ML equation discovery for velocity scale", &
                   units="nondim", default=.false.)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_VELOCITY_H", CS%eqdisc_v0h, &
                   "Logical flag for activating ML equation discovery for velocity scale with h as input", &
                   units="nondim", default=.false.)


  ! sets a  lower cap for abs_f (Coriolis parameter) required in equation for v_0.
  ! Small value, solution not sensitive below 1 deg Latitute
  ! Default value of 2.5384E-07 corresponds to 0.1 deg.
  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_CORIOLIS_LOWER_CAP", CS%f_lower, &
                       "value of lower limit cap for v0, default is for 0.1 deg, insensitive below 1deg", &
                       units="s-1", default=2.5384E-07, scale=US%T_to_S, &
                       do_not_log=.not.CS%eqdisc_v0)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_V0_LOWER_CAP", CS%v0_lower_cap, &
                       "value of lower limit cap for Coriolis in v0", &
                       units="m s-1", default=0.0001, scale=US%m_to_Z*US%T_to_s, &
                       do_not_log=.not.(CS%eqdisc_v0.or.CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_V0_UPPER_CAP", CS%v0_upper_cap, &
                       "value of upper limit cap for Coriolis in v0", &
                       units="m s-1", default=0.1, scale=US%m_to_Z*US%T_to_s, &
                       do_not_log=.not.(CS%eqdisc_v0.or.CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_BFLUX_LOWER_CAP", CS%bflux_lower_cap, &
                       "value of lower limit cap for Bflux used in setting in v0", &
                       units="m2 s-3", default=-7.0E-07, scale=(US%m_to_L**2)*(US%T_to_s**3), &
                       do_not_log=.not.(CS%eqdisc_v0.or.CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_BFLUX_UPPER_CAP", CS%bflux_upper_cap, &
                       "value of upper limit cap for Bflux used in setting in v0", &
                       units="m2 s-3", default=7.0E-07, scale=(US%m_to_L**2)*(US%T_to_s**3), &
                       do_not_log=.not.(CS%eqdisc_v0.or.CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_SIGMA_MAX_LOWER_CAP", CS%sigma_max_lower_cap, &
                       "value of lower limit cap for sigma coordinate of maximum for diffusivity", &
                       units="nondim", default=0.1, do_not_log=.not.CS%eqdisc)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_SIGMA_MAX_UPPER_CAP", CS%sigma_max_upper_cap, &
                       "value of upper limit cap for sigma coordinate of maximum for diffusivity", &
                       units="nondim", default=0.7, do_not_log=.not.CS%eqdisc)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_EH_UPPER_CAP", CS%Eh_upper_cap, &
                       "value of upper limit cap for boundary layer depth by Ekman depth hf/u", &
                       units="nondim", default=2.0, do_not_log=.not.CS%eqdisc)

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_LH_CAP", CS%Lh_cap, &
                       "value of upper limit cap for boundary layer depth by Monin-Obukhov depth hB/u^3", &
                       units="nondim", default=8.0, do_not_log=.not.CS%eqdisc)

  ! The coefficients used for machine learned diffusivity
  ! c1 to c6 used for sigma_m,
  !  7 to 9 v_0 surface heating, 10 to 14 v_0 surface cooling (ML velocity scale without h as input)
  ! 14, 15, & 16 for v_0h surface heating, 17, 18, & 14 for v_0h surface cooling (ML velocity scale with h as input)
  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_COEFFS", CS%ML_c, &
                 "Coefficient used for ML diffusivity 1 to 18 ", units="nondim", &
                  defaults=(/1.7908 , 0.6904, 0.0712, 0.4380, 2.6821, 1.5845, 0.1550,  1.1120,  0.8616, 0.0984, &
                             45.0,    2.8570, 3.290,  0.0785, 0.650,  0.0944, 6.0277, 15.7292 /), &
                  do_not_log=.not.(CS%eqdisc .or. CS%eqdisc_v0 .or. CS%eqdisc_v0h))

  call get_param(param_file, mdl, "EPBL_EQD_DIFFUSIVITY_SHAPE_FUNCTION_EPSILON", CS%shape_function_epsilon, &
                 "Constant value of OSBL shape function below the boundary layer", &
                 units="nondim", default=0.01, do_not_log=.not.CS%eqdisc)

  !/ options end for Machine Learning Equation Discovery

  !/ Options for documenting differences from parameter choices
  call get_param(param_file, mdl, "EPBL_OPTIONS_DIFF", CS%options_diff, &
                 "If positive, this is a coded integer indicating a pair of settings whose "//&
                 "differences are diagnosed in a passive diagnostic mode  via extra calls to "//&
                 "ePBL_column.  If this is 0 or negative no extra calls occur.", &
                 default=0)
  if (CS%options_diff > 0) then
    if (CS%options_diff == 1) then
      diff_text = "EPBL_ORIGINAL_PE_CALC settings"
    elseif (CS%options_diff == 2) then
      diff_text = "EPBL_ANSWER_DATE settings"
    elseif (CS%options_diff == 3) then
      diff_text = "DIRECT_EPBL_MIXING_CALC settings"
    elseif (CS%options_diff == 4) then
      diff_text = "BBL DIRECT_EPBL_MIXING_CALC settings"
    elseif (CS%options_diff == 5) then
      diff_text = "BBL DECAY_ADJUSTED_BBL_TKE settings"
    else
      diff_text = "unchanged settings"
    endif
  endif

!/ Logging parameters
  ! This gives a minimum decay scale that is typically much less than Angstrom.
  CS%ustar_min = 2e-4*CS%omega*(GV%Angstrom_Z + GV%dZ_subroundoff)
  call log_param(param_file, mdl, "!EPBL_USTAR_MIN", CS%ustar_min, &
                 "The (tiny) minimum friction velocity used within the "//&
                 "ePBL code, derived from OMEGA and ANGSTROM.", &
                 units="m s-1", unscale=US%Z_to_m*US%s_to_T, &
                 like_default=.true.)


!/ Checking output flags
  CS%id_Kd_ePBL_col_by_col = register_diag_field('ocean_model', 'Kd_ePBL_col_by_col', diag%axesTi, Time, &
      'ePBL diapycnal diffusivity at interfaces posted column by column', 'm2 s-1', conversion=GV%HZ_T_to_m2_s)
  CS%id_ML_depth = register_diag_field('ocean_model', 'ePBL_h_ML', diag%axesT1, &
      Time, 'Surface boundary layer depth', units='m', conversion=US%Z_to_m, &
      cmor_long_name='Ocean Mixed Layer Thickness Defined by Mixing Scheme')
  ! This is an alias for the same variable as ePBL_h_ML
  CS%id_hML_depth = register_diag_field('ocean_model', 'h_ML', diag%axesT1, &
      Time, 'Surface mixed layer depth based on active turbulence', units='m', conversion=US%Z_to_m)
  CS%id_ustar_ePBL = register_diag_field('ocean_model', 'ePBL_ustar', diag%axesT1, &
      Time, 'Surface friction in ePBL', units='m s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_bflx_ePBL = register_diag_field('ocean_model', 'ePBL_bflx', diag%axesT1, &
      Time, 'Surface buoyancy flux in ePBL', units='m2 s-3', conversion=US%Z_to_m**2*US%s_to_T**3)
  CS%id_TKE_wind = register_diag_field('ocean_model', 'ePBL_TKE_wind', diag%axesT1, &
      Time, 'Wind-stirring source of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_MKE = register_diag_field('ocean_model', 'ePBL_TKE_MKE', diag%axesT1, &
      Time, 'Mean kinetic energy source of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_conv = register_diag_field('ocean_model', 'ePBL_TKE_conv', diag%axesT1, &
      Time, 'Convective source of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_forcing = register_diag_field('ocean_model', 'ePBL_TKE_forcing', diag%axesT1, &
      Time, 'TKE consumed by mixing surface forcing or penetrative shortwave radation '//&
            'through model layers', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_mixing = register_diag_field('ocean_model', 'ePBL_TKE_mixing', diag%axesT1, &
      Time, 'TKE consumed by mixing that deepens the mixed layer', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_mech_decay = register_diag_field('ocean_model', 'ePBL_TKE_mech_decay', diag%axesT1, &
      Time, 'Mechanical energy decay sink of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_TKE_conv_decay = register_diag_field('ocean_model', 'ePBL_TKE_conv_decay', diag%axesT1, &
      Time, 'Convective energy decay sink of mixed layer TKE', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_Mixing_Length = register_diag_field('ocean_model', 'Mixing_Length', diag%axesTi, &
      Time, 'Mixing Length that is used', units='m', conversion=US%Z_to_m)
  CS%id_Velocity_Scale = register_diag_field('ocean_model', 'Velocity_Scale', diag%axesTi, &
      Time, 'Velocity Scale that is used.', units='m s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_mstar_sfc = register_diag_field('ocean_model', 'MSTAR', diag%axesT1, &
      Time, 'Total mstar that is used.', 'nondim')
  if ((CS%ePBL_BBL_effic > 0.0) .or. (CS%ePBL_tidal_effic > 0.0) .or. CS%ePBL_BBL_use_mstar) then
    CS%id_Kd_BBL = register_diag_field('ocean_model', 'Kd_ePBL_BBL', diag%axesTi, &
        Time, 'ePBL bottom boundary layer diffusivity', units='m2 s-1', conversion=GV%HZ_T_to_m2_s)
    CS%id_BBL_Mix_Length = register_diag_field('ocean_model', 'BBL_Mixing_Length', diag%axesTi, &
        Time, 'ePBL bottom boundary layer mixing length', units='m', conversion=US%Z_to_m)
    CS%id_BBL_Vel_Scale = register_diag_field('ocean_model', 'BBL_Velocity_Scale', diag%axesTi, &
        Time, 'ePBL bottom boundary layer velocity scale', units='m s-1', conversion=US%Z_to_m*US%s_to_T)
    CS%id_BBL_depth = register_diag_field('ocean_model', 'h_BBL', diag%axesT1, &
        Time, 'Bottom boundary layer depth based on active turbulence', units='m', conversion=US%Z_to_m)
    CS%id_ustar_BBL = register_diag_field('ocean_model', 'ePBL_ustar_BBL', diag%axesT1, &
        Time, 'The bottom boundary layer friction velocity', units='m s-1', conversion=GV%H_to_m*US%s_to_T)
    CS%id_BBL_decay_scale = register_diag_field('ocean_model', 'BBL_decay_scale', diag%axesT1, &
        Time, 'The bottom boundary layer TKE decay lengthscale', units='m', conversion=GV%H_to_m)
    CS%id_TKE_BBL = register_diag_field('ocean_model', 'ePBL_BBL_TKE', diag%axesT1, &
        Time, 'The source of TKE for the bottom boundary layer', units='W m-2', conversion=US%RZ3_T3_to_W_m2)
    CS%id_TKE_BBL_mixing = register_diag_field('ocean_model', 'ePBL_BBL_TKE_mixing', diag%axesT1, &
        Time, 'TKE consumed by mixing that thickens the bottom boundary layer', &
        units='W m-2', conversion=US%RZ3_T3_to_W_m2)
    CS%id_TKE_BBL_decay = register_diag_field('ocean_model', 'ePBL_BBL_TKE_decay', diag%axesT1, &
        Time, 'Energy decay sink of mixed layer TKE in the bottom boundary layer', &
        units='W m-2', conversion=US%RZ3_T3_to_W_m2)
    CS%id_mstar_BBL = register_diag_field('ocean_model', 'MSTAR_BBL', diag%axesT1, &
        Time, 'Total BBL mstar that is used.', 'nondim')
  endif
  if (CS%use_LT) then
    CS%id_LA = register_diag_field('ocean_model', 'LA', diag%axesT1, &
        Time, 'Langmuir number.', 'nondim')
    CS%id_LA_mod = register_diag_field('ocean_model', 'LA_MOD', diag%axesT1, &
        Time, 'Modified Langmuir number.', 'nondim')
    CS%id_mstar_LT = register_diag_field('ocean_model', 'MSTAR_LT', diag%axesT1, &
        Time, 'Increase in mstar due to Langmuir Turbulence.', 'nondim')
  endif

  if (CS%options_diff > 0) then
    CS%id_opt_diff_Kd_ePBL = register_diag_field('ocean_model', 'ePBL_opt_diff_Kd_ePBL', diag%axesTi, &
        Time, 'Change in ePBL diapycnal diffusivity at interfaces due to '//trim(diff_text), &
        units='m2 s-1', conversion=GV%HZ_T_to_m2_s)
    CS%id_opt_maxdiff_Kd_ePBL = register_diag_field('ocean_model', 'ePBL_opt_maxdiff_Kd_ePBL', diag%axesT1, &
        Time, 'Column maximum change in ePBL diapycnal diffusivity at interfaces due to '//trim(diff_text), &
        units='m2 s-1', conversion=GV%HZ_T_to_m2_s)
    CS%id_opt_diff_hML_depth = register_diag_field('ocean_model', 'ePBL_opt_diff_h_ML', diag%axesT1, Time, &
        'Change in surface or bottom boundary layer depth based on active turbulence due to '//trim(diff_text), &
        units='m', conversion=US%Z_to_m)
  endif

  if (report_avg_its) then
    CS%sum_its(1) = real_to_EFP(0.0) ; CS%sum_its(2) = real_to_EFP(0.0)
    CS%sum_its_BBL(1) = real_to_EFP(0.0) ; CS%sum_its_BBL(2) = real_to_EFP(0.0)
  endif

  CS%TKE_diagnostics = (max(CS%id_TKE_wind, CS%id_TKE_MKE, CS%id_TKE_conv, &
                            CS%id_TKE_mixing, CS%id_TKE_mech_decay, CS%id_TKE_forcing, &
                            CS%id_TKE_conv_decay) > 0)
  if ((CS%ePBL_BBL_effic > 0.0) .or. (CS%ePBL_tidal_effic > 0.0) .or. CS%ePBL_BBL_use_mstar) then
    CS%TKE_diagnostics = CS%TKE_diagnostics .or. &
        (max(CS%id_TKE_BBL, CS%id_TKE_BBL_mixing, CS%id_TKE_BBL_decay) > 0)
  endif

  call safe_alloc_alloc(CS%ML_depth, isd, ied, jsd, jed)
  call safe_alloc_alloc(CS%BBL_depth, isd, ied, jsd, jed)

end subroutine energetic_PBL_init

!> Clean up and deallocate memory associated with the energetic_PBL module.
subroutine energetic_PBL_end(CS)
  type(energetic_PBL_CS), intent(inout) :: CS !< Energetic_PBL control structure

  character(len=256) :: mesg
  real :: avg_its ! The averaged number of iterations used by ePBL [nondim]

  if (allocated(CS%ML_depth))            deallocate(CS%ML_depth)
  if (allocated(CS%BBL_depth))           deallocate(CS%BBL_depth)

  if (report_avg_its) then
    call EFP_sum_across_PEs(CS%sum_its, 2)
    avg_its = EFP_to_real(CS%sum_its(1)) / EFP_to_real(CS%sum_its(2))
    write (mesg,*) "Average ePBL iterations = ", avg_its
    call MOM_mesg(mesg)

    if ((CS%ePBL_BBL_effic > 0.0) .or. (CS%ePBL_tidal_effic > 0.0) .or. CS%ePBL_BBL_use_mstar) then
      call EFP_sum_across_PEs(CS%sum_its_BBL, 2)
      avg_its = EFP_to_real(CS%sum_its_BBL(1)) / EFP_to_real(CS%sum_its_BBL(2))
      write (mesg,*) "Average ePBL BBL iterations = ", avg_its
      call MOM_mesg(mesg)
    endif
  endif
end subroutine energetic_PBL_end


end module MOM_energetic_PBL
