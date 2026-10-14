! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

#include "omp_loop_bind_compat.h"
#include <MOM_memory.h>
#include "do_concurrent_compat.h"

!> Routines for calculating baroclinic wave speeds
module MOM_wave_speed

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h, interpolate_column
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density_derivs, calculate_specific_vol_derivs

implicit none ; private

public wave_speed, wave_speeds, wave_speed_init, wave_speed_set_param

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for MOM_wave_speed
type, public :: wave_speed_CS ; private
  logical :: initialized = .false.     !< True if this control structure has been initialized.
  logical :: use_ebt_mode = .false.    !< If true, calculate the equivalent barotropic wave speed instead
                                       !! of the first baroclinic wave speed.
                                       !! This parameter controls the default behavior of wave_speed() which
                                       !! can be overridden by optional arguments.
  logical :: better_cg1_est = .false.  !< If true, use an improved estimate of the first mode
                                       !! internal wave speed.
  real :: mono_N2_column_fraction = 0. !< The lower fraction of water column over which N2 is limited as
                                       !! monotonic for the purposes of calculating the equivalent barotropic
                                       !! wave speed [nondim]. This parameter controls the default behavior of
                                       !! wave_speed() which can be overridden by optional arguments.
  real :: mono_N2_depth = -1.          !< The depth below which N2 is limited as monotonic for the purposes of
                                       !! calculating the equivalent barotropic wave speed [H ~> m or kg m-2].
                                       !! If this parameter is negative, this limiting does not occur.
                                       !! This parameter controls the default behavior of wave_speed() which
                                       !! can be overridden by optional arguments.
  real :: min_speed2 = 0.              !< The minimum mode 1 internal wave speed squared [L2 T-2 ~> m2 s-2]
  real :: wave_speed_tol = 0.001       !< The fractional tolerance with which to solve for the wave
                                       !! speeds [nondim]
  real :: c1_thresh = -1.0             !< A minimal value of the first mode internal wave speed
                                       !! below which all higher mode speeds are not calculated but
                                       !! are simply reported as 0 [L T-1 ~> m s-1].  A non-negative
                                       !! value must be specified via a call to wave_speed_init for
                                       !! the subroutine wave_speeds to be used (but not wave_speed).
  type(remapping_CS) :: remap_2018_CS  !< Used for vertical remapping when calculating equivalent barotropic
                                       !! mode structure for answer dates below 20190101.
  type(remapping_CS) :: remap_CS       !< Used for vertical remapping when calculating equivalent barotropic
                                       !! mode structure.
  integer :: remap_answer_date = 99991231 !< The vintage of the order of arithmetic and expressions to use
                                       !! for remapping.  Values below 20190101 recover the remapping
                                       !! answers from 2018, while higher values use more robust
                                       !! forms of the same remapping expressions.
  integer :: niblock = 0               !< The i block size used in wave_speed, or 0 to work on the
                                       !! full i-extent of the arrays at once [nondim].
  integer :: njblock = 0               !< The j block size used in wave_speed, or 0 to work on the
                                       !! full j-extent of the arrays at once [nondim].
  type(diag_ctrl), pointer :: diag     !< Diagnostics control structure
end type wave_speed_CS

contains

!> Calculates the wave speed of the first baroclinic mode.
subroutine wave_speed(h, tv, G, GV, US, cg1, CS, halo_size, use_ebt_mode, mono_N2_column_fraction, &
                      mono_N2_depth, modal_structure)
  type(ocean_grid_type),            intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),          intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),            intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)  :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)  :: tv !< Thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: cg1 !< First mode internal wave speed [L T-1 ~> m s-1]
  type(wave_speed_CS),              intent(in)  :: CS !< Wave speed control struct
  integer,                optional, intent(in)  :: halo_size !< Width of halo within which to
                                                       !! calculate wave speeds
  logical,                optional, intent(in)  :: use_ebt_mode !< If true, use the equivalent
                                          !! barotropic mode instead of the first baroclinic mode.
  real,                   optional, intent(in)  :: mono_N2_column_fraction !< The lower fraction
                                          !! of water column over which N2 is limited as monotonic
                                          !! for the purposes of calculating vertical modal structure [nondim].
  real,                   optional, intent(in)  :: mono_N2_depth !< A depth below which N2 is limited as
                                          !! monotonic for the purposes of calculating vertical
                                          !! modal structure [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                          optional, intent(out) :: modal_structure !< Normalized model structure [nondim]

  ! Local variables
  real, dimension(merge(G%ied-G%isd+1, CS%niblock, CS%niblock==0), &
                  merge(G%jed-G%jsd+1, CS%njblock, CS%njblock==0), SZK_(GV)+1) :: &
    dRho_dT, &    ! Partial derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
    dRho_dS, &    ! Partial derivative of density with salinity [R S-1 ~> kg m-3 ppt-1]
    dSpV_dT, &    ! Partial derivative of specific volume with temperature [R-1 C-1 ~> m3 kg-1 degC-1]
    dSpV_dS, &    ! Partial derivative of specific volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1]
    pres, &       ! Interface pressure [R L2 T-2 ~> Pa]
    T_int, &      ! Temperature interpolated to interfaces [C ~> degC]
    S_int, &      ! Salinity interpolated to interfaces [S ~> ppt]
    gprime        ! The reduced gravity across each interface [L2 H-1 T-2 ~> m s-2 or m4 s-2 kg-1].
  real, dimension(SZK_(GV)) :: &
    Igl, Igu      ! The inverse of the reduced gravity across an interface times
                  ! the thickness of the layer below (Igl) or above (Igu) it, in [T2 L-2 ~> s2 m-2].
  real, dimension(SZK_(GV)+1) :: &
    H_top, &      ! The distance of each filtered (or merged) interface from the ocean surface,
                  ! recomputed as a column-local temporary each time it is needed
                  ! [H ~> m or kg m-2]
    H_bot         ! The distance of each filtered (or merged) interface from the bottom,
                  ! recomputed as a column-local temporary each time it is needed
                  ! [H ~> m or kg m-2]
  real, dimension(merge(G%ied-G%isd+1, CS%niblock, CS%niblock==0), &
                  merge(G%jed-G%jsd+1, CS%njblock, CS%njblock==0), SZK_(GV)+1) :: &
    Hf, &         ! Layer thicknesses after very thin layers are combined [H ~> m or kg m-2]
    Tf, &         ! Layer temperatures after very thin layers are combined [C ~> degC]
    Sf, &         ! Layer salinities after very thin layers are combined [S ~> ppt]
    Rf            ! Layer densities after very thin layers are combined [R ~> kg m-3]
  real, dimension(merge(G%ied-G%isd+1, CS%niblock, CS%niblock==0), &
                  merge(G%jed-G%jsd+1, CS%njblock, CS%njblock==0), SZK_(GV)) :: &
    Hc, &         ! A column of layer thicknesses after convective instabilities are removed [H ~> m or kg m-2]
    Tc, &         ! A column of layer temperatures after convective instabilities are removed [C ~> degC]
    Sc, &         ! A column of layer salinities after convective instabilities are removed [S ~> ppt]
    Rc            ! A column of layer densities after convective instabilities are removed [R ~> kg m-3]
  real, dimension(merge(G%ied-G%isd+1, CS%niblock, CS%niblock==0), &
                  merge(G%jed-G%jsd+1, CS%njblock, CS%njblock==0), SZK_(GV)) :: &
    mode_struct   ! The mode structure [nondim], but it is also temporarily
                  ! in units of [L2 T-2 ~> m2 s-2] after it is modified inside of tdma6_block.
  real, dimension(SZK_(GV)) :: &
    I_beta, &     ! A column-local temporary variable in the tridiagonal solve [L2 T-2 ~> m2 s-2]
    yy            ! A column-local temporary variable in the tridiagonal solve with the same
                  ! units as mode_struct on entry to the solve [nondim]
  integer, dimension(merge(G%ied-G%isd+1, CS%niblock, CS%niblock==0), &
                     merge(G%jed-G%jsd+1, CS%njblock, CS%njblock==0)) :: &
    nkc           ! The number of layers in each column after merging, or 0 if no wave speed
                  ! was calculated there.
  integer :: EOSdom(3,2)  ! The i-, j- and k-index bounds for the batched equation of state call
  real :: I_Htot  ! The inverse of the total filtered thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: det, ddet ! Determinant of the eigen system and its derivative with lam.  Because the
                  ! units of the eigenvalue change with the number of layers and because of the
                  ! dynamic rescaling that is used to keep det in a numerically representable range,
                  ! the units of of det are hard to interpret, but det/ddet is always in units
                  ! of [T2 L-2 ~> s2 m-2]
  real :: lam     ! The eigenvalue [T2 L-2 ~> s2 m-2]
  real :: dlam    ! The change in estimates of the eigenvalue [T2 L-2 ~> s2 m-2]
  real :: lam0    ! The first guess of the eigenvalue [T2 L-2 ~> s2 m-2]
  real :: H_to_pres  ! A conversion factor from thicknesses to pressure [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1]
  real, dimension(merge(G%ied-G%isd+1, CS%niblock, CS%niblock==0), &
                  merge(G%jed-G%jsd+1, CS%njblock, CS%njblock==0)) :: &
    htot           ! Thicknesses [H ~> m or kg m-2]
  real :: hmin     ! A thickness [H ~> m or kg m-2]
  real :: H_here   ! A thickness [H ~> m or kg m-2]
  real :: HxT_here ! A layer integrated temperature [C H ~> degC m or degC kg m-2]
  real :: HxS_here ! A layer integrated salinity [S H ~> ppt m or ppt kg m-2]
  real :: HxR_here ! A layer integrated density [R H ~> kg m-2 or kg2 m-5]
  real :: speed2_tot ! overestimate of the mode-1 speed squared [L2 T-2 ~> m2 s-2]
  real :: cg1_min2 ! A floor in the squared first mode speed below which 0 is returned [L2 T-2 ~> m2 s-2]
  real :: cg1_est  ! An initial estimate of the squared first mode speed [L2 T-2 ~> m2 s-2]
  real :: I_Hnew   ! The inverse of a new layer thickness [H-1 ~> m-1 or m2 kg-1]
  real :: drxh_sum ! The sum of density differences across interfaces times thicknesses [R H ~> kg m-2 or kg2 m-5]
  real :: dSpVxh_sum ! The sum of specific volume differences across interfaces times
                   ! thicknesses [H R-1 ~> m4 kg-1 or m], negative for stable stratification.
  real :: g_Rho0   ! G_Earth/Rho0 [L2 T-2 H-1 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2].
  real :: c2_scale ! A scaling factor for wave speeds to help control the growth of the determinant and
                   ! its derivative with lam between rows of the Thomas algorithm solver [L2 s2 T-2 m-2 ~> nondim].
                   ! The exact value should not matter for the final result if it is an even power of 2.
  real :: tol_Hfrac ! Layers that together are smaller than this fraction of
                    ! the total water column can be merged for efficiency [nondim].
  real :: min_h_frac ! tol_Hfrac divided by the total number of layers [nondim].
  real :: tol_solve ! The fractional tolerance with which to solve for the wave speeds [nondim]
  real :: tol_merge ! The fractional change in estimated wave speed that is allowed
                    ! when deciding to merge layers in the calculation [nondim]
  real :: rescale   ! A rescaling factor to control the magnitude of the determinant [nondim]
  real :: I_rescale ! The reciprocal of the rescaling factor to control the magnitude of the determinant [nondim]
  integer, dimension(merge(G%ied-G%isd+1, CS%niblock, CS%niblock==0), &
                     merge(G%jed-G%jsd+1, CS%njblock, CS%njblock==0)) :: &
    kf            ! The number of active layers after filtering.
  integer :: max_kf     ! The largest number of filtered layers in any column of the current
                        ! block, used to limit the vertical extent of the equation of state call.
  integer, parameter :: max_itt = 10
  logical :: use_EOS    ! If true, density or specific volume is calculated from T & S using an equation of state.
  logical :: nonBous    ! If true, do not make the Boussinesq approximation.
  logical :: better_est ! If true, use an improved estimate of the first mode internal wave speed.
  logical :: do_merge   ! If true, merge the current layer with the one above.
  integer :: kc         ! The number of layers in the column after merging
  integer :: i, j, k, k2, itt, is, ie, js, je, nz, halo
  integer :: nii, njj ! The declared horizontal extents of the blocked working arrays, which
                      ! tdma6_block needs in order to index mode_struct the same way that it is indexed
                      ! here.  These are also the strides of the block loops.
  integer :: isb, ieb ! The first and last i-indices of the current block.
  integer :: jsb, jeb ! The first and last j-indices of the current block.
  integer :: iie, jje ! The number of columns actually in the current block, which is smaller
                      ! than nii or njj for the last block of a row or column.
  integer :: ii, jj ! Block-local 1-based i- and j-indices.
  real :: hw      ! The mean of the adjacent layer thicknesses [H ~> m or kg m-2]
  real :: sum_hc  ! The sum of the layer thicknesses [H ~> m or kg m-2]
  real :: gp      ! A limited local copy of gprime [L2 H-1 T-2 ~> m s-2 or m4 s-2 kg-1]
  real :: N2min   ! A minimum buoyancy frequency, including a slope rescaling factor [L2 H-2 T-2 ~> s-2 or m6 kg-2 s-2]
  logical :: below_mono_N2_frac  ! True if an interface is below the fractional depth where N2 should not increase.
  logical :: below_mono_N2_depth ! True if an interface is below the absolute depth where N2 should not increase.
  logical :: l_use_ebt_mode, calc_modal_structure
  real :: l_mono_N2_column_fraction ! A local value of mono_N2_column_fraction [nondim]
  real :: l_mono_N2_depth ! A local value of mono_N2_column_depth [H ~> m or kg m-2]
  real :: ms_min, ms_max ! The minimum and maximum mode structure values returned from tdma6 [L2 T-2 ~> m2 s-2]
  real :: ms_sq          ! The sum of the square of the values returned from tdma6 [L4 T-4 ~> m4 s-4]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke ; halo = 0

  nii = CS%niblock ; if (nii == 0) nii = G%ied - G%isd + 1
  njj = CS%njblock ; if (njj == 0) njj = G%jed - G%jsd + 1

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_wave_speed / wave_speed: "// &
           "Module must be initialized before it is used.")

  if (present(halo_size)) then
    halo = halo_size
    is = G%isc - halo ; ie = G%iec + halo ; js = G%jsc - halo ; je = G%jec + halo
  endif

  !$omp target enter data map(alloc: Hf, Tf, Sf, Rf, pres, T_int, S_int, dRho_dT, dRho_dS)
  !$omp target enter data map(alloc: dSpV_dT, dSpV_dS, H_top, H_bot, gprime, Igl, Igu)
  !$omp target enter data map(alloc: Hc, Tc, Sc, Rc, mode_struct)
  !$omp target enter data map(alloc: htot, kf, nkc)

  l_use_ebt_mode = CS%use_ebt_mode
  if (present(use_ebt_mode)) l_use_ebt_mode = use_ebt_mode
  l_mono_N2_column_fraction = CS%mono_N2_column_fraction
  if (present(mono_N2_column_fraction)) l_mono_N2_column_fraction = mono_N2_column_fraction
  l_mono_N2_depth = CS%mono_N2_depth
  if (present(mono_N2_depth)) l_mono_N2_depth = mono_N2_depth
  calc_modal_structure = l_use_ebt_mode
  if (present(modal_structure)) calc_modal_structure = .true.
  if (calc_modal_structure) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      modal_structure(i,j,k) = 0.0
    enddo ; enddo ; enddo
  endif

  nonBous = .not.(GV%Boussinesq .or. GV%semi_Boussinesq)
  H_to_pres = GV%H_to_RZ * GV%g_Earth
  ! Note that g_Rho0 = H_to_pres / GV%Rho0**2
  if (.not.nonBous) g_Rho0 = GV%g_Earth*GV%H_to_Z / GV%Rho0
  use_EOS = associated(tv%eqn_of_state)

  better_est = CS%better_cg1_est

  if (better_est) then
    tol_solve = CS%wave_speed_tol
    tol_Hfrac  = 0.1*tol_solve ; tol_merge = tol_solve / real(nz)
  else
    tol_solve = 0.001 ; tol_Hfrac  = 0.0001 ; tol_merge = 0.001
  endif

  ! The rescaling below can control the growth of the determinant provided that
  ! (tol_merge*cg1_min2/c2_scale > I_rescale).  For default values, this suggests a stable lower
  ! bound on min_speed of sqrt(nz/(tol_solve*rescale)) or 3e2/1024**2 = 2.9e-4 m/s for 90 layers.
  ! The upper bound on the rate of increase in the determinant is g'H/c2_scale < rescale or in the
  ! worst possible oceanic case of g'H < 0.5*10m/s2*1e4m = 5.e4 m2/s2 < 1024**2*c2_scale, suggesting
  ! that c2_scale can safely be set to 1/(16*1024**2), which would decrease the stable floor on
  ! min_speed to ~6.9e-8 m/s for 90 layers or 2.33e-7 m/s for 1000 layers.
  cg1_min2 = CS%min_speed2
  rescale = 1024.0**4 ; I_rescale = 1.0/rescale
  c2_scale = US%m_s_to_L_T**2 / 4096.0**2 ! Other powers of 2 give identical results.

  min_h_frac = tol_Hfrac / real(nz)

  do jsb=js,je,njj ; do isb=is,ie,nii
    ieb = min(isb + nii - 1, ie) ; iie = ieb - isb + 1
    jeb = min(jsb + njj - 1, je) ; jje = jeb - jsb + 1
    max_kf = 0

    !   The loop below is skipped over land, but the equation of state is evaluated over the
    ! whole block, so seed its inputs there.  Everywhere else only the topmost filtered layer
    ! needs to be seeded, because the padding loop below fills in whatever the merging leaves
    ! short of the bottom of the array.
    do concurrent (jj=1:jje, ii=1:iie)
      Hf(ii,jj,1) = 0.0 ; Tf(ii,jj,1) = 0.0 ; Sf(ii,jj,1) = 0.0 ; Rf(ii,jj,1) = 0.0
    enddo
    do concurrent (K=1:nz+1, jj=1:jje, ii=1:iie, G%mask2dT(isb+ii-1,jsb+jj-1) <= 0.0)
      pres(ii,jj,K) = 0.0 ; T_int(ii,jj,K) = 0.0 ; S_int(ii,jj,K) = 0.0
    enddo

    !   First merge very thin layers with the one above (or below if they are
    ! at the top), and find the interface quantities that the equation of state
    ! needs.
    do concurrent (jj=1:jje, ii=1:iie, G%mask2dT(isb+ii-1,jsb+jj-1) > 0.0) &
      DO_LOCALITY(local(i, j, k, hmin, H_here, HxT_here, HxS_here, HxR_here)) &
      DO_LOCALITY(reduce(max:max_kf))
      i = isb + ii - 1
      j = jsb + jj - 1
      htot(ii,jj) = 0.0
      do k=1,nz ; htot(ii,jj) = htot(ii,jj) + h(i,j,k) ; enddo

      hmin = htot(ii,jj)*min_h_frac ; kf(ii,jj) = 1 ; H_here = 0.0
      HxT_here = 0.0 ; HxS_here = 0.0 ; HxR_here = 0.0
      if (use_EOS) then
        do k=1,nz
          if ((H_here > hmin) .and. (h(i,j,k) > hmin)) then
            Hf(ii,jj,kf(ii,jj)) = H_here
            Tf(ii,jj,kf(ii,jj)) = HxT_here / H_here
            Sf(ii,jj,kf(ii,jj)) = HxS_here / H_here
            kf(ii,jj) = kf(ii,jj) + 1

            ! Start a new layer
            H_here = h(i,j,k)
            HxT_here = h(i,j,k) * tv%T(i,j,k)
            HxS_here = h(i,j,k) * tv%S(i,j,k)
          else
            H_here = H_here + h(i,j,k)
            HxT_here = HxT_here + h(i,j,k) * tv%T(i,j,k)
            HxS_here = HxS_here + h(i,j,k) * tv%S(i,j,k)
          endif
        enddo
        if (H_here > 0.0) then
          Hf(ii,jj,kf(ii,jj)) = H_here
          Tf(ii,jj,kf(ii,jj)) = HxT_here / H_here
          Sf(ii,jj,kf(ii,jj)) = HxS_here / H_here
        endif
      else  ! .not. (use_EOS)
        do k=1,nz
          if ((H_here > hmin) .and. (h(i,j,k) > hmin)) then
            Hf(ii,jj,kf(ii,jj)) = H_here ; Rf(ii,jj,kf(ii,jj)) = HxR_here / H_here
            kf(ii,jj) = kf(ii,jj) + 1

            ! Start a new layer
            H_here = h(i,j,k)
            HxR_here = h(i,j,k)*GV%Rlay(k)
          else
            H_here = H_here + h(i,j,k)
            HxR_here = HxR_here + h(i,j,k)*GV%Rlay(k)
          endif
        enddo
        if (H_here > 0.0) then
          Hf(ii,jj,kf(ii,jj)) = H_here ; Rf(ii,jj,kf(ii,jj)) = HxR_here / H_here
        endif
      endif

      !   The equation of state is evaluated below with a single call spanning a rectangular
      ! range of interfaces, but the number of filtered layers varies between columns.  Extend
      ! each column with zero-thickness copies of its deepest filtered layer so that the extra
      ! interfaces hold valid temperatures and salinities. These extra values are never used.
      if (use_EOS) then
        do k=kf(ii,jj)+1,nz+1
          Hf(ii,jj,k) = 0.0
          Tf(ii,jj,k) = Tf(ii,jj,kf(ii,jj)) ; Sf(ii,jj,k) = Sf(ii,jj,kf(ii,jj))
        enddo
      else
        do k=kf(ii,jj)+1,nz+1
          Hf(ii,jj,k) = 0.0 ; Rf(ii,jj,k) = Rf(ii,jj,kf(ii,jj))
        enddo
      endif
      max_kf = max(max_kf, kf(ii,jj))

      if (use_EOS) then
        pres(ii,jj,1) = 0.0
        T_int(ii,jj,1) = Tf(ii,jj,1) ; S_int(ii,jj,1) = Sf(ii,jj,1)
        do K=2,nz+1
          pres(ii,jj,K) = pres(ii,jj,K-1) + H_to_pres*Hf(ii,jj,k-1)
          T_int(ii,jj,K) = 0.5*(Tf(ii,jj,k)+Tf(ii,jj,k-1))
          S_int(ii,jj,K) = 0.5*(Sf(ii,jj,k)+Sf(ii,jj,k-1))
        enddo
      endif
    enddo

    if (use_EOS) then
      if (nonBous) then
        !$omp target update from(T_int, S_int, pres)
        do jj=1,jje ; do ii=1,iie
          call calculate_specific_vol_derivs(T_int(ii,jj,:), S_int(ii,jj,:), pres(ii,jj,:), &
                                             dSpV_dT(ii,jj,:), dSpV_dS(ii,jj,:), &
                                             tv%eqn_of_state, (/2,max_kf+1/) )
        enddo ; enddo
        !$omp target update to(dSpV_dT, dSpV_dS)
      else
        EOSdom(1,1) = 1 ; EOSdom(1,2) = iie
        EOSdom(2,1) = 1 ; EOSdom(2,2) = jje
        EOSdom(3,1) = 2 ; EOSdom(3,2) = max_kf+1
        call calculate_density_derivs(T_int, S_int, pres, drho_dT, drho_dS, &
                                      tv%eqn_of_state, EOSdom)
      endif
    endif

    !$omp target teams loop LOOP_BIND_TEAMS_PARALLEL collapse(2) private( i, j, k, k2, itt, kc, do_merge, &
    !$omp &                                      I_Htot, I_Hnew, det, ddet, lam, dlam, lam0, cg1_est, &
    !$omp &                                      speed2_tot, drxh_sum, dSpVxh_sum, hw, sum_hc, gp, &
    !$omp &                                      N2min, ms_min, ms_max, ms_sq, below_mono_N2_frac, &
    !$omp &                                      below_mono_N2_depth, Igu, Igl, H_top, H_bot, &
    !$omp &                                      I_beta, yy)
    do jj=1,jje ; do ii=1,iie
      i = isb + ii - 1
      j = jsb + jj - 1
      nkc(ii,jj) = 0
      if (G%mask2dT(i,j) > 0.0) then
        if (use_EOS) then

          ! Sum the reduced gravities to find out how small a density difference is negligibly small.
          drxh_sum = 0.0 ; dSpVxh_sum = 0.0
          if (better_est) then
            ! This is an estimate that is correct for the non-EBT mode for 2 or 3 layers, or for
            ! clusters of massless layers at interfaces that can be grouped into 2 or 3 layers.
            ! For a uniform stratification and a huge number of layers uniformly distributed in
            ! density, this estimate is too large (as is desired) by a factor of pi^2/6 ~= 1.64.
            H_top(1) = 0.0
            do K=2,kf(ii,jj) ; H_top(K) = H_top(K-1) + Hf(ii,jj,k-1) ; enddo
            if (H_top(kf(ii,jj)) > 0.0) then
              ! This is 1.0 / (H_top(K) + H_bot(K)) for all K.
              I_Htot = 1.0 / (H_top(kf(ii,jj)) + Hf(ii,jj,kf(ii,jj)))
              H_bot(kf(ii,jj)+1) = 0.0
              if (nonBous) then
                do K=kf(ii,jj),2,-1
                  H_bot(K) = H_bot(K+1) + Hf(ii,jj,k)
                  dSpVxh_sum = dSpVxh_sum + ((H_top(K) * H_bot(K)) * I_Htot) * &
                      min(0.0, dSpV_dT(ii,jj,K)*(Tf(ii,jj,k)-Tf(ii,jj,k-1)) + &
                          dSpV_dS(ii,jj,K)*(Sf(ii,jj,k)-Sf(ii,jj,k-1)))
                enddo
              else
                do K=kf(ii,jj),2,-1
                  H_bot(K) = H_bot(K+1) + Hf(ii,jj,k)
                  drxh_sum = drxh_sum + ((H_top(K) * H_bot(K)) * I_Htot) * &
                      max(0.0, drho_dT(ii,jj,K)*(Tf(ii,jj,k)-Tf(ii,jj,k-1)) + &
                          drho_dS(ii,jj,K)*(Sf(ii,jj,k)-Sf(ii,jj,k-1)))
                enddo
              endif
            endif
          else
            ! This estimate is problematic in that it goes like 1/nz for a large number of layers,
            ! but it is an overestimate (as desired) for a small number of layers, by at a factor
            ! of (H1+H2)**2/(H1*H2) >= 4 for two thick layers.
            if (nonBous) then
              do K=2,kf(ii,jj)
                dSpVxh_sum = dSpVxh_sum + 0.5*(Hf(ii,jj,k-1)+Hf(ii,jj,k)) * &
                    min(0.0, dSpV_dT(ii,jj,K)*(Tf(ii,jj,k)-Tf(ii,jj,k-1)) + &
                        dSpV_dS(ii,jj,K)*(Sf(ii,jj,k)-Sf(ii,jj,k-1)))
              enddo
            else
              do K=2,kf(ii,jj)
                drxh_sum = drxh_sum + 0.5*(Hf(ii,jj,k-1)+Hf(ii,jj,k)) * &
                    max(0.0, drho_dT(ii,jj,K)*(Tf(ii,jj,k)-Tf(ii,jj,k-1)) + &
                        drho_dS(ii,jj,K)*(Sf(ii,jj,k)-Sf(ii,jj,k-1)))
              enddo
            endif
          endif
        else  ! .not. (use_EOS)
          drxh_sum = 0.0 ; dSpVxh_sum = 0.0
          if (better_est) then
            H_top(1) = 0.0
            do K=2,kf(ii,jj) ; H_top(K) = H_top(K-1) + Hf(ii,jj,k-1) ; enddo
            if (H_top(kf(ii,jj)) > 0.0) then
              ! This is 1.0 / (H_top(K) + H_bot(K)) for all K.
              I_Htot = 1.0 / (H_top(kf(ii,jj)) + Hf(ii,jj,kf(ii,jj)))
              H_bot(kf(ii,jj)+1) = 0.0
              if (nonBous) then
                do K=kf(ii,jj),2,-1
                  H_bot(K) = H_bot(K+1) + Hf(ii,jj,k)
                  dSpVxh_sum = dSpVxh_sum + ((H_top(K) * H_bot(K)) * I_Htot) * &
                      min(0.0, (Rf(ii,jj,k-1)-Rf(ii,jj,k)) / (Rf(ii,jj,k)*Rf(ii,jj,k-1)))
                enddo
              else
                do K=kf(ii,jj),2,-1
                  H_bot(K) = H_bot(K+1) + Hf(ii,jj,k)
                  drxh_sum = drxh_sum + ((H_top(K) * H_bot(K)) * I_Htot) * &
                             max(0.0,Rf(ii,jj,k)-Rf(ii,jj,k-1))
                enddo
              endif
            endif
          else
            if (nonBous) then
              do K=2,kf(ii,jj)
                dSpVxh_sum = dSpVxh_sum + 0.5*(Hf(ii,jj,k-1)+Hf(ii,jj,k)) * &
                      min(0.0, (Rf(ii,jj,k-1)-Rf(ii,jj,k)) / (Rf(ii,jj,k)*Rf(ii,jj,k-1)))
              enddo
            else
              do K=2,kf(ii,jj)
                drxh_sum = drxh_sum + 0.5*(Hf(ii,jj,k-1)+Hf(ii,jj,k)) * max(0.0,Rf(ii,jj,k)-Rf(ii,jj,k-1))
              enddo
            endif
          endif
        endif ! use_EOS

        if (nonBous) then
          ! Note that dSpVxh_sum is negative for stable stratification.
          cg1_est = H_to_pres * abs(dSpVxh_sum)
        else
          cg1_est = g_Rho0 * drxh_sum
        endif

        !   Find gprime across each internal interface, taking care of convective instabilities by
        ! merging layers.  If the estimated wave speed is too small, simply return zero.
        if (cg1_est <= cg1_min2) then
          cg1(i,j) = 0.0
          if (present(modal_structure)) then
            do k=1,nz ; modal_structure(i,j,k) = 0. ; enddo
          endif
        else
          ! Merge layers to eliminate convective instabilities or exceedingly
          ! small reduced gravities.  Merging layers reduces the estimated wave speed by
          ! (rho(2)-rho(1))*h(1)*h(2) / H_tot.
          if (use_EOS) then
            kc = 1
            Hc(ii,jj,1) = Hf(ii,jj,1) ; Tc(ii,jj,1) = Tf(ii,jj,1) ; Sc(ii,jj,1) = Sf(ii,jj,1)
            do k=2,kf(ii,jj)
              if (better_est .and. nonBous) then
                do_merge = ((dSpV_dT(ii,jj,K)*(Tc(ii,jj,kc)-Tf(ii,jj,k)) + &
                             dSpV_dS(ii,jj,K)*(Sc(ii,jj,kc)-Sf(ii,jj,k))) * &
                             ((Hc(ii,jj,kc) * Hf(ii,jj,k))*I_Htot) < abs(2.0 * tol_merge * dSpVxh_sum))
              elseif (better_est) then
                do_merge = ((drho_dT(ii,jj,K)*(Tf(ii,jj,k)-Tc(ii,jj,kc)) + &
                             drho_dS(ii,jj,K)*(Sf(ii,jj,k)-Sc(ii,jj,kc))) * &
                             ((Hc(ii,jj,kc) * Hf(ii,jj,k))*I_Htot) < 2.0 * tol_merge*drxh_sum)
              elseif (nonBous) then
                do_merge = ((dSpV_dT(ii,jj,K)*(Tc(ii,jj,kc)-Tf(ii,jj,k)) + &
                             dSpV_dS(ii,jj,K)*(Sc(ii,jj,kc)-Sf(ii,jj,k))) * &
                             (Hc(ii,jj,kc) + Hf(ii,jj,k)) < abs(2.0 * tol_merge * dSpVxh_sum))
              else
                do_merge = ((drho_dT(ii,jj,K)*(Tf(ii,jj,k)-Tc(ii,jj,kc)) + &
                             drho_dS(ii,jj,K)*(Sf(ii,jj,k)-Sc(ii,jj,kc))) * &
                             (Hc(ii,jj,kc) + Hf(ii,jj,k)) < 2.0 * tol_merge*drxh_sum)
              endif
              if (do_merge) then
                ! Merge this layer with the one above and backtrack.
                I_Hnew = 1.0 / (Hc(ii,jj,kc) + Hf(ii,jj,k))
                Tc(ii,jj,kc) = (Hc(ii,jj,kc)*Tc(ii,jj,kc) + Hf(ii,jj,k)*Tf(ii,jj,k)) * I_Hnew
                Sc(ii,jj,kc) = (Hc(ii,jj,kc)*Sc(ii,jj,kc) + Hf(ii,jj,k)*Sf(ii,jj,k)) * I_Hnew
                Hc(ii,jj,kc) = (Hc(ii,jj,kc) + Hf(ii,jj,k))
                ! Backtrack to remove any convective instabilities above...  Note
                ! that the tolerance is a factor of two larger, to avoid limit how
                ! far back we go.
                do K2=kc,2,-1
                  if (better_est .and. nonBous) then
                    do_merge = ( (dSpV_dT(ii,jj,K2)*(Tc(ii,jj,k2-1)-Tc(ii,jj,k2)) + &
                               dSpV_dS(ii,jj,K2)*(Sc(ii,jj,k2-1)-Sc(ii,jj,k2))) * &
                              ((Hc(ii,jj,k2) * Hc(ii,jj,k2-1))*I_Htot) < abs(tol_merge * dSpVxh_sum) )
                  elseif (better_est) then
                    do_merge = ((drho_dT(ii,jj,K2)*(Tc(ii,jj,k2)-Tc(ii,jj,k2-1)) + &
                              drho_dS(ii,jj,K2)*(Sc(ii,jj,k2)-Sc(ii,jj,k2-1))) * &
                             ((Hc(ii,jj,k2) * Hc(ii,jj,k2-1))*I_Htot) < tol_merge*drxh_sum)
                  elseif (nonBous) then
                    do_merge = ( (dSpV_dT(ii,jj,K2)*(Tc(ii,jj,k2-1)-Tc(ii,jj,k2)) + &
                               dSpV_dS(ii,jj,K2)*(Sc(ii,jj,k2-1)-Sc(ii,jj,k2))) * &
                              (Hc(ii,jj,k2) + Hc(ii,jj,k2-1)) < abs(tol_merge * dSpVxh_sum) )
                  else
                    do_merge = ((drho_dT(ii,jj,K2)*(Tc(ii,jj,k2)-Tc(ii,jj,k2-1)) + &
                              drho_dS(ii,jj,K2)*(Sc(ii,jj,k2)-Sc(ii,jj,k2-1))) * &
                             (Hc(ii,jj,k2) + Hc(ii,jj,k2-1)) < tol_merge*drxh_sum)
                  endif
                  if (do_merge) then
                    ! Merge the two bottommost layers.  At this point kc = k2.
                    I_Hnew = 1.0 / (Hc(ii,jj,kc) + Hc(ii,jj,kc-1))
                    Tc(ii,jj,kc-1) = (Hc(ii,jj,kc)*Tc(ii,jj,kc) + Hc(ii,jj,kc-1)*Tc(ii,jj,kc-1)) * I_Hnew
                    Sc(ii,jj,kc-1) = (Hc(ii,jj,kc)*Sc(ii,jj,kc) + Hc(ii,jj,kc-1)*Sc(ii,jj,kc-1)) * I_Hnew
                    Hc(ii,jj,kc-1) = (Hc(ii,jj,kc) + Hc(ii,jj,kc-1))
                    kc = kc - 1
                  else ; exit ; endif
                enddo
              else
                ! Add a new layer to the column.
                kc = kc + 1
                if (nonBous) then
                  dSpV_dS(ii,jj,Kc) = dSpV_dS(ii,jj,K) ; dSpV_dT(ii,jj,Kc) = dSpV_dT(ii,jj,K)
                else
                  drho_dS(ii,jj,Kc) = drho_dS(ii,jj,K) ; drho_dT(ii,jj,Kc) = drho_dT(ii,jj,K)
                endif
                Tc(ii,jj,kc) = Tf(ii,jj,k) ; Sc(ii,jj,kc) = Sf(ii,jj,k) ; Hc(ii,jj,kc) = Hf(ii,jj,k)
              endif
            enddo
            ! At this point there are kc layers and the gprimes should be positive.
            if (nonBous) then
              do K=2,kc
                gprime(ii,jj,K) = H_to_pres * (dSpV_dT(ii,jj,K)*(Tc(ii,jj,k-1)-Tc(ii,jj,k)) + &
                                               dSpV_dS(ii,jj,K)*(Sc(ii,jj,k-1)-Sc(ii,jj,k)))
              enddo
            else
              do K=2,kc
                gprime(ii,jj,K) = g_Rho0 * (drho_dT(ii,jj,K)*(Tc(ii,jj,k)-Tc(ii,jj,k-1)) + &
                                            drho_dS(ii,jj,K)*(Sc(ii,jj,k)-Sc(ii,jj,k-1)))
              enddo
            endif
          else  ! .not. (use_EOS)
            ! Do the same with density directly...
            kc = 1
            Hc(ii,jj,1) = Hf(ii,jj,1) ; Rc(ii,jj,1) = Rf(ii,jj,1)
            do k=2,kf(ii,jj)
              if (nonBous .and. better_est) then
                do_merge = ((Rf(ii,jj,k) - Rc(ii,jj,kc)) *  ((Hc(ii,jj,kc) * Hf(ii,jj,k))*I_Htot) < &
                         (Rc(ii,jj,kc)*Rf(ii,jj,k)) * abs(2.0 * tol_merge * dSpVxh_sum))
              elseif (nonBous) then
                do_merge = ((Rf(ii,jj,k) - Rc(ii,jj,kc)) * (Hc(ii,jj,kc) + Hf(ii,jj,k)) < &
                         (Rc(ii,jj,kc)*Rf(ii,jj,k)) * abs(2.0 * tol_merge * dSpVxh_sum))
              elseif (better_est) then
                do_merge = ((Rf(ii,jj,k) - Rc(ii,jj,kc)) * &
                            ((Hc(ii,jj,kc) * Hf(ii,jj,k))*I_Htot) < 2.0*tol_merge*drxh_sum)
              else
                do_merge = ((Rf(ii,jj,k) - Rc(ii,jj,kc)) * (Hc(ii,jj,kc) + Hf(ii,jj,k)) < 2.0*tol_merge*drxh_sum)
              endif
              if (do_merge) then
                ! Merge this layer with the one above and backtrack.
                Rc(ii,jj,kc) = (Hc(ii,jj,kc)*Rc(ii,jj,kc) + Hf(ii,jj,k)*Rf(ii,jj,k)) / (Hc(ii,jj,kc) + Hf(ii,jj,k))
                Hc(ii,jj,kc) = (Hc(ii,jj,kc) + Hf(ii,jj,k))
                ! Backtrack to remove any convective instabilities above...  Note
                ! that the tolerance is a factor of two larger, to avoid limit how
                ! far back we go.
                do k2=kc,2,-1
                  if (nonBous .and. better_est) then
                    do_merge = ((Rc(ii,jj,k2) - Rc(ii,jj,k2-1)) *  ((Hc(ii,jj,kc) * Hf(ii,jj,k))*I_Htot) < &
                            (Rc(ii,jj,k2-1)*Rc(ii,jj,k2)) * abs(2.0 * tol_merge * dSpVxh_sum))
                  elseif (nonBous) then
                    do_merge = ((Rc(ii,jj,k2) - Rc(ii,jj,k2-1)) * (Hc(ii,jj,kc) + Hf(ii,jj,k)) < &
                             (Rc(ii,jj,k2-1)*Rc(ii,jj,k2)) * abs(2.0 * tol_merge * dSpVxh_sum))
                  elseif (better_est) then
                    do_merge = ((Rc(ii,jj,k2)-Rc(ii,jj,k2-1)) * ((Hc(ii,jj,k2) * Hc(ii,jj,k2-1))*I_Htot) < &
                             tol_merge*drxh_sum)
                  else
                    do_merge = ((Rc(ii,jj,k2)-Rc(ii,jj,k2-1)) * (Hc(ii,jj,k2)+Hc(ii,jj,k2-1)) < tol_merge*drxh_sum)
                  endif
                  if (do_merge) then
                    ! Merge the two bottommost layers.  At this point kc = k2.
                    Rc(ii,jj,kc-1) = (Hc(ii,jj,kc)*Rc(ii,jj,kc) + Hc(ii,jj,kc-1)*Rc(ii,jj,kc-1)) / (Hc(ii,jj,kc) + &
                                                                                                    Hc(ii,jj,kc-1))
                    Hc(ii,jj,kc-1) = (Hc(ii,jj,kc) + Hc(ii,jj,kc-1))
                    kc = kc - 1
                  else ; exit ; endif
                enddo
              else
                ! Add a new layer to the column.
                kc = kc + 1
                Rc(ii,jj,kc) = Rf(ii,jj,k) ; Hc(ii,jj,kc) = Hf(ii,jj,k)
              endif
            enddo
            ! At this point there are kc layers and the gprimes should be positive.
            if (nonBous) then
              do K=2,kc
                gprime(ii,jj,K) = H_to_pres * (Rc(ii,jj,k) - Rc(ii,jj,k-1)) / (Rc(ii,jj,k) * Rc(ii,jj,k-1))
              enddo
            else
              do K=2,kc
                gprime(ii,jj,K) = g_Rho0 * (Rc(ii,jj,k)-Rc(ii,jj,k-1))
              enddo
            endif
          endif  ! use_EOS

          ! Sum the contributions from all of the interfaces to give an over-estimate
          ! of the first-mode wave speed.  Also populate Igl and Igu which are the
          ! non-leading diagonals of the tridiagonal matrix.
          if (kc >= 2) then
            speed2_tot = 0.0
            if (better_est) then
              H_top(1) = 0.0 ; H_bot(kc+1) = 0.0
              do K=2,kc+1 ; H_top(K) = H_top(K-1) + Hc(ii,jj,k-1) ; enddo
              do K=kc,2,-1 ; H_bot(K) = H_bot(K+1) + Hc(ii,jj,k) ; enddo
              I_Htot = 0.0 ; if (H_top(kc+1) > 0.0) I_Htot = 1.0 / H_top(kc+1)
            endif

            if (l_use_ebt_mode) then
              Igu(1) = 0. ! Neumann condition for pressure modes
              sum_hc = Hc(ii,jj,1)
              N2min = gprime(ii,jj,2)/Hc(ii,jj,1)

              below_mono_N2_frac = .false.
              below_mono_N2_depth = .false.
              do k=2,kc
                hw = 0.5*(Hc(ii,jj,k-1)+Hc(ii,jj,k))
                gp = gprime(ii,jj,K)

                if (l_mono_N2_column_fraction>0. .or. l_mono_N2_depth>=0.) then
                  ! Determine whether N2 estimates should not be allowed to increase with depth.
                  if (l_mono_N2_column_fraction>0.) then
                    if (GV%Boussinesq .or. GV%semi_Boussinesq) then
                      below_mono_N2_frac = &
                          (max(G%meanSL(i,j) + G%bathyT(i,j), 0.0) - GV%H_to_Z * sum_hc < &
                           l_mono_N2_column_fraction * max(G%meanSL(i,j) + G%bathyT(i,j), 0.0))
                    else
                      below_mono_N2_frac = (htot(ii,jj) - sum_hc < l_mono_N2_column_fraction*htot(ii,jj))
                    endif
                  endif
                  if (l_mono_N2_depth >= 0.) below_mono_N2_depth = (sum_hc > l_mono_N2_depth)

                  if ( (gp > N2min*hw) .and. (below_mono_N2_frac .or. below_mono_N2_depth) ) then
                    ! Filters out regions where N2 increases with depth, but only in a lower fraction
                    ! of the water column or below a certain depth.
                    gp = N2min * hw
                  else
                    N2min = gp / hw
                  endif
                endif

                Igu(k) = 1.0/(gp*Hc(ii,jj,k))
                Igl(k-1) = 1.0/(gp*Hc(ii,jj,k-1))
                sum_hc = sum_hc + Hc(ii,jj,k)

                if (better_est) then
                  ! Estimate that the ebt_mode is sqrt(2) times the speed of the flat bottom modes.
                  speed2_tot = speed2_tot + 2.0 * gprime(ii,jj,K)*((H_top(K) * H_bot(K)) * I_Htot)
                else ! The ebt_mode wave should be faster than the flat-bottom mode, so 0.707 should be > 1?
                  speed2_tot = speed2_tot + gprime(ii,jj,K)*(Hc(ii,jj,k-1)+Hc(ii,jj,k))*0.707
                endif
              enddo
             !Igl(kc) = 0. ! Neumann condition for pressure modes
              Igl(kc) = 2.*Igu(kc) ! Dirichlet condition for pressure modes
            else ! .not. l_use_ebt_mode
              do K=2,kc
                Igl(K) = 1.0/(gprime(ii,jj,K)*Hc(ii,jj,k)) ; Igu(K) = 1.0/(gprime(ii,jj,K)*Hc(ii,jj,k-1))
                if (better_est) then
                  speed2_tot = speed2_tot + gprime(ii,jj,K)*((H_top(K) * H_bot(K)) * I_Htot)
                else
                  speed2_tot = speed2_tot + gprime(ii,jj,K)*(Hc(ii,jj,k-1)+Hc(ii,jj,k))
                endif
              enddo
            endif

            if (calc_modal_structure) then
              do k=1,nz ; mode_struct(ii,jj,k) = 0. ; enddo
              do k=1,kc ; mode_struct(ii,jj,k) = 1. ; enddo ! Uniform flow, first guess
            endif

            ! Under estimate the first eigenvalue (overestimate the speed) to start with.
            if (calc_modal_structure) then
              lam0 = 0.5 / speed2_tot ; lam = lam0
            else
              lam0 = 1.0 / speed2_tot ; lam = lam0
            endif
            ! Find the determinant and its derivative with lam.
            do itt=1,max_itt
              if (l_use_ebt_mode) then
                ! This initialization of det,ddet imply Neumann boundary conditions for horizontal
                ! velocity or pressure modes, so that first 3 rows of the matrix are
                !    /   b(1)-lam  igl(1)      0        0     0  ...  \
                !    |  igu(2)    b(2)-lam   igl(2)     0     0  ...  |
                !    |    0        igu(3)   b(3)-lam  igl(3)  0  ...  |
                ! The last two rows of the pressure equation matrix are
                !    |    ...  0  igu(kc-1)  b(kc-1)-lam  igl(kc-1)  |
                !    \    ...  0     0        igu(kc)     b(kc)-lam  /
                call tridiag_det(Igu, Igl, 1, kc, lam, det, ddet, row_scale=c2_scale)
              else
                ! This initialization of det,ddet imply Dirichlet boundary conditions for vertical
                ! velocity modes, so that first 3 rows of the matrix are
                !    /  b(2)-lam  igl(2)      0       0     0  ...  |
                !    |  igu(3)  b(3)-lam   igl(3)     0     0  ...  |
                !    |    0       igu(4)  b(4)-lam  igl(4)  0  ...  |
                ! The last three rows of the w equation matrix are
                !    |    ...   0  igu(kc-2)  b(kc-2)-lam  igl(kc-2)     0       |
                !    |    ...   0     0        igu(kc-1)  b(kc-1)-lam  igl(kc-1) |
                !    \    ...   0     0           0        igu(kc)    b(kc)-lam  /
                call tridiag_det(Igu, Igl, 2, kc, lam, det, ddet, row_scale=c2_scale)
              endif
              ! Use Newton's method iteration to find a new estimate of lam.

              if ((ddet >= 0.0) .or. (-det > -0.5*lam*ddet)) then
                ! lam was not an under-estimate, as intended, so Newton's method
                ! may not be reliable; lam must be reduced, but not by more
                ! than half.
                lam = 0.5 * lam
                dlam = -lam
              else  ! Newton's method is OK.
                dlam = - det / ddet
                lam = lam + dlam
              endif

              if (calc_modal_structure) then
                call tdma6_block(kc, Igu, Igl, nii, njj, nz, ii, jj, lam, mode_struct, I_beta, yy)
                ! Note that tdma6_block changes the units of mode_struct to [L2 T-2 ~> m2 s-2]
                ms_min = mode_struct(ii,jj,1)
                ms_max = mode_struct(ii,jj,1)
                ms_sq = mode_struct(ii,jj,1)**2
                do k = 2,kc
                  ms_min = min(ms_min, mode_struct(ii,jj,k))
                  ms_max = max(ms_max, mode_struct(ii,jj,k))
                  ms_sq = ms_sq + mode_struct(ii,jj,k)**2
                enddo
                if (ms_min<0. .and. ms_max>0.) then ! Any zero crossings => lam is too high
                  lam = 0.5 * ( lam - dlam )
                  dlam = -lam
                  do k=1,kc ; mode_struct(ii,jj,k) = abs(mode_struct(ii,jj,k)) / sqrt( ms_sq ) ; enddo
                else
                  do k=1,kc ; mode_struct(ii,jj,k) = mode_struct(ii,jj,k) / sqrt( ms_sq ) ; enddo
                endif
                ! After the nondimensionalization above, mode_struct is once again [nondim]
              endif

              if (abs(dlam) < tol_solve*lam) exit
            enddo

            cg1(i,j) = 0.0
            if (lam > 0.0) cg1(i,j) = 1.0 / sqrt(lam)

            if (present(modal_structure)) then
              if (mode_struct(ii,jj,1)/=0.) then ! Normalize
                do k=kc,1,-1 ; mode_struct(ii,jj,k) = mode_struct(ii,jj,k) / mode_struct(ii,jj,1) ; enddo
              else
                do k=1,kc ; mode_struct(ii,jj,k) = 0. ; enddo
              endif
            endif

            !   Record that this column has a mode structure that still needs to be remapped
            ! onto the model grid.  That remapping is done on the host, after this loop.
            nkc(ii,jj) = kc
          else
            cg1(i,j) = 0.0
            if (present(modal_structure)) then
              do k=1,nz ; modal_structure(i,j,k) = 0. ; enddo
            endif
          endif
        endif ! cg1 /= 0.0
      else
        cg1(i,j) = 0.0 ! This is a land point.
        if (present(modal_structure)) then
          do k=1,nz ; modal_structure(i,j,k) = 0. ; enddo
        endif
      endif
    enddo ; enddo ! column loop

    !   Remapping done on the host to avoid accessing polymorphic Recon1d type on device
    if (present(modal_structure)) then
      !$omp target update from(nkc, Hc, mode_struct)
      do jj=1,jje ; do ii=1,iie ; if (nkc(ii,jj) >= 2) then
        i = isb + ii - 1
        j = jsb + jj - 1
        if (CS%remap_answer_date < 20190101) then
          call remapping_core_h(CS%remap_2018_CS, nkc(ii,jj), Hc(ii,jj,:), mode_struct(ii,jj,:), &
                                nz, h(i,j,:), modal_structure(i,j,:))
        else
          call remapping_core_h(CS%remap_CS, nkc(ii,jj), Hc(ii,jj,:), mode_struct(ii,jj,:), &
                                nz, h(i,j,:), modal_structure(i,j,:))
        endif
      endif ; enddo ; enddo
      !$omp target update to(modal_structure)
    endif

  enddo ; enddo ! end of the i- and j-block loops

  !$omp target exit data map(delete: Hf, Tf, Sf, Rf, pres, T_int, S_int, dRho_dT, dRho_dS)
  !$omp target exit data map(delete: dSpV_dT, dSpV_dS, H_top, H_bot, gprime, Igl, Igu)
  !$omp target exit data map(delete: Hc, Tc, Sc, Rc, mode_struct)
  !$omp target exit data map(delete: htot, kf, nkc)

end subroutine wave_speed

!> Solve a non-symmetric tridiagonal problem with the sum of the upper and lower diagonals minus a
!! scalar contribution as the leading diagonal.
!! This uses the Thomas algorithm rather than the Hallberg algorithm since the matrix is not symmetric.
pure subroutine tdma6(n, a, c, lam, y)
  integer,            intent(in)    :: n !< Number of rows of matrix
  real, dimension(:), intent(in)    :: a !< Lower diagonal   [T2 L-2 ~> s2 m-2]
  real, dimension(:), intent(in)    :: c !< Upper diagonal   [T2 L-2 ~> s2 m-2]
  real,               intent(in)    :: lam !< Scalar subtracted from leading diagonal [T2 L-2 ~> s2 m-2]
  real, dimension(:), intent(inout) :: y !< RHS on entry [A ~> a], result on exit [A L2 T-2 ~> a m2 s-2]

  ! Local variables
  real :: lambda     ! A temporary variable in [T2 L-2 ~> s2 m-2]
  real :: beta(n)    ! A temporary variable in [T2 L-2 ~> s2 m-2]
  real :: I_beta(n)  ! A temporary variable in [L2 T-2 ~> m2 s-2]
  real :: yy(n)      ! A temporary variable with the same units as y on entry [A ~> a]
  integer :: k, m

  lambda = lam
  beta(1) = (a(1)+c(1)) - lambda
  if (beta(1)==0.) then ! lam was chosen too perfectly
    ! Change lambda and redo this first row
    lambda = (1. + 1.e-5) * lambda
    beta(1) = (a(1)+c(1)) - lambda
  endif
  I_beta(1) = 1. / beta(1)
  yy(1) = y(1)
  do k = 2, n
    beta(k) = ( (a(k)+c(k)) - lambda ) - a(k) * c(k-1) * I_beta(k-1)
    ! Perhaps the following 0 needs to become a tolerance to handle underflow?
    if (beta(k)==0.) then ! lam was chosen too perfectly
      ! Change lambda and redo everything up to row k
      lambda = (1. + 1.e-5) * lambda
      I_beta(1) = 1. / ( (a(1)+c(1)) - lambda )
      do m = 2, k
        I_beta(m) = 1. / ( ( (a(m)+c(m)) - lambda ) - a(m) * c(m-1) * I_beta(m-1) )
        yy(m) = y(m) + a(m) * yy(m-1) * I_beta(m-1)
      enddo
    else
      I_beta(k) = 1. / beta(k)
    endif
    yy(k) = y(k) + a(k) * yy(k-1) * I_beta(k-1)
  enddo
  ! The units of y change by a factor of [L2 T-2 ~> m2 s-2] in the following lines.
  y(n) = yy(n) * I_beta(n)
  do k = n-1, 1, -1
    y(k) = ( yy(k) + c(k) * y(k+1) ) * I_beta(k)
  enddo

end subroutine tdma6

! TODO: This routine can probably be merged with tdma6 once wave_speeds is also blocked
!> Solve the same non-symmetric tridiagonal problem as tdma6 for a single column of a
!! 3-dimensional array of right hand sides, addressing that column by its i- and j-indices.
!! The matrix diagonals and the elimination workspace are column-local arrays.
!!
pure subroutine tdma6_block(n, a, c, nii, njj, nk, i, j, lam, y, I_beta, yy)
  integer, intent(in)    :: n   !< Number of rows of matrix
  integer, intent(in)    :: nii !< Upper i-bound of the 3-d right hand side array
  integer, intent(in)    :: njj !< Upper j-bound of the 3-d right hand side array
  integer, intent(in)    :: nk  !< Vertical extent of the arrays
  real, dimension(nk), intent(in) :: a !< Lower diagonal   [T2 L-2 ~> s2 m-2]
  real, dimension(nk), intent(in) :: c !< Upper diagonal   [T2 L-2 ~> s2 m-2]
  integer, intent(in)    :: i   !< The i-index of the column to work on
  integer, intent(in)    :: j   !< The j-index of the column to work on
  real,    intent(in)    :: lam !< Scalar subtracted from leading diagonal [T2 L-2 ~> s2 m-2]
  real, dimension(1:nii,1:njj,nk), intent(inout) :: y !< RHS on entry [A ~> a], result on exit
                                                     !! [A L2 T-2 ~> a m2 s-2]
  real, dimension(nk), intent(inout) :: I_beta !< Workspace holding the inverse of the leading
                                               !! diagonal after elimination [L2 T-2 ~> m2 s-2]
  real, dimension(nk), intent(inout) :: yy !< Workspace holding the forward-eliminated right hand
                                           !! side, with the same units as y on entry [A ~> a]
  !$omp declare target

  ! Local variables
  real :: lambda     ! A temporary variable in [T2 L-2 ~> s2 m-2]
  real :: beta       ! A temporary variable in [T2 L-2 ~> s2 m-2].  Unlike in tdma6 this is a scalar,
                     ! because it is only ever read back at the row where it was set.
  integer :: k, m

  lambda = lam
  beta = (a(1)+c(1)) - lambda
  if (beta==0.) then ! lam was chosen too perfectly
    ! Change lambda and redo this first row
    lambda = (1. + 1.e-5) * lambda
    beta = (a(1)+c(1)) - lambda
  endif
  I_beta(1) = 1. / beta
  yy(1) = y(i,j,1)
  do k = 2, n
    beta = ( (a(k)+c(k)) - lambda ) - a(k) * c(k-1) * I_beta(k-1)
    ! Perhaps the following 0 needs to become a tolerance to handle underflow?
    if (beta==0.) then ! lam was chosen too perfectly
      ! Change lambda and redo everything up to row k
      lambda = (1. + 1.e-5) * lambda
      I_beta(1) = 1. / ( (a(1)+c(1)) - lambda )
      do m = 2, k
        I_beta(m) = 1. / ( ( (a(m)+c(m)) - lambda ) - &
                           a(m) * c(m-1) * I_beta(m-1) )
        yy(m) = y(i,j,m) + a(m) * yy(m-1) * I_beta(m-1)
      enddo
    else
      I_beta(k) = 1. / beta
    endif
    yy(k) = y(i,j,k) + a(k) * yy(k-1) * I_beta(k-1)
  enddo
  ! The units of y change by a factor of [L2 T-2 ~> m2 s-2] in the following lines.
  y(i,j,n) = yy(n) * I_beta(n)
  do k = n-1, 1, -1
    y(i,j,k) = ( yy(k) + c(k) * y(i,j,k+1) ) * I_beta(k)
  enddo

end subroutine tdma6_block

!> Calculates the wave speeds for the first few barolinic modes.
subroutine wave_speeds(h, tv, G, GV, US, nmodes, cn, CS, w_struct, u_struct, u_struct_max, u_struct_bot, Nb, int_w2, &
                       int_U2, int_N2w2, halo_size)
  type(ocean_grid_type),                           intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),                         intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),                           intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),       intent(in)  :: h  !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                           intent(in)  :: tv !< Thermodynamic variables
  integer,                                         intent(in)  :: nmodes !< Number of modes
  type(wave_speed_CS),                             intent(in)  :: CS !< Wave speed control struct
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1,nmodes),intent(out) :: w_struct !< Wave vertical velocity profile [nondim]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV),nmodes),intent(out) :: u_struct !< Wave horizontal velocity profile
                                                                     !! [Z-1 ~> m-1]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: cn !< Waves speeds [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: u_struct_max !< Maximum of wave horizontal velocity
                                                                     !! profile [Z-1 ~> m-1]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: u_struct_bot !< Bottom value of wave horizontal
                                                                     !! velocity profile [Z-1 ~> m-1]
  real, dimension(SZI_(G),SZJ_(G)),                intent(out) :: Nb !< Bottom value of buoyancy freqency
                                                                     !! [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: int_w2 !< depth-integrated vertical velocity
                                                                     !! profile squared [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: int_U2 !< depth-integrated horizontal velocity
                                                                     !! profile squared [H Z-2 ~> m-1 or kg m-4]
  real, dimension(SZI_(G),SZJ_(G),nmodes),         intent(out) :: int_N2w2 !< depth-integrated buoyancy frequency
                                                                     !! times vertical velocity profile
                                                                     !! squared [H T-2 ~> m s-2 or kg m-2 s-2]
  integer,                               optional, intent(in)  :: halo_size !< Width of halo within which to
                                                                     !! calculate wave speeds

  ! Local variables
  real, dimension(SZK_(GV)+1) :: &
    dRho_dT, &    ! Partial derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
    dRho_dS, &    ! Partial derivative of density with salinity [R S-1 ~> kg m-3 ppt-1]
    dSpV_dT, &    ! Partial derivative of specific volume with temperature [R-1 C-1 ~> m3 kg-1 degC-1]
    dSpV_dS, &    ! Partial derivative of specific volume with salinity [R-1 S-1 ~> m3 kg-1 ppt-1]
    pres, &       ! Interface pressure [R L2 T-2 ~> Pa]
    T_int, &      ! Temperature interpolated to interfaces [C ~> degC]
    S_int, &      ! Salinity interpolated to interfaces [S ~> ppt]
    H_top, &      ! The distance of each filtered interface from the ocean surface [H ~> m or kg m-2]
    H_bot, &      ! The distance of each filtered interface from the bottom [H ~> m or kg m-2]
    gprime, &     ! The reduced gravity across each interface [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
    N2            ! The buoyancy freqency squared [T-2 ~> s-2]
  real, dimension(SZK_(GV),SZI_(G)) :: &
    Hf, &         ! Layer thicknesses after very thin layers are combined [H ~> m or kg m-2]
    dzf, &        ! Layer vertical extents after very thin layers are combined [Z ~> m]
    Tf, &         ! Layer temperatures after very thin layers are combined [C ~> degC]
    Sf, &         ! Layer salinities after very thin layers are combined [S ~> ppt]
    Rf            ! Layer densities after very thin layers are combined [R ~> kg m-3]
  real, dimension(SZI_(G),SZK_(GV)) :: &
    dz_2d         ! Height change across layers [Z ~> m]
  real, dimension(SZK_(GV)) :: &
    Igl, Igu, &   ! The inverse of the reduced gravity across an interface times
                  ! the thickness of the layer below (Igl) or above (Igu) it, in [T2 L-2 ~> s2 m-2].
    Hc, &         ! A column of layer thicknesses after convective instabilities are removed [H ~> m or kg m-2]
    dzc, &        ! A column of layer vertical extents after convective instabilities are removed [Z ~> m]
    Tc, &         ! A column of layer temperatures after convective instabilities are removed [C ~> degC]
    Sc, &         ! A column of layer salinities after convective instabilities are removed [S ~> ppt]
    Rc            ! A column of layer densities after convective instabilities are removed [R ~> kg m-3]
  real :: I_Htot  ! The inverse of the total filtered thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: c2_scale ! A scaling factor for wave speeds to help control the growth of the determinant and its
                   ! derivative with lam between rows of the Thomas algorithm solver [L2 s2 T-2 m-2 ~> nondim].
                   ! The exact value should not matter for the final result if it is an even power of 2.
  real :: det, ddet       ! Determinant of the eigen system and its derivative with lam.  Because the
                          ! units of the eigenvalue change with the number of layers and because of the
                          ! dynamic rescaling that is used to keep det in a numerically representable range,
                          ! the units of of det are hard to interpret, but det/ddet is always in units
                          ! of [T2 L-2 ~> s2 m-2]
  real :: lam_1           ! approximate mode-1 eigenvalue [T2 L-2 ~> s2 m-2]
  real :: lam_n           ! approximate mode-n eigenvalue [T2 L-2 ~> s2 m-2]
  real :: dlam            ! The change in estimates of the eigenvalue [T2 L-2 ~> s2 m-2]
  real :: lamMin          ! minimum lam value for root searching range [T2 L-2 ~> s2 m-2]
  real :: lamMax          ! maximum lam value for root searching range [T2 L-2 ~> s2 m-2]
  real :: lamInc          ! width of moving window for root searching [T2 L-2 ~> s2 m-2]
  real :: det_l, ddet_l   ! determinant of the eigensystem and its derivative with lam at the lower
                          ! end of the range of values bracketing a particular root, in dynamically
                          ! rescaled units that may differ from the other det variables, but such
                          ! that the units of det_l/ddet_l are [T2 L-2 ~> s2 m-2]
  real :: det_r, ddet_r   ! determinant and its derivative with lam at the lower end of the
                          ! bracket in arbitrarily rescaled units, but such that the units of
                          ! det_r/ddet_r are [T2 L-2 ~> s2 m-2]
  real :: det_sub, ddet_sub ! determinant and its derivative with lam at a subinterval endpoint that
                          ! is a candidate for a new bracket endpoint in arbitrarily rescaled units,
                          ! but such that the units of det_sub/ddet_sub are [T2 L-2 ~> s2 m-2]
  real :: xl, xr          ! lam guesses at left and right of window [T2 L-2 ~> s2 m-2]
  real :: xl_sub          ! lam guess at left of subinterval window [T2 L-2 ~> s2 m-2]
  real, dimension(nmodes) :: &
          xbl, xbr        ! lam guesses bracketing a zero-crossing (root) [T2 L-2 ~> s2 m-2]
  integer :: numint       ! number of widows (intervals) in root searching range
  integer :: nrootsfound  ! number of extra roots found (not including 1st root)
  real :: H_to_pres ! A conversion factor from thicknesses to pressure [R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1]
  real, dimension(SZI_(G)) :: &
    htot, hmin, &  ! Thicknesses [H ~> m or kg m-2]
    H_here, &      ! A layer thickness [H ~> m or kg m-2]
    dz_here, &     ! A layer vertical extent [Z ~> m]
    HxT_here, &    ! A layer integrated temperature [C H ~> degC m or degC kg m-2]
    HxS_here, &    ! A layer integrated salinity [S H ~> ppt m or ppt kg m-2]
    HxR_here       ! A layer integrated density [R H ~> kg m-2 or kg2 m-5]
  real :: speed2_tot ! overestimate of the mode-1 speed squared [L2 T-2 ~> m2 s-2]
  real :: speed2_min ! minimum mode speed (squared) to consider in root searching [L2 T-2 ~> m2 s-2]
  real :: cg1_min2   ! A floor in the squared first mode speed below which 0 is returned [L2 T-2 ~> m2 s-2]
  real :: cg1_est    ! An initial estimate of the squared first mode speed [L2 T-2 ~> m2 s-2]
  real, parameter :: reduct_factor = 0.5  ! A factor used in setting speed2_min [nondim]
  real :: I_Hnew     ! The inverse of a new layer thickness [H-1 ~> m-1 or m2 kg-1]
  real :: drxh_sum   ! The sum of density differences across interfaces times thicknesses [R H ~> kg m-2 or kg2 m-5]
  real :: dSpVxh_sum ! The sum of specific volume differences across interfaces times
                     ! thicknesses [H R-1 ~> m4 kg-1 or m], negative for stable stratification.
  real :: g_Rho0     ! G_Earth/Rho0 [L2 T-2 H-1 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2].
  real :: tol_Hfrac  ! Layers that together are smaller than this fraction of
                     ! the total water column can be merged for efficiency [nondim].
  real :: min_h_frac ! tol_Hfrac divided by the total number of layers [nondim].
  real :: tol_solve  ! The fractional tolerance with which to solve for the wave speeds [nondim].
  real :: tol_merge  ! The fractional change in estimated wave speed that is allowed
                     ! when deciding to merge layers in the calculation [nondim]
  integer :: kf(SZI_(G)) ! The number of active layers after filtering.
  integer, parameter :: max_itt = 30
  logical :: use_EOS    ! If true, density or specific volume is calculated from T & S using the equation of state.
  logical :: nonBous    ! If true, do not make the Boussinesq approximation.
  logical :: better_est ! If true, use an improved estimate of the first mode internal wave speed.
  logical :: merge      ! If true, merge the current layer with the one above.
  integer :: nsub       ! number of subintervals used for root finding
  integer, parameter :: sub_it_max = 4
                        ! maximum number of times to subdivide interval
                        ! for root finding (# intervals = 2**sub_it_max)
  logical :: sub_rootfound ! if true, subdivision has located root
  integer :: kc         ! The number of layers in the column after merging
  integer :: sub, sub_it
  integer :: i, j, k, k2, itt, is, ie, js, je, nz, iint, m, halo
  real, dimension(SZK_(GV)+1) ::  modal_structure !< Normalized model structure [nondim]
  real, dimension(SZK_(GV)) ::  modal_structure_fder !< Normalized model structure [Z-1 ~> m-1]
  real :: mode_struct(SZK_(GV)+1) ! The mode structure [nondim], but it is also temporarily
                         ! in units of [L2 T-2 ~> m2 s-2] after it is modified inside of tdma6.
  real :: mode_struct_fder(SZK_(GV)) ! The mode structure 1st derivative [Z-1 ~> m-1], but it is also temporarily
                         ! in units of [L2 Z-1 T-2 ~> m s-2] after it is modified inside of tdma6.
  real :: mode_struct_sq(SZK_(GV)+1) ! The square of mode structure [nondim]
  real :: mode_struct_fder_sq(SZK_(GV)) ! The square of mode structure 1st derivative [Z-2 ~> m-2]

  real :: w2avg          ! A total for renormalization [H L4 T-4 ~> m5 s-4 or kg m2 s-4]
  real, parameter :: a_int = 0.5 ! Integral total for normalization [nondim]
  real :: renorm         ! Normalization factor [T2 L-2 ~> s2 m-2]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke ; halo = 0

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_wave_speed / wave_speeds: "// &
         "Module must be initialized before it is used.")

  if (present(halo_size)) then
    halo = halo_size
    is = G%isc - halo ; ie = G%iec + halo ; js = G%jsc - halo ; je = G%jec + halo
  endif

  nonBous = .not.(GV%Boussinesq .or. GV%semi_Boussinesq)
  H_to_pres = GV%H_to_RZ * GV%g_Earth
  if (.not.nonBous) g_Rho0 = GV%g_Earth * GV%H_to_Z / GV%Rho0
  use_EOS = associated(tv%eqn_of_state)

  if (CS%c1_thresh < 0.0) &
    call MOM_error(FATAL, "INTERNAL_WAVE_CG1_THRESH must be set to a non-negative "//&
                          "value via wave_speed_init for wave_speeds to be used.")
  c2_scale = US%m_s_to_L_T**2 / 4096.0**2 ! Other powers of 2 give identical results.

  better_est = CS%better_cg1_est
  if (better_est) then
    tol_solve = CS%wave_speed_tol
    tol_Hfrac  = 0.1*tol_solve ; tol_merge = tol_solve / real(nz)
  else
    tol_solve = 0.001 ; tol_Hfrac  = 0.0001 ; tol_merge = 0.001
  endif
  cg1_min2 = CS%min_speed2

  ! Zero out all local values.  Values over land or for columns that are too weakly stratified
  ! are not changed from this zero value.
  cn(:,:,:) = 0.0
  u_struct_max(:,:,:) = 0.0
  u_struct_bot(:,:,:) = 0.0
  Nb(:,:) = 0.0
  int_w2(:,:,:) = 0.0
  int_N2w2(:,:,:) = 0.0
  int_U2(:,:,:) = 0.0
  u_struct(:,:,:,:) = 0.0
  w_struct(:,:,:,:) = 0.0

  min_h_frac = tol_Hfrac / real(nz)
  !$OMP parallel do default(private) shared(is,ie,js,je,nz,h,G,GV,US,CS,use_EOS,nonBous, &
  !$OMP                                     min_h_frac,H_to_pres,tv,cn,g_Rho0,nmodes,cg1_min2, &
  !$OMP                                     better_est,tol_solve,tol_merge,c2_scale)
  do j=js,je
    !   First merge very thin layers with the one above (or below if they are
    ! at the top).  This also transposes the row order so that columns can
    ! be worked upon one at a time.
    do i=is,ie ; htot(i) = 0.0 ; enddo
    do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo

    call thickness_to_dz(h, tv, dz_2d, j, G, GV, halo_size=halo)

    do i=is,ie
      hmin(i) = htot(i)*min_h_frac ; kf(i) = 1 ; H_here(i) = 0.0 ; dz_here(i) = 0.0
      HxT_here(i) = 0.0 ; HxS_here(i) = 0.0 ; HxR_here(i) = 0.0
    enddo
    if (use_EOS) then
      do k=1,nz ; do i=is,ie
        if ((H_here(i) > hmin(i)) .and. (h(i,j,k) > hmin(i))) then
          Hf(kf(i),i) = H_here(i)
          dzf(kf(i),i) = dz_here(i)
          Tf(kf(i),i) = HxT_here(i) / H_here(i)
          Sf(kf(i),i) = HxS_here(i) / H_here(i)
          kf(i) = kf(i) + 1

          ! Start a new layer
          H_here(i) = h(i,j,k)
          dz_here(i) = dz_2d(i,k)
          HxT_here(i) = h(i,j,k)*tv%T(i,j,k)
          HxS_here(i) = h(i,j,k)*tv%S(i,j,k)
        else
          H_here(i) = H_here(i) + h(i,j,k)
          dz_here(i) = dz_here(i) + dz_2d(i,k)
          HxT_here(i) = HxT_here(i) + h(i,j,k)*tv%T(i,j,k)
          HxS_here(i) = HxS_here(i) + h(i,j,k)*tv%S(i,j,k)
        endif
      enddo ; enddo
      do i=is,ie ; if (H_here(i) > 0.0) then
        Hf(kf(i),i) = H_here(i)
        dzf(kf(i),i) = dz_here(i)
        Tf(kf(i),i) = HxT_here(i) / H_here(i)
        Sf(kf(i),i) = HxS_here(i) / H_here(i)
      endif ; enddo
    else  ! .not. (use_EOS)
      do k=1,nz ; do i=is,ie
        if ((H_here(i) > hmin(i)) .and. (h(i,j,k) > hmin(i))) then
          Hf(kf(i),i) = H_here(i) ; Rf(kf(i),i) = HxR_here(i) / H_here(i)
          dzf(kf(i),i) = dz_here(i)
          kf(i) = kf(i) + 1

          ! Start a new layer
          H_here(i) = h(i,j,k)
          dz_here(i) = dz_2d(i,k)
          HxR_here(i) = h(i,j,k)*GV%Rlay(k)
        else
          H_here(i) = H_here(i) + h(i,j,k)
          dz_here(i) = dz_here(i) + dz_2d(i,k)
          HxR_here(i) = HxR_here(i) + h(i,j,k)*GV%Rlay(k)
        endif
      enddo ; enddo
      do i=is,ie ; if (H_here(i) > 0.0) then
        Hf(kf(i),i) = H_here(i) ; Rf(kf(i),i) = HxR_here(i) / H_here(i)
        dzf(kf(i),i) = dz_here(i)
      endif ; enddo
    endif

    ! From this point, we can work on individual columns without causing memory to have page faults.
    do i=is,ie
      if (G%mask2dT(i,j) > 0.0) then
        if (use_EOS) then
          pres(1) = 0.0 ; H_top(1) = 0.0
          do K=2,kf(i)
            pres(K) = pres(K-1) + H_to_pres*Hf(k-1,i)
            T_int(K) = 0.5*(Tf(k,i)+Tf(k-1,i))
            S_int(K) = 0.5*(Sf(k,i)+Sf(k-1,i))
            H_top(K) = H_top(K-1) + Hf(k-1,i)
          enddo
          if (nonBous) then
            call calculate_specific_vol_derivs(T_int, S_int, pres, dSpV_dT, dSpV_dS, &
                                               tv%eqn_of_state, (/2,kf(i)/) )
          else
            call calculate_density_derivs(T_int, S_int, pres, drho_dT, drho_dS, &
                                          tv%eqn_of_state, (/2,kf(i)/) )
          endif

          ! Sum the reduced gravities to find out how small a density difference is negligibly small.
          drxh_sum = 0.0 ; dSpVxh_sum = 0.0
          if (better_est) then
            ! This is an estimate that is correct for the non-EBT mode for 2 or 3 layers, or for
            ! clusters of massless layers at interfaces that can be grouped into 2 or 3 layers.
            ! For a uniform stratification and a huge number of layers uniformly distributed in
            ! density, this estimate is too large (as is desired) by a factor of pi^2/6 ~= 1.64.
            if (H_top(kf(i)) > 0.0) then
              I_Htot = 1.0 / (H_top(kf(i)) + Hf(kf(i),i))  ! = 1.0 / (H_top(K) + H_bot(K)) for all K.
              H_bot(kf(i)+1) = 0.0
              if (nonBous) then
                do K=kf(i),2,-1
                  H_bot(K) = H_bot(K+1) + Hf(k,i)
                  dSpVxh_sum = dSpVxh_sum + ((H_top(K) * H_bot(K)) * I_Htot) * &
                      min(0.0, dSpV_dT(K)*(Tf(k,i)-Tf(k-1,i)) + dSpV_dS(K)*(Sf(k,i)-Sf(k-1,i)))
                enddo
              else
                do K=kf(i),2,-1
                  H_bot(K) = H_bot(K+1) + Hf(k,i)
                  drxh_sum = drxh_sum + ((H_top(K) * H_bot(K)) * I_Htot) * &
                      max(0.0, drho_dT(K)*(Tf(k,i)-Tf(k-1,i)) + drho_dS(K)*(Sf(k,i)-Sf(k-1,i)))
                enddo
              endif
            endif
          else
            ! This estimate is problematic in that it goes like 1/nz for a large number of layers,
            ! but it is an overestimate (as desired) for a small number of layers, by at a factor
            ! of (H1+H2)**2/(H1*H2) >= 4 for two thick layers.
            if (nonBous) then
              do K=2,kf(i)
                dSpVxh_sum = dSpVxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * &
                    min(0.0, dSpV_dT(K)*(Tf(k,i)-Tf(k-1,i)) + dSpV_dS(K)*(Sf(k,i)-Sf(k-1,i)))
              enddo
            else
              do K=2,kf(i)
                drxh_sum = drxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * &
                    max(0.0, drho_dT(K)*(Tf(k,i)-Tf(k-1,i)) + drho_dS(K)*(Sf(k,i)-Sf(k-1,i)))
              enddo
            endif
          endif
        else  ! Not use_EOS
          drxh_sum = 0.0 ; dSpVxh_sum = 0.0
          if (better_est) then
            H_top(1) = 0.0
            do K=2,kf(i) ; H_top(K) = H_top(K-1) + Hf(k-1,i) ; enddo
              if (H_top(kf(i)) > 0.0) then
              I_Htot = 1.0 / (H_top(kf(i)) + Hf(kf(i),i))  ! = 1.0 / (H_top(K) + H_bot(K)) for all K.
              H_bot(kf(i)+1) = 0.0
              if (nonBous) then
                do K=kf(i),2,-1
                  H_bot(K) = H_bot(K+1) + Hf(k,i)
                  dSpVxh_sum = dSpVxh_sum + ((H_top(K) * H_bot(K)) * I_Htot) * &
                      min(0.0, (Rf(k-1,i)-Rf(k,i)) / (Rf(k,i)*Rf(k-1,i)))
                enddo
              else
                do K=kf(i),2,-1
                  H_bot(K) = H_bot(K+1) + Hf(k,i)
                  drxh_sum = drxh_sum + ((H_top(K) * H_bot(K)) * I_Htot) * max(0.0,Rf(k,i)-Rf(k-1,i))
                enddo
              endif
            endif
          else
            do K=2,kf(i)
              drxh_sum = drxh_sum + 0.5*(Hf(k-1,i)+Hf(k,i)) * max(0.0,Rf(k,i)-Rf(k-1,i))
            enddo
          endif
        endif

        if (nonBous) then
          ! Note that dSpVxh_sum is negative for stable stratification.
          cg1_est = H_to_pres * abs(dSpVxh_sum)
        else
          cg1_est = g_Rho0 * drxh_sum
        endif

        !   Find gprime across each internal interface, taking care of convective
        ! instabilities by merging layers.
        if (cg1_est > cg1_min2) then
          ! Merge layers to eliminate convective instabilities or exceedingly
          ! small reduced gravities.  Merging layers reduces the estimated wave speed by
          ! (rho(2)-rho(1))*h(1)*h(2) / H_tot.
          if (use_EOS) then
            kc = 1
            Hc(1) = Hf(1,i) ; dzc(1) = dzf(1,i) ; Tc(1) = Tf(1,i) ; Sc(1) = Sf(1,i)
            do k=2,kf(i)
              if (better_est .and. nonBous) then
                merge = ((dSpV_dT(K)*(Tc(kc)-Tf(k,i)) + dSpV_dS(K)*(Sc(kc)-Sf(k,i))) * &
                         ((Hc(kc) * Hf(k,i))*I_Htot) < abs(2.0 * tol_merge * dSpVxh_sum))
              elseif (better_est) then
                merge = ((drho_dT(K)*(Tf(k,i)-Tc(kc)) + drho_dS(K)*(Sf(k,i)-Sc(kc))) * &
                         ((Hc(kc) * Hf(k,i))*I_Htot) < 2.0 * tol_merge*drxh_sum)
              elseif (nonBous) then
                merge = ((dSpV_dT(K)*(Tc(kc)-Tf(k,i)) + dSpV_dS(K)*(Sc(kc)-Sf(k,i))) * &
                         (Hc(kc) + Hf(k,i)) < abs(2.0 * tol_merge * dSpVxh_sum))
              else
                merge = ((drho_dT(K)*(Tf(k,i)-Tc(kc)) + drho_dS(K)*(Sf(k,i)-Sc(kc))) * &
                         (Hc(kc) + Hf(k,i)) < 2.0 * tol_merge*drxh_sum)
              endif
              if (merge) then
                ! Merge this layer with the one above and backtrack.
                I_Hnew = 1.0 / (Hc(kc) + Hf(k,i))
                Tc(kc) = (Hc(kc)*Tc(kc) + Hf(k,i)*Tf(k,i)) * I_Hnew
                Sc(kc) = (Hc(kc)*Sc(kc) + Hf(k,i)*Sf(k,i)) * I_Hnew
                Hc(kc) = Hc(kc) + Hf(k,i)
                dzc(kc) = dzc(kc) + dzf(k,i)
                ! Backtrack to remove any convective instabilities above...  Note
                ! that the tolerance is a factor of two larger, to avoid limit how
                ! far back we go.
                do K2=kc,2,-1
                  if (better_est .and. nonBous) then
                    merge = ( (dSpV_dT(K2)*(Tc(k2-1)-Tc(k2)) + dSpV_dS(K2)*(Sc(k2-1)-Sc(k2))) * &
                              ((Hc(k2) * Hc(k2-1))*I_Htot) < abs(tol_merge * dSpVxh_sum) )
                  elseif (better_est) then
                    merge = ((drho_dT(K2)*(Tc(k2)-Tc(k2-1)) + drho_dS(K2)*(Sc(k2)-Sc(k2-1))) * &
                             ((Hc(k2) * Hc(k2-1))*I_Htot) < tol_merge*drxh_sum)
                  elseif (nonBous) then
                    merge = ( (dSpV_dT(K2)*(Tc(k2-1)-Tc(k2)) + dSpV_dS(K2)*(Sc(k2-1)-Sc(k2))) * &
                              (Hc(k2) + Hc(k2-1)) < abs(tol_merge * dSpVxh_sum) )
                  else
                    merge = ((drho_dT(K2)*(Tc(k2)-Tc(k2-1)) + drho_dS(K2)*(Sc(k2)-Sc(k2-1))) * &
                             (Hc(k2) + Hc(k2-1)) < tol_merge*drxh_sum)
                  endif
                  if (merge) then
                    ! Merge the two bottommost layers.  At this point kc = k2.
                    I_Hnew = 1.0 / (Hc(kc) + Hc(kc-1))
                    Tc(kc-1) = (Hc(kc)*Tc(kc) + Hc(kc-1)*Tc(kc-1)) * I_Hnew
                    Sc(kc-1) = (Hc(kc)*Sc(kc) + Hc(kc-1)*Sc(kc-1)) * I_Hnew
                    Hc(kc-1) = Hc(kc) + Hc(kc-1)
                    dzc(kc-1) = dzc(kc) + dzc(kc-1)
                    kc = kc - 1
                  else ; exit ; endif
                enddo
              else
                ! Add a new layer to the column.
                kc = kc + 1
                if (nonBous) then
                  dSpV_dS(Kc) = dSpV_dS(K) ; dSpV_dT(Kc) = dSpV_dT(K)
                else
                  drho_dS(Kc) = drho_dS(K) ; drho_dT(Kc) = drho_dT(K)
                endif
                Tc(kc) = Tf(k,i) ; Sc(kc) = Sf(k,i) ; Hc(kc) = Hf(k,i) ; dzc(kc) = dzf(k,i)
              endif
            enddo
            ! At this point there are kc layers and the gprimes should be positive.
            if (nonBous) then
              do K=2,kc
                gprime(K) = H_to_pres * (dSpV_dT(K)*(Tc(k-1)-Tc(k)) + dSpV_dS(K)*(Sc(k-1)-Sc(k)))
              enddo
            else
              do K=2,kc
                gprime(K) = g_Rho0 * (drho_dT(K)*(Tc(k)-Tc(k-1)) + drho_dS(K)*(Sc(k)-Sc(k-1)))
              enddo
            endif
          else  ! .not. (use_EOS)
            ! Do the same with density directly...
            kc = 1
            Hc(1) = Hf(1,i) ; dzc(1) = dzf(1,i) ; Rc(1) = Rf(1,i)
            do k=2,kf(i)
              if (nonBous .and. better_est) then
                merge = ((Rf(k,i) - Rc(kc)) *  ((Hc(kc) * Hf(k,i))*I_Htot) < &
                         (Rc(kc)*Rf(k,i)) * abs(2.0 * tol_merge * dSpVxh_sum))
              elseif (nonBous) then
                merge = ((Rf(k,i) - Rc(kc)) * (Hc(kc) + Hf(k,i)) < &
                         (Rc(kc)*Rf(k,i)) * abs(2.0 * tol_merge * dSpVxh_sum))
              elseif (better_est) then
                merge = ((Rf(k,i) - Rc(kc)) * ((Hc(kc) * Hf(k,i))*I_Htot) < 2.0*tol_merge*drxh_sum)
              else
                merge = ((Rf(k,i) - Rc(kc)) * (Hc(kc) + Hf(k,i)) < 2.0*tol_merge*drxh_sum)
              endif
              if (merge) then
                ! Merge this layer with the one above and backtrack.
                Rc(kc) = (Hc(kc)*Rc(kc) + Hf(k,i)*Rf(k,i)) / (Hc(kc) + Hf(k,i))
                Hc(kc) = Hc(kc) + Hf(k,i)
                dzc(kc) = dzc(kc) + dzf(k,i)
                ! Backtrack to remove any convective instabilities above...  Note
                ! that the tolerance is a factor of two larger, to avoid limit how
                ! far back we go.
                do k2=kc,2,-1
                  if (better_est) then
                    merge = ((Rc(k2)-Rc(k2-1)) * ((Hc(k2) * Hc(k2-1))*I_Htot) < tol_merge*drxh_sum)
                  else
                    merge = ((Rc(k2)-Rc(k2-1)) * (Hc(k2)+Hc(k2-1)) < tol_merge*drxh_sum)
                  endif
                  if (merge) then
                    ! Merge the two bottommost layers.  At this point kc = k2.
                    Rc(kc-1) = (Hc(kc)*Rc(kc) + Hc(kc-1)*Rc(kc-1)) / (Hc(kc) + Hc(kc-1))
                    Hc(kc-1) = Hc(kc) + Hc(kc-1)
                    dzc(kc-1) = dzc(kc) + dzc(kc-1)
                    kc = kc - 1
                  else ; exit ; endif
                enddo
              else
                ! Add a new layer to the column.
                kc = kc + 1
                Rc(kc) = Rf(k,i) ; Hc(kc) = Hf(k,i) ; dzc(kc) = dzf(k,i)
              endif
            enddo
            ! At this point there are kc layers and the gprimes should be positive.
            if (nonBous) then
              do K=2,kc
                gprime(K) = H_to_pres * (Rc(k) - Rc(k-1)) / (Rc(k) * Rc(k-1))
              enddo
            else
              do K=2,kc
                gprime(K) = g_Rho0 * (Rc(k)-Rc(k-1))
              enddo
            endif
          endif  ! use_EOS

          !-----------------NOW FIND WAVE SPEEDS---------------------------------------
          ! ig = i + G%idg_offset ; jg = j + G%jdg_offset
          !   Sum the contributions from all of the interfaces to give an over-estimate
          ! of the first-mode wave speed.  Also populate Igl and Igu which are the
          ! non-leading diagonals of the tridiagonal matrix.
          if (kc >= 2) then
            ! initialize speed2_tot
            speed2_tot = 0.0
            if (better_est) then
              H_top(1) = 0.0 ; H_bot(kc+1) = 0.0
              do K=2,kc+1 ; H_top(K) = H_top(K-1) + Hc(k-1) ; enddo
              do K=kc,2,-1 ; H_bot(K) = H_bot(K+1) + Hc(k) ; enddo
              I_Htot = 0.0 ; if (H_top(kc+1) > 0.0) I_Htot = 1.0 / H_top(kc+1)
            endif

            ! Calculate Igu, Igl, depth, and N2 at each interior interface
            ! [excludes surface (K=1) and bottom (K=kc+1)]
            Igl(:) = 0.
            Igu(:) = 0.
            N2(:) = 0.

            do K=2,kc
              Igl(K) = 1.0 / (gprime(K)*Hc(k)) ; Igu(K) = 1.0 / (gprime(K)*Hc(k-1))
              if (nonBous) then
                N2(K) = 2.0*US%L_to_Z**2*gprime(K) * (Hc(k) + Hc(k-1)) / &  ! Units are [T-2 ~> s-2]
                       (dzc(k) + dzc(k-1))**2
              else
                N2(K) = 2.0*US%L_to_Z**2*GV%Z_to_H*gprime(K) / (dzc(k) + dzc(k-1))  ! Units are [T-2 ~> s-2]
              endif
              if (better_est) then
                speed2_tot = speed2_tot + gprime(K)*((H_top(K) * H_bot(K)) * I_Htot)
              else
                speed2_tot = speed2_tot + gprime(K)*(Hc(k-1)+Hc(k))
              endif
            enddo

            ! Set stratification for surface and bottom (setting equal to nearest interface for now)
            N2(1) = N2(2) ; N2(kc+1) = N2(kc)
            ! set bottom stratification
            Nb(i,j) = sqrt(N2(kc+1))

            ! Under estimate the first eigenvalue (overestimate the speed) to start with.
            lam_1 = 1.0 / speed2_tot

            ! init and first guess for mode structure
            mode_struct(:) = 0.
            mode_struct_fder(:) = 0.
            mode_struct(2:kc) = 1. ! Uniform flow, first guess
            modal_structure(:) = 0.
            modal_structure_fder(:) = 0.

            ! Find the first eigen value
            do itt=1,max_itt
              ! calculate the determinant of (A-lam_1*I)
              call tridiag_det(Igu, Igl, 2, kc, lam_1, det, ddet, row_scale=c2_scale)

              ! If possible, use Newton's method iteration to find a new estimate of lam_1
              !det = det_it(itt) ; ddet = ddet_it(itt)
              if ((ddet >= 0.0) .or. (-det > -0.5*lam_1*ddet)) then
                ! lam_1 was not an under-estimate, as intended, so Newton's method
                ! may not be reliable; lam_1 must be reduced, but not by more than half.
                lam_1 = 0.5 * lam_1
                dlam = -lam_1
              else  ! Newton's method is OK.
                dlam = - det / ddet
                lam_1 = lam_1 + dlam
              endif

              call tdma6(kc-1, Igu(2:kc), Igl(2:kc), lam_1, mode_struct(2:kc))
              ! Note that tdma6 changes the units of mode_struct to [L2 T-2 ~> m2 s-2]
              ! apply BC
              mode_struct(1) = 0.
              mode_struct(kc+1) = 0.

              ! renormalization of the integral of the profile
              w2avg = 0.0
              do k=1,kc
                w2avg = w2avg + 0.5*(mode_struct(K)**2+mode_struct(K+1)**2)*Hc(k) ! [H L4 T-4 ~> m5 s-4 or kg m2 s-4]
              enddo
              renorm = sqrt(htot(i)*a_int/w2avg) ! [T2 L-2 ~> s2 m-2]
              do K=1,kc+1 ; mode_struct(K) = renorm * mode_struct(K) ; enddo
              ! after renorm, mode_struct is again [nondim]
              if (abs(dlam) < tol_solve*lam_1) exit
            enddo

            if (lam_1 > 0.0) cn(i,j,1) = 1.0 / sqrt(lam_1)

            ! sign of wave structure is irrelevant, flip to positive if needed
            if (mode_struct(2)<0.) then
              mode_struct(2:kc) = -1. * mode_struct(2:kc)
            endif

            ! vertical derivative of w at interfaces lives on the layer points
            do k=1,kc
              mode_struct_fder(k) = (mode_struct(k) - mode_struct(k+1)) / dzc(k)
            enddo

            ! boundary condition for derivative is no-gradient
            do k=kc+1,nz
              mode_struct_fder(k) = mode_struct_fder(kc)
            enddo

            ! now save maximum value and bottom value
            u_struct_bot(i,j,1) = mode_struct_fder(kc)
            u_struct_max(i,j,1) = maxval(abs(mode_struct_fder(1:kc)))

            ! Calculate terms for vertically integrated energy equation
            do k=1,kc
              mode_struct_fder_sq(k) = mode_struct_fder(k)**2
            enddo
            do K=1,kc+1
              mode_struct_sq(K) = mode_struct(K)**2
            enddo

            ! sum over layers for quantities defined on layer
            do k=1,kc
              int_U2(i,j,1) = int_U2(i,j,1) + mode_struct_fder_sq(k) * Hc(k)
            enddo

            ! vertical integration with Trapezoidal rule for values at interfaces
            do K=1,kc
              int_w2(i,j,1) = int_w2(i,j,1) + 0.5*(mode_struct_sq(K)+mode_struct_sq(K+1)) * Hc(k)
              int_N2w2(i,j,1) = int_N2w2(i,j,1) + 0.5*(mode_struct_sq(K)*N2(K) + &
                                                       mode_struct_sq(K+1)*N2(K+1)) * Hc(k)
            enddo

            ! for w (diag) interpolate onto all interfaces
            call interpolate_column(kc, Hc(1:kc), mode_struct(1:kc+1), &
                                    nz, h(i,j,:), modal_structure(:), .false.)

            ! for u (remap) onto all layers
            call remapping_core_h(CS%remap_CS, kc, Hc(1:kc), mode_struct_fder(1:kc), &
                                  nz, h(i,j,:), modal_structure_fder(:))

            ! write the wave structure
            do k=1,nz+1
              w_struct(i,j,k,1) = modal_structure(k)
            enddo

            do k=1,nz
              u_struct(i,j,k,1) = modal_structure_fder(k)
            enddo

            ! Find other eigen values if c1 is of significant magnitude, > cn_thresh
            nrootsfound = 0    ! number of extra roots found (not including 1st root)
            if ((nmodes > 1) .and. (kc >= nmodes+1) .and. (cn(i,j,1) > CS%c1_thresh)) then
              ! Set the range to look for the other desired eigen values
              ! set min value just greater than the 1st root (found above)
              lamMin = lam_1*(1.0 + tol_solve)
              ! set max value based on a low guess at wavespeed for highest mode
              speed2_min = (reduct_factor*cn(i,j,1)/real(nmodes))**2
              lamMax = 1.0 / speed2_min
              ! set width of interval (not sure about this - BDM)
              lamInc = 0.5*lam_1
              ! set number of intervals within search range
              numint = nint((lamMax - lamMin)/lamInc)

              !   Find intervals containing zero-crossings (roots) of the determinant
              ! that are beyond the first root

              ! find det_l of first interval (det at left endpoint)
              call tridiag_det(Igu, Igl, 2, kc, lamMin, det_l, ddet_l, row_scale=c2_scale)
              ! move interval window looking for zero-crossings************************
              do iint=1,numint
                xr = lamMin + lamInc * iint
                xl = xr - lamInc
                call tridiag_det(Igu, Igl, 2, kc, xr, det_r, ddet_r, row_scale=c2_scale)
                if (det_l*det_r < 0.0) then  ! if function changes sign
                  if (det_l*ddet_l < 0.0) then ! if function at left is headed to zero
                    nrootsfound = nrootsfound + 1
                    xbl(nrootsfound) = xl
                    xbr(nrootsfound) = xr
                  else
                    !   function changes sign but has a local max/min in interval,
                    ! try subdividing interval as many times as necessary (or sub_it_max).
                    ! loop that increases number of subintervals:
                    !call MOM_error(WARNING, "determinant changes sign "// &
                    !            "but has a local max/min in interval; "//&
                    !            "reduce increment in lam.")
                    ! begin subdivision loop -------------------------------------------
                    sub_rootfound = .false. ! initialize
                    do sub_it=1,sub_it_max
                      nsub = 2**sub_it ! number of subintervals; nsub=2,4,8,...
                      ! loop over each subinterval:
                      do sub=1,nsub-1,2 ! only check odds; sub = 1; 1,3; 1,3,5,7; ...
                        xl_sub = xl + lamInc/(nsub)*sub
                        call tridiag_det(Igu, Igl, 2, kc, xl_sub, det_sub, ddet_sub, &
                                         row_scale=c2_scale)
                        if (det_sub*det_r < 0.0) then  ! if function changes sign
                          if (det_sub*ddet_sub < 0.0) then ! if function at left is headed to zero
                            sub_rootfound = .true.
                            nrootsfound = nrootsfound + 1
                            xbl(nrootsfound) = xl_sub
                            xbr(nrootsfound) = xr
                            exit ! exit sub loop
                          endif ! headed toward zero
                        endif ! sign change
                      enddo ! sub-loop
                      if (sub_rootfound) exit ! root has been found, exit sub_it loop
                      !   Otherwise, function changes sign but has a local max/min in one of the
                      ! sub intervals, try subdividing again unless sub_it_max has been reached.
                      if (sub_it == sub_it_max) then
                        call MOM_error(WARNING, "wave_speed: root not found "// &
                                       "after sub_it_max subdivisions of original "// &
                                       "interval.")
                      endif ! sub_it == sub_it_max
                    enddo ! sub_it-loop-------------------------------------------------
                  endif ! det_l*ddet_l < 0.0
                endif ! det_l*det_r < 0.0
                ! exit iint-loop if all desired roots have been found
                if (nrootsfound >= nmodes-1) then
                  ! exit if all additional roots found
                  exit
                elseif (iint == numint) then
                  ! oops, lamMax not large enough - could add code to increase (BDM)
                  ! set unfound modes to zero for now (BDM)
                  !   cn(i,j,nrootsfound+2:nmodes) = 0.0
                else
                  ! else shift interval and keep looking until nmodes or numint is reached
                  det_l = det_r
                  ddet_l = ddet_r
                endif
              enddo ! iint-loop

              ! Use Newton's method to find the roots within the identified windows
              do m=1,nrootsfound ! loop over the root-containing widows (excluding 1st mode)
                lam_n = xbl(m) ! first guess is left edge of window

                ! init and first guess for mode structure
                mode_struct(:) = 0.
                mode_struct_fder(:) = 0.
                mode_struct(2:kc) = 1. ! Uniform flow, first guess
                modal_structure(:) = 0.
                modal_structure_fder(:) = 0.

                do itt=1,max_itt
                  ! calculate the determinant of (A-lam_n*I)
                  call tridiag_det(Igu, Igl, 2, kc, lam_n, det, ddet, row_scale=c2_scale)
                  ! Use Newton's method to find a new estimate of lam_n
                  dlam = - det / ddet
                  lam_n = lam_n + dlam

                  call tdma6(kc-1, Igu(2:kc), Igl(2:kc), lam_n, mode_struct(2:kc))
                  ! Note that tdma6 changes the units of mode_struct to [L2 T-2 ~> m2 s-2]
                  ! apply BC
                  mode_struct(1) = 0.
                  mode_struct(kc+1) = 0.

                  ! renormalization of the integral of the profile
                  w2avg = 0.0
                  do k=1,kc
                    w2avg = w2avg + 0.5*(mode_struct(K)**2+mode_struct(K+1)**2)*Hc(k)
                  enddo
                  renorm = sqrt(htot(i)*a_int/w2avg)
                  do K=1,kc+1 ; mode_struct(K) = renorm * mode_struct(K) ; enddo

                  if (abs(dlam) < tol_solve*lam_1)  exit
                enddo ! itt-loop

                ! calculate nth mode speed
                if (lam_n > 0.0) cn(i,j,m+1) = 1.0 / sqrt(lam_n)

                ! sign is irrelevant, flip to positive if needed
                if (mode_struct(2)<0.) then
                   mode_struct(2:kc) = -1. * mode_struct(2:kc)
                endif

                ! derivative of vertical profile (i.e. dw/dz) is evaluated at the layer point
                do k=1,kc
                  mode_struct_fder(k) = (mode_struct(k) - mode_struct(k+1)) / dzc(k)
                enddo

                ! boundary condition for 1st derivative is no-gradient
                do k=kc+1,nz
                  mode_struct_fder(k) = mode_struct_fder(kc)
                enddo

                ! now save maximum value and bottom value
                u_struct_bot(i,j,m) = mode_struct_fder(kc)
                u_struct_max(i,j,m) = maxval(abs(mode_struct_fder(1:kc)))

                ! Calculate terms for vertically integrated energy equation
                do k=1,kc
                  mode_struct_fder_sq(k) = mode_struct_fder(k)**2
                enddo
                do K=1,kc+1
                  mode_struct_sq(K) = mode_struct(K)**2
                enddo

                ! sum over layers for integral of quantities defined at layer points
                do k=1,kc
                  int_U2(i,j,m) = int_U2(i,j,m) + mode_struct_fder_sq(k) * Hc(k)
                enddo

                ! vertical integration with Trapezoidal rule for quantities on interfaces
                do K=1,kc
                  int_w2(i,j,m) = int_w2(i,j,m) + 0.5*(mode_struct_sq(K)+mode_struct_sq(K+1)) * Hc(k)
                  int_N2w2(i,j,m) = int_N2w2(i,j,m) + 0.5*(mode_struct_sq(K)*N2(K) + &
                                                           mode_struct_sq(K+1)*N2(K+1)) * Hc(k)
                enddo

                ! for w (diag) interpolate onto all interfaces
                call interpolate_column(kc, Hc(1:kc), mode_struct(1:kc+1), &
                                        nz, h(i,j,:), modal_structure(:), .false.)

                ! for u (remap) onto all layers
                call remapping_core_h(CS%remap_CS, kc, Hc(1:kc), mode_struct_fder(1:kc), &
                                      nz, h(i,j,:), modal_structure_fder(:))

                ! write the wave structure
                ! note that m=1 solves for 2nd mode,...
                do k=1,nz+1
                  w_struct(i,j,k,m+1) = modal_structure(k)
                enddo

                do k=1,nz
                  u_struct(i,j,k,m+1) = modal_structure_fder(k)
                enddo

              enddo ! n-loop
            endif ! if nmodes>1 .and. kc>nmodes .and. c1>c1_thresh
          endif ! if more than 2 layers
        endif ! if drxh_sum < 0
      endif ! if not land
    enddo ! i-loop
  enddo ! j-loop

end subroutine wave_speeds

!> Calculate the determinant of a tridiagonal matrix with diagonals a,b-lam,c and its derivative
!! with lam, where lam is constant across rows.  Only the ratio of det to its derivative and their
!! signs are typically used, so internal rescaling by consistent factors are used to avoid
!! over- or underflow.
pure subroutine tridiag_det(a, c, ks, ke, lam, det, ddet, row_scale)
  real, dimension(:), intent(in) :: a     !< Lower diagonal of matrix (first entry unused) [T2 L-2 ~> s2 m-2]
  real, dimension(:), intent(in) :: c     !< Upper diagonal of matrix (last entry unused) [T2 L-2 ~> s2 m-2]
  integer,            intent(in) :: ks    !< Starting index to use in determinant
  integer,            intent(in) :: ke    !< Ending index to use in determinant
  real,               intent(in) :: lam   !< Value subtracted from b [T2 L-2 ~> s2 m-2]
  real,               intent(out):: det   !< Determinant of the matrix in dynamically rescaled units that
                                          !! depend on the number of rows and the cumulative magnitude of
                                          !! det and are therefore difficult to interpret, but the units
                                          !! of det/ddet are always in [T2 L-2 ~> s2 m-2]
  real,               intent(out):: ddet  !< Derivative of determinant with lam in units that are dynamically
                                          !! rescaled along with those of det, such that the units of
                                          !! det/ddet are always in [T2 L-2 ~> s2 m-2]
  real,               intent(in) :: row_scale !< A scaling factor of the rows of the matrix to
                                          !! limit the growth of the determinant [L2 s2 T-2 m-2 ~> 1]

  ! Local variables
  real :: detKm1, detKm2   ! Cumulative value of the determinant for the previous two layers in units
                           ! that vary with the number of layers that have been worked on [various]
  real :: ddetKm1, ddetKm2 ! Derivative of the cumulative determinant with lam for the previous two
                           ! layers [various], but the units of detKm1/ddetKm1 are [T2 L-2 ~> s2 m-2]
  real, parameter :: rescale = 1024.0**4 ! max value of determinant allowed before rescaling [nondim]
  real :: I_rescale ! inverse of rescale [nondim]
  integer :: k      ! row (layer interface) index

  I_rescale = 1.0 / rescale

  detKm1 = 1.0 ; ddetKm1 = 0.0
  det = (a(ks)+c(ks)) - lam ; ddet = -1.0
  do k=ks+1,ke
    ! Shift variables and rescale rows to avoid over- or underflow.
    detKm2 = row_scale*detKm1 ; ddetKm2 = row_scale*ddetKm1
    detKm1 = row_scale*det    ; ddetKm1 = row_scale*ddet

    det =  ((a(k)+c(k))-lam)*detKm1  - (a(k)*c(k-1))*detKm2
    ddet = ((a(k)+c(k))-lam)*ddetKm1 - (a(k)*c(k-1))*ddetKm2 - detKm1

    ! Rescale det & ddet if det is getting too large or too small.
    if (abs(det) > rescale) then
      det = I_rescale*det ; detKm1 = I_rescale*detKm1
      ddet = I_rescale*ddet ; ddetKm1 = I_rescale*ddetKm1
    elseif (abs(det) < I_rescale) then
      det = rescale*det ; detKm1 = rescale*detKm1
      ddet = rescale*ddet ; ddetKm1 = rescale*ddetKm1
    endif
  enddo

end subroutine tridiag_det

!> Initialize control structure for MOM_wave_speed
subroutine wave_speed_init(CS, GV, param_file, use_ebt_mode, mono_N2_column_fraction, mono_N2_depth, &
                           remap_answers_2018, remap_answer_date, better_speed_est, om4_remap_via_sub_cells, &
                           min_speed, wave_speed_tol, c1_thresh)
  type(wave_speed_CS), intent(inout) :: CS  !< Wave speed control struct
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  logical, optional, intent(in) :: use_ebt_mode  !< If true, use the equivalent
                                     !! barotropic mode instead of the first baroclinic mode.
  real,    optional, intent(in) :: mono_N2_column_fraction !< The lower fraction of water column over
                                     !! which N2 is limited as monotonic for the purposes of
                                     !! calculating the vertical modal structure [nondim].
  real,    optional, intent(in) :: mono_N2_depth !< The depth below which N2 is limited
                                     !! as monotonic for the purposes of calculating the
                                     !! vertical modal structure [H ~> m or kg m-2].
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                     !! that recover the remapping answers from 2018.  Otherwise
                                     !! use more robust but mathematically equivalent expressions.
  integer, optional, intent(in) :: remap_answer_date  !< The vintage of the order of arithmetic and expressions
                                      !! to use for remapping.  Values below 20190101 recover the remapping
                                      !! answers from 2018, while higher values use more robust
                                      !! forms of the same remapping expressions.
  logical, optional, intent(in) :: better_speed_est !< If true, use a more robust estimate of the first
                                     !! mode speed as the starting point for iterations.
  logical, optional, intent(in) :: om4_remap_via_sub_cells !< Use the OM4-era ramap_via_sub_cells
                                     !! for calculating the EBT structure
  real,    optional, intent(in) :: min_speed !< If present, set a floor in the first mode speed
                                     !! below which 0 is returned [L T-1 ~> m s-1].
  real,    optional, intent(in) :: wave_speed_tol !< The fractional tolerance for finding the
                                     !! wave speeds [nondim]
  real,    optional, intent(in) :: c1_thresh !< A minimal value of the first mode internal wave speed
                                       !! below which all higher mode speeds are not calculated but are
                                       !! simply reported as 0 [L T-1 ~> m s-1].  A non-negative value
                                       !! must be specified for wave_speeds to be used (but not wave_speed).

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_wave_speed"  ! This module's name.
  !   Do not use an apostrophe in the comments in the block below.  The doxygen preprocessor
  ! treats one in a Fortran comment as the start of a character literal that spans lines, and
  ! if that literal closes between the #ifdef and the #endif, doxygen reports the #endif as
  ! having no matching #if.
#ifdef __NVCOMPILER_OPENMP_GPU
  ! The whole computational domain is worked on at once so that the offloaded loops have as
  ! much parallelism as possible.
  integer, parameter :: default_niblock = 0 ! Default i block size for the wave speed calculations
  integer, parameter :: default_njblock = 0 ! Default j block size for the wave speed calculations
#else
  ! These are a starting point for tuning, chosen so that the columns of a block share the
  ! cache lines that the vertical sweep of each column touches.
  integer, parameter :: default_niblock = 0 ! Default i block size for the wave speed calculations
  integer, parameter :: default_njblock = 1 ! Default j block size for the wave speed calculations
#endif

  CS%initialized = .true.

  ! Write all relevant parameters to the model log.
  call log_version(mdl, version)

  call get_param(param_file, mdl, "WAVE_SPEED_NIBLOCK", CS%niblock, &
                 "The i-direction block size used in the first baroclinic mode wave speed "//&
                 "calculations. If 0, the full computational extent in i is used.", &
                 default=default_niblock, layoutParam=.true.)
  call get_param(param_file, mdl, "WAVE_SPEED_NJBLOCK", CS%njblock, &
                 "The j-direction block size used in the first baroclinic mode wave speed "//&
                 "calculations. If 0, the full computational extent in j is used.", &
                 default=default_njblock, layoutParam=.true.)
  if (CS%niblock < 0) &
    call MOM_error(FATAL, "WAVE_SPEED_NIBLOCK must be nonnegative; "//&
                          "use 0 to select the default block size.")
  if (CS%njblock < 0) &
    call MOM_error(FATAL, "WAVE_SPEED_NJBLOCK must be nonnegative; "//&
                          "use 0 to select the default block size.")

  call wave_speed_set_param(CS, use_ebt_mode=use_ebt_mode, mono_N2_column_fraction=mono_N2_column_fraction, &
                            mono_N2_depth=mono_N2_depth, better_speed_est=better_speed_est, &
                            min_speed=min_speed, wave_speed_tol=wave_speed_tol, &
                            remap_answers_2018=remap_answers_2018, remap_answer_date=remap_answer_date, &
                            c1_thresh=c1_thresh)

  ! The following remapping is only used for wave_speed with pre-2019 answers.
  if (CS%remap_answer_date < 20190101) &
    call initialize_remapping(CS%remap_2018_CS, 'PLM', boundary_extrapolation=.false., &
                              om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                              answer_date=CS%remap_answer_date, &
                              h_neglect=1.0e-30*GV%m_to_H, h_neglect_edge=1.0e-10*GV%m_to_H)

  ! This is used in wave_speeds in all cases, and in wave_speed with newer answers.
  call initialize_remapping(CS%remap_CS, 'PLM', boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                            answer_date=CS%remap_answer_date, &
                            h_neglect=GV%H_subroundoff, h_neglect_edge=GV%H_subroundoff)

end subroutine wave_speed_init

!> Sets internal parameters for MOM_wave_speed
subroutine wave_speed_set_param(CS, use_ebt_mode, mono_N2_column_fraction, mono_N2_depth, remap_answers_2018, &
                                remap_answer_date, better_speed_est, min_speed, wave_speed_tol, c1_thresh)
  type(wave_speed_CS), intent(inout)  :: CS
                                      !< Control structure for MOM_wave_speed
  logical, optional, intent(in) :: use_ebt_mode  !< If true, use the equivalent
                                      !! barotropic mode instead of the first baroclinic mode.
  real,    optional, intent(in) :: mono_N2_column_fraction !< The lower fraction of water column over
                                      !! which N2 is limited as monotonic for the purposes of
                                      !! calculating the vertical modal structure [nondim].
  real,    optional, intent(in) :: mono_N2_depth !< The depth below which N2 is limited
                                      !! as monotonic for the purposes of calculating the
                                      !! vertical modal structure [H ~> m or kg m-2].
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                      !! that recover the remapping answers from 2018.  Otherwise
                                      !! use more robust but mathematically equivalent expressions.
  integer, optional, intent(in) :: remap_answer_date  !< The vintage of the order of arithmetic and expressions
                                      !! to use for remapping.  Values below 20190101 recover the remapping
                                      !! answers from 2018, while higher values use more robust
                                      !! forms of the same remapping expressions.
  logical, optional, intent(in) :: better_speed_est !< If true, use a more robust estimate of the first
                                     !! mode speed as the starting point for iterations.
  real,    optional, intent(in) :: min_speed !< If present, set a floor in the first mode speed
                                     !! below which 0 is returned [L T-1 ~> m s-1].
  real,    optional, intent(in) :: wave_speed_tol !< The fractional tolerance for finding the
                                     !! wave speeds [nondim]
  real,    optional, intent(in) :: c1_thresh !< A minimal value of the first mode internal wave speed
                                       !! below which all higher mode speeds are not calculated but are
                                       !! simply reported as 0 [L T-1 ~> m s-1].  A non-negative value
                                       !! must be specified for wave_speeds to be used (but not wave_speed).

  if (present(use_ebt_mode)) CS%use_ebt_mode = use_ebt_mode
  if (present(mono_N2_column_fraction)) CS%mono_N2_column_fraction = mono_N2_column_fraction
  if (present(mono_N2_depth)) CS%mono_N2_depth = mono_N2_depth
  if (present(remap_answers_2018)) then
    if (remap_answers_2018) then
      CS%remap_answer_date = 20181231
    else
      CS%remap_answer_date = 20190101
    endif
  endif
  if (present(remap_answer_date)) CS%remap_answer_date = remap_answer_date
  if (present(better_speed_est)) CS%better_cg1_est = better_speed_est
  if (present(min_speed)) CS%min_speed2 = min_speed**2
  if (present(wave_speed_tol)) CS%wave_speed_tol = wave_speed_tol
  if (present(c1_thresh)) CS%c1_thresh = c1_thresh

end subroutine wave_speed_set_param

!> \namespace mom_wave_speed

!!
!! Subroutine wave_speed() solves for the first baroclinic mode wave speed.  (It could
!! solve for all the wave speeds, but the iterative approach taken here means
!! that this is not particularly efficient.)
!!
!! If `e(k)` is the perturbation interface height, this means solving for the
!! smallest eigenvalue (`lam` = 1/c^2) of the system
!!
!! \verbatim
!!   -Igu(k)*e(k-1) + (Igu(k)+Igl(k)-lam)*e(k) - Igl(k)*e(k+1) = 0.0
!! \endverbatim
!!
!! with rigid lid boundary conditions e(1) = e(nz+1) = 0.0 giving
!!
!! \verbatim
!!   (Igu(2)+Igl(2)-lam)*e(2) - Igl(2)*e(3) = 0.0
!!   -Igu(nz)*e(nz-1) + (Igu(nz)+Igl(nz)-lam)*e(nz) = 0.0
!! \endverbatim
!!
!! Here
!! \verbatim
!!   Igl(k) = 1.0/(gprime(K)*h(k)) ; Igu(k) = 1.0/(gprime(K)*h(k-1))
!! \endverbatim
!!
!! Alternately, these same eigenvalues can be found from the second smallest
!! eigenvalue of the Montgomery potential (M(k)) calculation:
!!
!! \verbatim
!!   -Igl(k)*M(k-1) + (Igl(k)+Igu(k+1)-lam)*M(k) - Igu(k+1)*M(k+1) = 0.0
!! \endverbatim
!!
!! with rigid lid and flat bottom boundary conditions
!!
!! \verbatim
!!   (Igu(2)-lam)*M(1) - Igu(2)*M(2) = 0.0
!!   -Igl(nz)*M(nz-1) + (Igl(nz)-lam)*M(nz) = 0.0
!! \endverbatim
!!
!! Note that the barotropic mode has been eliminated from the rigid lid
!! interface height equations, hence the matrix is one row smaller.  Without
!! the rigid lid, the top boundary condition is simpler to implement with
!! the M equations.

end module MOM_wave_speed
