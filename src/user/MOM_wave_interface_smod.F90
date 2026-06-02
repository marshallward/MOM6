! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

#include <MOM_memory.h>

!> Interface for surface waves
submodule (MOM_wave_interface) MOM_wave_interface_impl

implicit none

contains

! TODO: I should not need to define this in module and submodule, but NVIDIA
! seems to need both.
!> Return the value of (1 - exp(-x))/x [nondim], using an accurate expression for small values of x.
real function one_minus_exp_x_impl(x)
  real, intent(in) :: x !< The argument of the function ((1 - exp(-x))/x) [nondim]
  real, parameter :: C1_6 = 1.0/6.0  ! A rational fraction [nondim]
  if (abs(x) <= 2.0e-5) then
    ! The Taylor series expression for exp(-x) gives a more accurate expression for 64-bit reals.
    one_minus_exp_x_impl = 1.0 - x * (0.5 - C1_6*x)
  else
    one_minus_exp_x_impl = (1.0 - exp(-x)) / x
  endif
end function one_minus_exp_x_impl


!> Interface to get Langmuir number based on options stored in wave structure
!!
!! Note this can be called with an unallocated Waves pointer, which is okay if we
!!  want the wind-speed only dependent Langmuir number.  Therefore, we need to be
!!  careful about what we try to access here.
module subroutine get_Langmuir_Number( LA, G, GV, US, HBL, ustar, i, j, dz, Waves, &
                                      U_H, V_H, Override_MA)
  !$omp declare target
  type(ocean_grid_type),     intent(in)  :: G     !< Ocean grid structure
  type(verticalGrid_type),   intent(in)  :: GV    !< Ocean vertical grid structure
  real,                      intent(out) :: LA    !< Langmuir number [nondim]
  type(unit_scale_type),     intent(in)  :: US    !< A dimensional unit scaling type
  real,                      intent(in)  :: HBL   !< (Positive) thickness of boundary layer [Z ~> m]
  real,                      intent(in)  :: ustar !< Friction velocity [Z T-1 ~> m s-1]
  integer,                   intent(in)  :: i     !< Meridional index of h-point
  integer,                   intent(in)  :: j     !< Zonal index of h-point
  real, dimension(SZK_(GV)), intent(in)  :: dz    !< Grid layer thickness [Z ~> m]
  type(Wave_parameters_CS),  pointer     :: Waves !< Surface wave control structure.
  real, dimension(SZK_(GV)), &
                   optional, intent(in)  :: U_H   !< Zonal velocity at H point [L T-1 ~> m s-1] or [m s-1]
  real, dimension(SZK_(GV)), &
                   optional, intent(in)  :: V_H   !< Meridional velocity at H point [L T-1 ~> m s-1] or [m s-1]
  logical,         optional, intent(in)  :: Override_MA !< Override to use misalignment in LA
                                                  !! calculation. This can be used if diagnostic
                                                  !! LA outputs are desired that are different than
                                                  !! those used by the dynamical model.

!Local Variables
  real :: Top, Bottom, MidPoint  ! Positions within each layer [Z ~> m]
  real :: Dpt_LASL         ! Averaging depth for Stokes drift [Z ~> m]
  real :: ShearDirection   ! Shear angular direction from atan2 [radians]
  real :: WaveDirection    ! Wave angular direction from atan2 [radians]
  real :: LA_STKx, LA_STKy, LA_STK ! Stokes velocities in [L T-1 ~> m s-1]
  logical :: ContinueLoop, USE_MA
  real, dimension(SZK_(GV)) :: US_H, VS_H ! Profiles of Stokes velocities [L T-1 ~> m s-1]
  real, allocatable :: StkBand_X(:), StkBand_Y(:) ! Stokes drifts by band [L T-1 ~> m s-1]
  integer :: k, BB

  ! Compute averaging depth for Stokes drift (negative)
  Dpt_LASL = -1.0*max(Waves%LA_FracHBL*HBL, Waves%LA_HBL_min)

  USE_MA = Waves%LA_Misalignment
  if (present(Override_MA)) USE_MA = Override_MA

  ! If requesting to use misalignment in the Langmuir number compute the Shear Direction
  if (USE_MA) then
    !*!if (.not.(present(U_H).and.present(V_H))) call MOM_error(FATAL, &
    !*!    "Get_LA_waves requested to consider misalignment, but velocities were not provided.")
    ContinueLoop = .true.
    bottom = 0.0
    do k = 1,GV%ke
      Top = Bottom
      MidPoint = Bottom + 0.5*dz(k)
      Bottom = Bottom + dz(k)

      if (Waves%LA_Misalign_bug) then
        ! Given the sign convention that Dpt_LASL is negative, the next line has a bug.
        if (MidPoint > Dpt_LASL .and. k > 1 .and. ContinueLoop) then
          ShearDirection = atan2(V_H(1)-V_H(k), U_H(1)-U_H(k))
          ContinueLoop = .false.
        endif
      else ! This version avoids the bug in the version above.
        if (MidPoint > abs(Dpt_LASL) .and. (k > 1) .and. ContinueLoop) then
          ShearDirection = atan2(V_H(1)-V_H(k), U_H(1)-U_H(k))
          ContinueLoop = .false.
        endif
      endif
    enddo
  endif

  if (Waves%WaveMethod==TESTPROF) then
    !**!do k = 1,GV%ke
    !**!  US_H(k) = 0.5*(Waves%US_X(I,j,k)+Waves%US_X(I-1,j,k))
    !**!  VS_H(k) = 0.5*(Waves%US_Y(i,J,k)+Waves%US_Y(i,J-1,k))
    !**!enddo
    !**!call Get_SL_Average_Prof( GV, Dpt_LASL, dz, US_H, LA_STKx)
    !**!call Get_SL_Average_Prof( GV, Dpt_LASL, dz, VS_H, LA_STKy)
    !**!LA_STK = sqrt((LA_STKX*LA_STKX) + (LA_STKY*LA_STKY))
  elseif (Waves%WaveMethod==SURFBANDS) then
    !**!allocate(StkBand_X(Waves%NumBands), StkBand_Y(Waves%NumBands))
    !**!do bb = 1,Waves%NumBands
    !**!  StkBand_X(bb) = 0.5*(Waves%STKx0(I,j,bb)+Waves%STKx0(I-1,j,bb))
    !**!  StkBand_Y(bb) = 0.5*(Waves%STKy0(i,J,bb)+Waves%STKy0(i,J-1,bb))
    !**!enddo
    !**!call Get_SL_Average_Band(GV, Dpt_LASL, Waves%NumBands, Waves%WaveNum_Cen, StkBand_X, LA_STKx )
    !**!call Get_SL_Average_Band(GV, Dpt_LASL, Waves%NumBands, Waves%WaveNum_Cen, StkBand_Y, LA_STKy )
    !**!LA_STK = sqrt((LA_STKX**2) + (LA_STKY**2))
    !**!deallocate(StkBand_X, StkBand_Y)
  elseif (Waves%WaveMethod==DHH85) then
    ! Temporarily integrating profile rather than spectrum for simplicity
    do k = 1,GV%ke
      US_H(k) = 0.5*(Waves%US_X(I,j,k)+Waves%US_X(I-1,j,k))
      VS_H(k) = 0.5*(Waves%US_Y(i,J,k)+Waves%US_Y(i,J-1,k))
    enddo
    !***** is this it? *****!
    !**!call Get_SL_Average_Prof( GV, Dpt_LASL, dz, US_H, LA_STKx)
    !**!call Get_SL_Average_Prof( GV, Dpt_LASL, dz, VS_H, LA_STKy)
    !***** is this it? *****!
    LA_STK = sqrt((LA_STKX**2) + (LA_STKY**2))
  elseif (Waves%WaveMethod==LF17) then
    !*!call get_StokesSL_LiFoxKemper(ustar, HBL*Waves%LA_FracHBL, GV, US, Waves, LA_STK, LA)
  elseif (Waves%WaveMethod==Null_WaveMethod) then
    !*!call MOM_error(FATAL, "Get_Langmuir_number called without defining a WaveMethod. "//&
    !*!                      "Suggest to make sure USE_LT is set/overridden to False or choose "//&
    !*!                      "a wave method (or set USE_LA_LI2016 to use statistical waves).")
  endif

  !*!if (.not.(Waves%WaveMethod==LF17)) then
  !*!  ! This expression uses an arbitrary lower bound on Langmuir number.
  !*!  ! We shouldn't expect values lower than this, but there is also no good reason to cap it here
  !*!  ! other than to prevent large enhancements in unconstrained parts of the curve fit parameterizations.
  !*!  LA = max(Waves%La_min, sqrt(US%Z_to_L*ustar / (LA_STK + Waves%La_Stk_backgnd)))
  !*!endif

  !*!if (Use_MA) then
  !*!  WaveDirection = atan2(LA_STKy, LA_STKx)
  !*!  LA = LA / sqrt(max(1.e-8, cos( WaveDirection - ShearDirection)))
  !*!endif

end subroutine get_Langmuir_Number


!> Get SL averaged Stokes drift from Li/FK 17 method
!!
!! Original description:
!! - This function returns the enhancement factor, given the 10-meter
!!   wind [m s-1], friction velocity [m s-1] and the boundary layer depth [m].
!!
!! Update (Jan/25):
!! - Converted from function to subroutine, now returns Langmuir number.
!! - Compute 10m wind internally, so only ustar and hbl need passed to
!!   subroutine.
!!
!! Qing Li, 160606
!! - BGR port from CVMix to MOM6 Jan/25/2017
!! - BGR change output to LA from Efactor
!! - BGR remove u10 input
!! - BGR note: fixed parameter values should be changed to "get_params"
subroutine get_StokesSL_LiFoxKemper(ustar, hbl, GV, US, CS, UStokes_SL, LA)
  real, intent(in)  :: ustar !< water-side surface friction velocity [Z T-1 ~> m s-1].
  real, intent(in)  :: hbl   !< boundary layer depth [Z ~> m].
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer   :: CS  !< Wave parameter Control structure
  real, intent(out) :: UStokes_SL !< Surface layer averaged Stokes drift [L T-1 ~> m s-1]
  real, intent(out) :: LA    !< Langmuir number [nondim]
  ! Local variables
  ! parameters
  real, parameter :: u19p5_to_u10 = 1.075 ! ratio of U19.5 to U10 (Holthuijsen, 2007) [nondim]
  real, parameter :: fm_into_fp = 1.296 ! ratio of mean frequency to peak frequency for
                   ! Pierson-Moskowitz spectrum (Webb, 2011) [nondim]
  real, parameter :: us_to_u10 = 0.0162 ! ratio of surface Stokes drift to U10 [nondim]
  real, parameter :: r_loss = 0.667 ! loss ratio of Stokes transport [nondim]
  real :: UStokes  ! The surface Stokes drift [L T-1 ~> m s-1]
  real :: hm0      ! The significant wave height [Z ~> m]
  real :: fm       ! The mean wave frequency [T-1 ~> s-1]
  real :: fp       ! The peak wave frequency [T-1 ~> s-1]
  real :: kphil    ! A peak wavenumber in the Phillips spectrum [Z-1 ~> m-1]
  real :: kstar    ! A rescaled wavenumber? [Z-1 ~> m-1]
  real :: vstokes  ! The total Stokes transport [Z L T-1 ~> m2 s-1]
  real :: z0       ! The boundary layer depth [Z ~> m]
  real :: z0i      ! The inverse of the boundary layer depth [Z-1 ~> m-1]
  real :: r1, r2, r3, r4  ! Nondimensional ratios [nondim]
  real :: r5       ! A single expression that combines r2 and r4 [nondim]
  real :: root_2kz ! The square root of twice the peak wavenumber times the
                   ! boundary layer depth [nondim]
  real :: u10      ! The 10 m wind speed [L T-1 ~> m s-1]
  real :: PI       ! 3.1415926535... [nondim]

  PI = 4.0*atan(1.0)
  UStokes_sl = 0.0
  LA = 1.e8
  if (ustar > 0.0) then
    ! This code should be revised to minimize the number of divisions and cancel out common factors.

    ! Computing u10 based on u_star and COARE 3.5 relationships
    call ust_2_u10_coare3p5(ustar*sqrt(CS%rho_ocn/CS%rho_air), u10, GV, US, CS)
    ! surface Stokes drift
    UStokes = us_to_u10*u10
    !
    ! significant wave height from Pierson-Moskowitz spectrum (Bouws, 1998)
    hm0 = CS%SWH_from_u10sq * u10**2
    !
    ! peak frequency (PM, Bouws, 1998)
    fp = 0.877 * (US%L_to_Z*GV%g_Earth) / (2.0 * PI * u19p5_to_u10 * u10)
    !
    ! mean frequency
    fm = fm_into_fp * fp
    !
    ! total Stokes transport (a factor r_loss is applied to account
    !  for the effect of directional spreading, multidirectional waves
    !  and the use of PM peak frequency and PM significant wave height
    !  on estimating the Stokes transport)
    vstokes = 0.125 * PI * r_loss * US%Z_to_L * fm * hm0**2
    !
    ! the general peak wavenumber for Phillips' spectrum
    ! (Breivik et al., 2016) with correction of directional spreading
    kphil = 0.176 * UStokes / vstokes

    ! Combining all of the expressions above gives kPhil as the following
    ! where the first two lines are just a constant:
    ! kphil = ((0.176 * us_to_u10 * u19p5_to_u10) / &
    !          (0.5*0.125 * r_loss * fm_into_fp * 0.877 * CS%SWH_from_u10sq**2)) / &
    !         (GV%g_Earth * u10**2)

    ! surface layer
    z0 = abs(hbl)

    if (CS%answer_date < 20230102) then
      z0i = 1.0 / z0

      ! Surface layer averaged Stokes drift with Stokes drift profile
      ! estimated from Phillips' spectrum (Breivik et al., 2016)
      ! The directional spreading effect from Webb and Fox-Kemper, 2015 is also included.
      kstar = kphil * 2.56

      ! Terms 1 to 4, as written in the appendix of Li et al. (2017)
      r1 = ( 0.151 / kphil * z0i - 0.84 ) * &
           ( 1.0 - exp(-2.0 * kphil * z0) )
      r2 = -( 0.84 + 0.0591 / kphil * z0i ) * &
           sqrt( 2.0 * PI * kphil * z0 ) * &
           erfc( sqrt( 2.0 * kphil * z0 ) )
      r3 = ( 0.0632 / kstar * z0i + 0.125 ) * &
           (1.0 - exp(-2.0 * kstar * z0) )
      r4 = ( 0.125 + 0.0946 / kstar * z0i ) * &
           sqrt( 2.0 * PI * kstar * z0) * &
           erfc( sqrt( 2.0 * kstar * z0 ) )
      UStokes_sl = UStokes * (0.715 + r1 + r2 + r3 + r4)
    else
      ! The following is equivalent to the code above, but avoids singularities
      r1 = ( 0.302 - 1.68*(kphil*z0) ) * one_minus_exp_x_impl(2.0 * (kphil * z0))
      r3 = ( 0.1264 + 0.64*(kphil*z0) ) * one_minus_exp_x_impl(5.12 * (kphil * z0))

      root_2kz = sqrt(2.0 * kphil * z0)
      ! r2 = -( 0.84 + 0.0591*2.0 / (root_2kz**2) ) * sqrt(PI) * root_2kz * erfc( root_2kz )
      ! r4 = ( 0.2 + 0.059125*2.0 / (root_2kz**2) ) * sqrt(PI) * root_2kz * erfc( 1.6 * root_2kz )

      ! r5 = r2 + r4 (with a small correction to one coefficient to avoid a singularity when z0 = 0):
      ! The correction leads to <1% relative differences in (r2+r4) for root_2kz > 0.05, but without
      ! it the values of r2 + r4 are qualitatively wrong (>50% errors) for root_2kz < 0.0015 .
      !   It has been verified that these two expressions for r5 are the same to 6 decimal places for
      ! root_2kz  between 1e-10 and 1e-3, but that the first one degrades for smaller values.
      if (root_2kz > 1e-3) then
        r5 = sqrt(PI) * (root_2kz * (-0.84 * erfc(root_2kz) + 0.2 * erfc(1.6*root_2kz)) + &
                         0.1182 * (erfc(1.6*root_2kz) - erfc(root_2kz)) / root_2kz)
      else
        ! It is more accurate to replace erf with the first two terms of its Taylor series
        !  erf(z) = (2/sqrt(pi)) * z * (1. - (1/3)*z**2 + (1/10)*z**4 - (1/42)*z**6 + ...)
        ! and then cancel or combine common terms and drop negligibly small terms.
        r5 = -0.64*sqrt(PI)*root_2kz + (-0.14184 + 1.0839648 * root_2kz**2)
      endif
      UStokes_sl = UStokes * (0.715 + ((r1 + r3) + r5))
    endif

    if (UStokes_sl /= 0.0) LA = sqrt(US%Z_to_L*ustar / UStokes_sl)
  endif

end subroutine Get_StokesSL_LiFoxKemper

!> Get SL Averaged Stokes drift from a Stokes drift Profile
subroutine Get_SL_Average_Prof( GV, AvgDepth, dz, Profile, Average )
  type(verticalGrid_type),  &
       intent(in)   :: GV       !< Ocean vertical grid structure
  real, intent(in)  :: AvgDepth !< Depth to average over (negative) [Z ~> m]
  real, dimension(SZK_(GV)), &
       intent(in)   :: dz       !< Grid thickness [Z ~> m]
  real, dimension(SZK_(GV)), &
       intent(in)   :: Profile  !< Profile of quantity to be averaged in arbitrary units [A]
  real, intent(out) :: Average  !< Output quantity averaged over depth AvgDepth [A]
                                !! (used here for Stokes drift)
  !Local variables
  real :: Top, Bottom ! Depths, negative downward [Z ~> m]
  real :: vsum  ! The depth weighted vertical sum of a quantity [A Z ~> A m]
  integer :: k

  ! Initializing sum
  vsum = 0.0

  ! Integrate
  bottom = 0.0
  do k = 1, GV%ke
    Top = Bottom
    Bottom = Bottom - dz(k)
    if (AvgDepth < Bottom) then ! The whole cell is within H_LA
      vsum = vsum + Profile(k) * dz(k)
    elseif (AvgDepth < Top) then ! A partial cell is within H_LA
      vsum = vsum + Profile(k) * (Top-AvgDepth)
      exit
    else
      exit
    endif
  enddo

  ! Divide by AvgDepth or the depth in the column, whichever is smaller.
  if (abs(AvgDepth) <= abs(Bottom)) then
    Average = vsum / abs(AvgDepth)
  elseif (abs(Bottom) > 0.0) then
    Average = vsum / abs(Bottom)
  else
    Average = 0.0
  endif

end subroutine Get_SL_Average_Prof

!> Get SL averaged Stokes drift from the banded Spectrum method
subroutine Get_SL_Average_Band( GV, AvgDepth, NB, WaveNumbers, SurfStokes, Average )
  type(verticalGrid_type),  &
       intent(in)     :: GV          !< Ocean vertical grid
  real, intent(in)    :: AvgDepth    !< Depth to average over [Z ~> m].
  integer, intent(in) :: NB          !< Number of bands used
  real, dimension(NB), &
       intent(in)     :: WaveNumbers !< Wavenumber corresponding to each band [Z-1 ~> m-1]
  real, dimension(NB), &
       intent(in)     :: SurfStokes  !< Surface Stokes drift for each band [L T-1 ~> m s-1]
  real, intent(out)   :: Average     !< Output average Stokes drift over depth AvgDepth [L T-1 ~> m s-1]

  ! Local variables
  integer :: bb

  ! Loop over bands
  Average = 0.0
  do bb = 1, NB
    ! Factor includes analytical integration of e(2kz)
    !  - divided by (-H_LA) to get average from integral.
    Average = Average + SurfStokes(BB) * &
              (1.-EXP(-abs(AvgDepth * 2.0 * WaveNumbers(BB)))) / &
              abs(AvgDepth * 2.0 * WaveNumbers(BB))

    ! For accuracy when AvgDepth is small change the above to:
    ! Average = Average + SurfStokes(BB) * one_minus_exp_x_impl(abs(AvgDepth * 2.0 * WaveNumbers(BB)))
  enddo

end subroutine Get_SL_Average_Band


!> Computes wind speed from ustar_air based on COARE 3.5 Cd relationship
!! Probably doesn't belong in this module, but it is used here to estimate
!! wind speed for wind-wave relationships.  Should be a fine way to estimate
!! the neutral wind-speed as written here.
subroutine ust_2_u10_coare3p5(USTair, U10, GV, US, CS)
  !$omp declare target
  real, intent(in)                    :: USTair !< Wind friction velocity [Z T-1 ~> m s-1]
  real, intent(out)                   :: U10    !< 10-m neutral wind speed [L T-1 ~> m s-1]
  type(verticalGrid_type), intent(in) :: GV     !< vertical grid type
  type(unit_scale_type),   intent(in) :: US     !< A dimensional unit scaling type
  type(wave_parameters_CS), pointer   :: CS     !< Wave parameter Control structure

  ! Local variables
  real :: z0sm, z0, z0rough  ! Roughness lengths [Z ~> m]
  real :: ten_m_scale ! The 10 m reference height, in rescaled units [Z ~> m]
  real :: I_ten_m_scale ! The inverse of the 10 m reference height, in rescaled units [Z-1 ~> m-1]
  real :: u10a  ! The previous guess for u10 [L T-1 ~> m s-1]
  real :: alpha ! The Charnock coeffient relating the wind friction velocity squared to the
                ! roughness length [nondim]
  real :: Cd    ! The drag coefficient [nondim]
  real :: I_sqrtCd  ! The inverse of the square root of the drag coefficient [nondim]
  real :: I_vonKar  ! The inverse of the von Karman coefficient [nondim]
  integer :: CT

  ! Uses empirical formula for z0 to convert ustar_air to u10 based on the
  !  COARE 3.5 paper (Edson et al., 2013)
  ! alpha=m*U10+b
  ! Note in Edson et al. 2013, eq. 13 m is given as 0.017.  However,
  ! m=0.0017 reproduces the curve in their figure 6.

  !*!if (CS%vonKar < 0.0) call MOM_error(FATAL, &
  !*!  "ust_2_u10_coare3p5 called with a negative value of Waves%vonKar")

  z0sm = 0.11 * CS%nu_air / USTair ! Compute z0smooth from ustar guess
  u10a = 1000.0*US%m_s_to_L_T ! An insanely large upper bound for u10.

  if (CS%answer_date < 20230103) then
    u10 = US%Z_to_L*USTair / sqrt(0.001)  ! Guess for u10
    ten_m_scale = 10.0*US%m_to_Z
    CT=0
    do while (abs(u10a/u10 - 1.) > 0.001)
      CT=CT+1
      u10a = u10
      alpha = min(CS%Charnock_min, CS%Charnock_slope_U10 * u10 + CS%Charnock_intercept)
      z0rough = alpha * (US%Z_to_L*USTair)**2 / GV%g_Earth ! Compute z0rough from ustar guess
      z0 = z0sm + z0rough
      Cd = ( CS%vonKar / log(ten_m_scale / z0) )**2 ! Compute Cd from derived roughness
      u10 = US%Z_to_L*USTair/sqrt(Cd)  ! Compute new u10 from derived Cd, while loop
                             ! ends and checks for convergence...CT counter
                             ! makes sure loop doesn't run away if function
                             ! doesn't converge.  This code was produced offline
                             ! and converged rapidly (e.g. 2 cycles)
                             ! for ustar=0.0001:0.0001:10.
      if (CT>20) then
        u10 = US%Z_to_L*USTair/sqrt(0.0015) ! I don't expect to get here, but just
                                !  in case it will output a reasonable value.
        exit
      endif
    enddo

  else ! Use more efficient expressions that are mathematically equivalent to those above.
    u10 = US%Z_to_L*USTair * sqrt(1000.0) ! First guess for u10.
    ! In the line above 1000 is the inverse of a plausible first guess of the drag coefficient.
    I_vonKar = 1.0 / CS%vonKar
    I_ten_m_scale = 0.1*US%Z_to_m

    do CT=1,20
      if (abs(u10a - u10) <= 0.001*u10) exit   ! Check for convergence.
      u10a = u10
      alpha = min(CS%Charnock_min, CS%Charnock_slope_U10 * u10 + CS%Charnock_intercept)
      z0rough = alpha * (CS%I_g_Earth * USTair**2) ! Compute z0rough from ustar guess
      z0 = z0sm + z0rough
      I_sqrtCd = abs(log(z0 * I_ten_m_scale)) * I_vonKar ! Compute Cd from derived roughness
      u10 = US%Z_to_L*USTair * I_sqrtCd  ! Compute new u10 from the derived Cd.
    enddo

    ! Output a reasonable estimate of u10 if the iteration has not converged. The hard-coded
    ! number 25.82 is 1/sqrt(0.0015) to 4 decimal places, but the exact value should not matter.
    if (abs(u10a - u10) > 0.001*u10) u10 = US%Z_to_L*USTair * 25.82
  endif

end subroutine ust_2_u10_coare3p5


!> \namespace  mom_wave_interface
!!
!! \author Brandon Reichl, 2018.
!!
!! This module should be moved as wave coupling progresses and
!! likely will should mirror the iceberg or sea-ice model set-up.
!!
!! This module is meant to contain the routines to read in and
!! interpret surface wave data for MOM6. In its original form, the
!! capabilities include setting the Stokes drift in the model (from a
!! variety of sources including prescribed, empirical, and input
!! files).  In short order, the plan is to also amend the subroutine
!! to accept Stokes drift information from an external coupler.
!! Eventually, it will be necessary to break this file apart so that
!! general wave information may be stored in the control structure
!! and the Stokes drift effect can be isolated from processes such as
!! sea-state dependent momentum fluxes, gas fluxes, and other wave
!! related air-sea interaction and boundary layer phenomenon.
!!
!! The Stokes drift are stored on the C-grid with the conventional
!! protocol to interpolate to the h-grid to compute Langmuir number,
!! the primary quantity needed for Langmuir turbulence
!! parameterizations in both the ePBL and KPP approach.  This module
!! also computes full 3d Stokes drift profiles, which will be useful
!! if second-order type boundary layer parameterizations are
!! implemented (perhaps via GOTM, work in progress).

end submodule MOM_wave_interface_impl
