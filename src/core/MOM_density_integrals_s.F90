! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Provides integrals of density
submodule (MOM_density_integrals) MOM_density_integrals_s

implicit none

#include <MOM_memory.h>

contains

!> Calls the appropriate subroutine to calculate analyti
!> Compute pressure gradient force integrals by quadrature for the case where
!! T and S are linear profiles.
module subroutine int_density_dz_generic_plm(k, tv, T_t, T_b, S_t, S_b, e, rho_ref, &
                                      rho_0, G_e, dz_subroundoff, bathyT, HI, GV, EOS, US, use_stanley_eos, dpa, &
                                      intz_dpa, intx_dpa, inty_dpa, MassWghtInterp, &
                                      use_inaccurate_form, Z_0p, MassWghtInterpVanOnly, h_nv)
  integer,              intent(in)  :: k   !< Layer index to calculate integrals for
  type(hor_index_type), intent(in)  :: HI  !< Ocean horizontal index structures for the input arrays
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(thermo_var_ptrs), intent(in) :: tv  !< Thermodynamic variables
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_t !< Potential temperature at the cell top [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_b !< Potential temperature at the cell bottom [C ~> degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_t !< Salinity at the cell top [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_b !< Salinity at the cell bottom [S ~> ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)+1), &
                        intent(in)  :: e   !< Height of interfaces [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is subtracted
                                           !! out to reduce the magnitude of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3], that is used to calculate
                                           !! the pressure (as p~=-z*rho_0*G_e) used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,                 intent(in)  :: dz_subroundoff !< A minuscule thickness change [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  logical,              intent(in) :: use_stanley_eos !< If true, turn on Stanley SGS T variance parameterization
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dpa !< The change in the pressure anomaly across the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer of
                                           !! the pressure anomaly relative to the anomaly at the
                                           !! top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the y grid spacing [R L2 T-2 ~> Pa]
  integer,    optional, intent(in)  :: MassWghtInterp !< A flag indicating whether and how to use
                                           !! mass weighting to interpolate T/S in integrals
  logical,    optional, intent(in)  :: use_inaccurate_form !< If true, uses an inaccurate form of
                                           !! density anomalies, as was used prior to March 2018.
  logical,    optional, intent(in)  :: MassWghtInterpVanOnly !< If true, does not do mass weighting
                                           !! of T/S unless one side smaller than h_nv (i.e. vanished)
  real,       optional, intent(in)  :: h_nv !< Nonvanished height [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

! This subroutine calculates (by numerical quadrature) integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumption here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! It is assumed that the salinity and temperature profiles are linear in the
! vertical. The top and bottom values within each layer are provided and
! a linear interpolation is used to compute intermediate values.

  ! Local variables
  real :: T5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Temperatures along a line of subgrid locations [C ~> degC]
  real :: S5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Salinities along a line of subgrid locations [S ~> ppt]
  real :: T25((5*HI%iscB+1):(5*(HI%iecB+2))) ! SGS temperature variance along a line of subgrid
                                             ! locations [C2 ~> degC2]
  real :: TS5((5*HI%iscB+1):(5*(HI%iecB+2))) ! SGS temp-salt covariance along a line of subgrid
                                             ! locations [C S ~> degC ppt]
  real :: S25((5*HI%iscB+1):(5*(HI%iecB+2))) ! SGS salinity variance along a line of subgrid locations [S2 ~> ppt2]
  real :: p5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Pressures along a line of subgrid locations [R L2 T-2 ~> Pa]
  real :: r5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Densities anomalies along a line of subgrid
                                             ! locations [R ~> kg m-3]
  real :: u5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Densities anomalies along a line of subgrid locations
                                             ! (used for inaccurate form) [R ~> kg m-3]
  real :: T15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Temperatures at an array of subgrid locations [C ~> degC]
  real :: S15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Salinities at an array of subgrid locations [S ~> ppt]
  real :: T215((15*HI%iscB+1):(15*(HI%iecB+1))) ! SGS temperature variance along a line of subgrid
                                                ! locations [C2 ~> degC2]
  real :: TS15((15*HI%iscB+1):(15*(HI%iecB+1))) ! SGS temp-salt covariance along a line of subgrid
                                                ! locations [C S ~> degC ppt]
  real :: S215((15*HI%iscB+1):(15*(HI%iecB+1))) ! SGS salinity variance along a line of subgrid
                                                ! locations [S2 ~> ppt2]
  real :: p15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Pressures at an array of subgrid locations [R L2 T-2 ~> Pa]
  real :: r15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Densities at an array of subgrid locations [R ~> kg m-3]
  real :: wt_t(5), wt_b(5)          ! Top and bottom weights [nondim]
  real :: rho_anom                  ! A density anomaly [R ~> kg m-3]
  real :: w_left, w_right           ! Left and right weights [nondim]
  real :: intz(5)    ! The gravitational acceleration times the integrals of density
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa]
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant [nondim]
  real :: GxRho      ! The product of the gravitational acceleration and reference density [R L2 Z-1 T-2 ~> Pa m-1]
  real :: dz(HI%iscB:HI%iecB+1)   ! Layer thicknesses at tracer points [Z ~> m]
  real :: dz_x(5,HI%iscB:HI%iecB) ! Layer thicknesses along an x-line of subgrid locations [Z ~> m]
  real :: dz_y(5,HI%isc:HI%iec)   ! Layer thicknesses along a y-line of subgrid locations [Z ~> m]
  real :: massWeightToggle          ! A non-dimensional toggle factor for near-bottom mass weighting (0 or 1) [nondim]
  real :: TopWeightToggle           ! A non-dimensional toggle factor for near-surface mass weighting (0 or 1) [nondim]
  real :: massWeightNVonlyToggle    ! A non-dimensional toggle factor for only using mass weighting
                                    ! if at least one side vanished (0 or 1) [nondim]
  real :: Ttl, Tbl, Ttr, Tbr        ! Temperatures at the velocity cell corners [C ~> degC]
  real :: Stl, Sbl, Str, Sbr        ! Salinities at the velocity cell corners [S ~> ppt]
  real :: z0pres(HI%isd:HI%ied,HI%jsd:HI%jed) ! The height at which the pressure is zero [Z ~> m]
  real :: hWght                     ! A topographically limited thickness weight [Z ~> m]
  real :: hWghtTop                  ! An ice draft limited thickness weight [Z ~> m]
  real :: hL, hR                    ! Thicknesses to the left and right [Z ~> m]
  real :: iDenom                    ! The denominator of the thickness weight expressions [Z-2 ~> m-2]
  real :: h_nonvanished             ! nonvanished height [Z ~> m]
  logical :: use_rho_ref ! Pass rho_ref to the equation of state for more accurate calculation
                         ! of density anomalies.
  logical :: use_varT, use_varS, use_covarTS ! Logicals for SGS variances fields
  integer, dimension(2) :: EOSdom_h5  ! The 5-point h-point i-computational domain for the equation of state
  integer, dimension(2) :: EOSdom_q15 ! The 3x5-point q-point i-computational domain for the equation of state
  integer, dimension(2) :: EOSdom_h15 ! The 3x5-point h-point i-computational domain for the equation of state
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n, pos

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  GxRho = G_e * rho_0
  if (present(Z_0p)) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      z0pres(i,j) = Z_0p(i,j)
    enddo ; enddo
  else
    z0pres(:,:) = 0.0
  endif
  massWeightToggle = 0. ; TopWeightToggle = 0.
  if (present(MassWghtInterp)) then
    if (BTEST(MassWghtInterp, 0)) massWeightToggle = 1.
    if (BTEST(MassWghtInterp, 1)) TopWeightToggle = 1.
  endif
  massWeightNVonlyToggle = 1.
  if (present(MassWghtInterpVanOnly)) then
    if (MassWghtInterpVanOnly) massWeightNVonlyToggle = 0.
  endif
  h_nonvanished = 0.
  if (present(h_nv)) then
    h_nonvanished = h_nv
  endif
  use_rho_ref = .true.
  if (present(use_inaccurate_form)) use_rho_ref = .not. use_inaccurate_form

  use_varT = .false. !ensure initialized
  use_covarTS = .false.
  use_varS = .false.
  if (use_stanley_eos) then
    use_varT = associated(tv%varT)
    use_covarTS = associated(tv%covarTS)
    use_varS = associated(tv%varS)
  endif

  T25(:) = 0.
  TS5(:) = 0.
  S25(:) = 0.
  T215(:) = 0.
  TS15(:) = 0.
  S215(:) = 0.

  do n = 1, 5
    wt_t(n) = 0.25 * real(5-n)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  ! Set the loop ranges for equation of state calculations at various points.
  EOSdom_h5(1) = 1 ; EOSdom_h5(2) = 5*(Ieq-Isq+2)
  EOSdom_q15(1) = 1 ; EOSdom_q15(2) = 15*(Ieq-Isq+1)
  EOSdom_h15(1) = 1 ; EOSdom_h15(2) = 15*(HI%iec-HI%isc+1)

  ! 1. Compute vertical integrals
  do j=Jsq,Jeq+1
    do i = Isq,Ieq+1
      dz(i) = e(i,j,K) - e(i,j,K+1)
      do n=1,5
        p5(i*5+n) = -GxRho*((e(i,j,K) - z0pres(i,j)) - 0.25*real(n-1)*dz(i))
        ! Salinity and temperature points are linearly interpolated
        S5(i*5+n) = wt_t(n) * S_t(i,j,k) + wt_b(n) * S_b(i,j,k)
        T5(i*5+n) = wt_t(n) * T_t(i,j,k) + wt_b(n) * T_b(i,j,k)
      enddo
      if (use_varT) T25(i*5+1:i*5+5) = tv%varT(i,j,k)
      if (use_covarTS) TS5(i*5+1:i*5+5) = tv%covarTS(i,j,k)
      if (use_varS) S25(i*5+1:i*5+5) = tv%varS(i,j,k)
    enddo
    if (use_Stanley_eos) then
      call calculate_density(T5, S5, p5, T25, TS5, S25, r5, EOS, EOSdom_h5, rho_ref=rho_ref)
    else
      if (use_rho_ref) then
        call calculate_density(T5, S5, p5, r5, EOS, EOSdom_h5, rho_ref=rho_ref)
      else
        call calculate_density(T5, S5, p5, r5, EOS, EOSdom_h5)
        u5(:) = r5(:) - rho_ref
      endif
    endif

    if (use_rho_ref) then
      do i=Isq,Ieq+1
        ! Use Boole's rule to estimate the pressure anomaly change.
        rho_anom = C1_90*(7.0*(r5(i*5+1)+r5(i*5+5)) + 32.0*(r5(i*5+2)+r5(i*5+4)) + 12.0*r5(i*5+3))
        dpa(i,j) = G_e*dz(i)*rho_anom
        if (present(intz_dpa)) then
          ! Use a Boole's-rule-like fifth-order accurate estimate of
          ! the double integral of the pressure anomaly.
          intz_dpa(i,j) = 0.5*G_e*dz(i)**2 * &
                  (rho_anom - C1_90*(16.0*(r5(i*5+4)-r5(i*5+2)) + 7.0*(r5(i*5+5)-r5(i*5+1))) )
        endif
      enddo
    else
      do i=Isq,Ieq+1
        ! Use Boole's rule to estimate the pressure anomaly change.
        rho_anom = C1_90*(7.0*(r5(i*5+1)+r5(i*5+5)) + 32.0*(r5(i*5+2)+r5(i*5+4)) + 12.0*r5(i*5+3)) &
                   - rho_ref
        dpa(i,j) = G_e*dz(i)*rho_anom
        if (present(intz_dpa)) then
          ! Use a Boole's-rule-like fifth-order accurate estimate of
          ! the double integral of the pressure anomaly.
          intz_dpa(i,j) = 0.5*G_e*dz(i)**2 * &
                  (rho_anom - C1_90*(16.0*(u5(i*5+4)-u5(i*5+2)) + 7.0*(u5(i*5+5)-u5(i*5+1))) )
        endif
      enddo
    endif
  enddo ! end loops on j

  ! 2. Compute horizontal integrals in the x direction
  if (present(intx_dpa)) then ; do j=HI%jsc,HI%jec
    do I=Isq,Ieq
      ! Corner values of T and S
      ! hWght is the distance measure by which the cell is violation of
      ! hydrostatic consistency. For large hWght we bias the interpolation
      ! of T,S along the top and bottom integrals, almost like thickness
      ! weighting.
      ! Note: To work in terrain following coordinates we could offset
      ! this distance by the layer thickness to replicate other models.
      hWght = massWeightToggle * &
              max(0., -bathyT(i,j)-e(i+1,j,K), -bathyT(i+1,j)-e(i,j,K))
      ! CY: The below code just uses top interface, which may be bad in high res open ocean
      ! We want something like if (pa(i+1,k+1)<pa(i,1)) or (pa(i+1,1) <pa(i,k+1)) then...
      ! but pressures are not passed through to this submodule, and tv just has surface press.
      !if ((p(i+1,j,k+1)<p(i,j,1)).or.(tv%p(i+1,j,k+1)<tv%p(i,j,1))) then
      hWghtTop = TopWeightToggle * &
              max(0., e(i+1,j,K+1)-e(i,j,1), e(i,j,K+1)-e(i+1,j,1))
      !else ! pressure criteria not activated
      !  hWghtTop = 0.
      !endif
      ! Set it to be max of the bottom and top hWghts:
      hWght = max(hWght, hWghtTop)
      ! If both sides are nonvanished, then set it back to zero.
      if (((e(i,j,K) - e(i,j,K+1)) > h_nonvanished) .and. ((e(i+1,j,K) - e(i+1,j,K+1)) > h_nonvanished)) then
        hWght = massWeightNVonlyToggle * hWght
      endif
      if (hWght > 0.) then
        hL = (e(i,j,K) - e(i,j,K+1)) + dz_subroundoff
        hR = (e(i+1,j,K) - e(i+1,j,K+1)) + dz_subroundoff
        hWght = hWght * ( (hL-hR)/(hL+hR) )**2
        iDenom = 1./( hWght*(hR + hL) + hL*hR )
        Ttl = ( (hWght*hR)*T_t(i+1,j,k) + (hWght*hL + hR*hL)*T_t(i,j,k) ) * iDenom
        Ttr = ( (hWght*hL)*T_t(i,j,k) + (hWght*hR + hR*hL)*T_t(i+1,j,k) ) * iDenom
        Tbl = ( (hWght*hR)*T_b(i+1,j,k) + (hWght*hL + hR*hL)*T_b(i,j,k) ) * iDenom
        Tbr = ( (hWght*hL)*T_b(i,j,k) + (hWght*hR + hR*hL)*T_b(i+1,j,k) ) * iDenom
        Stl = ( (hWght*hR)*S_t(i+1,j,k) + (hWght*hL + hR*hL)*S_t(i,j,k) ) * iDenom
        Str = ( (hWght*hL)*S_t(i,j,k) + (hWght*hR + hR*hL)*S_t(i+1,j,k) ) * iDenom
        Sbl = ( (hWght*hR)*S_b(i+1,j,k) + (hWght*hL + hR*hL)*S_b(i,j,k) ) * iDenom
        Sbr = ( (hWght*hL)*S_b(i,j,k) + (hWght*hR + hR*hL)*S_b(i+1,j,k) ) * iDenom
      else
        Ttl = T_t(i,j,k) ; Tbl = T_b(i,j,k) ; Ttr = T_t(i+1,j,k) ; Tbr = T_b(i+1,j,k)
        Stl = S_t(i,j,k) ; Sbl = S_b(i,j,k) ; Str = S_t(i+1,j,k) ; Sbr = S_b(i+1,j,k)
      endif

      do m=2,4
        w_left = wt_t(m) ; w_right = wt_b(m)
        dz_x(m,i) = (w_left*(e(i,j,K) - e(i,j,K+1))) + (w_right*(e(i+1,j,K) - e(i+1,j,K+1)))

        ! Salinity and temperature points are linearly interpolated in
        ! the horizontal. The subscript (1) refers to the top value in
        ! the vertical profile while subscript (5) refers to the bottom
        ! value in the vertical profile.
        pos = i*15+(m-2)*5
        T15(pos+1) = (w_left*Ttl) + (w_right*Ttr)
        T15(pos+5) = (w_left*Tbl) + (w_right*Tbr)

        S15(pos+1) = (w_left*Stl) + (w_right*Str)
        S15(pos+5) = (w_left*Sbl) + (w_right*Sbr)

        p15(pos+1) = -GxRho * ((w_left*(e(i,j,K)-z0pres(i,j))) + (w_right*(e(i+1,j,K)-z0pres(i+1,j))))

        ! Pressure
        do n=2,5
          p15(pos+n) = p15(pos+n-1) + GxRho*0.25*dz_x(m,i)
        enddo

        ! Salinity and temperature (linear interpolation in the vertical)
        do n=2,4
          S15(pos+n) = wt_t(n) * S15(pos+1) + wt_b(n) * S15(pos+5)
          T15(pos+n) = wt_t(n) * T15(pos+1) + wt_b(n) * T15(pos+5)
        enddo
        if (use_varT) T215(pos+1:pos+5) = (w_left*tv%varT(i,j,k)) + (w_right*tv%varT(i+1,j,k))
        if (use_covarTS) TS15(pos+1:pos+5) = (w_left*tv%covarTS(i,j,k)) + (w_right*tv%covarTS(i+1,j,k))
        if (use_varS) S215(pos+1:pos+5) = (w_left*tv%varS(i,j,k)) + (w_right*tv%varS(i+1,j,k))
      enddo
    enddo

    if (use_stanley_eos) then
      call calculate_density(T15, S15, p15, T215, TS15, S215, r15, EOS, EOSdom_q15, rho_ref=rho_ref)
    else
      if (use_rho_ref) then
        call calculate_density(T15, S15, p15, r15, EOS, EOSdom_q15, rho_ref=rho_ref)
      else
        call calculate_density(T15, S15, p15, r15, EOS, EOSdom_q15)
      endif
    endif

    do I=Isq,Ieq
      intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)

      ! Use Boole's rule to estimate the pressure anomaly change.
      if (use_rho_ref) then
        do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = (G_e*dz_x(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + 32.0*(r15(pos+2)+r15(pos+4)) + &
                            12.0*r15(pos+3)) ))
        enddo
      else
        do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = (G_e*dz_x(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + 32.0*(r15(pos+2)+r15(pos+4)) + &
                            12.0*r15(pos+3)) - rho_ref ))
        enddo
      endif
      ! Use Boole's rule to integrate the bottom pressure anomaly values in x.
      intx_dpa(I,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    enddo
  enddo ; endif

  ! 3. Compute horizontal integrals in the y direction
  if (present(inty_dpa)) then ; do J=Jsq,Jeq
    do i=HI%isc,HI%iec
    ! Corner values of T and S
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting.
    ! Note: To work in terrain following coordinates we could offset
    ! this distance by the layer thickness to replicate other models.
      hWght = massWeightToggle * &
              max(0., -bathyT(i,j)-e(i,j+1,K), -bathyT(i,j+1)-e(i,j,K))
      ! CY: The below code just uses top interface, which may be bad in high res open ocean
      ! We want something like if (pa(j+1,k+1)<pa(j,1)) or (pa(j+1,1) <pa(i,j,k+1)) then...
      ! but pressures are not passed through to this submodule, and tv just has surface press.
      !if ((p(i,j+1,k+1)<p(i,j,1)).or.(tv%p(i,j+1,k+1)<tv%p(i,j,1))) then
      hWghtTop = TopWeightToggle * &
              max(0., e(i,j+1,K+1)-e(i,j,1), e(i,j,K+1)-e(i,j+1,1))
      !else ! pressure criteria not activated
      !  hWghtTop = 0.
      !endif
      ! Set it to be max of the bottom and top hWghts:
      hWght = max(hWght, hWghtTop)
      ! If both sides are nonvanished, then set it back to zero.
      if (((e(i,j,K) - e(i,j,K+1)) > h_nonvanished) .and. ((e(i,j+1,K) - e(i,j+1,K+1)) > h_nonvanished)) then
        hWght = massWeightNVonlyToggle * hWght
      endif

      if (hWght > 0.) then
        hL = (e(i,j,K) - e(i,j,K+1)) + dz_subroundoff
        hR = (e(i,j+1,K) - e(i,j+1,K+1)) + dz_subroundoff
        hWght = hWght * ( (hL-hR)/(hL+hR) )**2
        iDenom = 1./( hWght*(hR + hL) + hL*hR )
        Ttl = ( (hWght*hR)*T_t(i,j+1,k) + (hWght*hL + hR*hL)*T_t(i,j,k) ) * iDenom
        Ttr = ( (hWght*hL)*T_t(i,j,k) + (hWght*hR + hR*hL)*T_t(i,j+1,k) ) * iDenom
        Tbl = ( (hWght*hR)*T_b(i,j+1,k) + (hWght*hL + hR*hL)*T_b(i,j,k) ) * iDenom
        Tbr = ( (hWght*hL)*T_b(i,j,k) + (hWght*hR + hR*hL)*T_b(i,j+1,k) ) * iDenom
        Stl = ( (hWght*hR)*S_t(i,j+1,k) + (hWght*hL + hR*hL)*S_t(i,j,k) ) * iDenom
        Str = ( (hWght*hL)*S_t(i,j,k) + (hWght*hR + hR*hL)*S_t(i,j+1,k) ) * iDenom
        Sbl = ( (hWght*hR)*S_b(i,j+1,k) + (hWght*hL + hR*hL)*S_b(i,j,k) ) * iDenom
        Sbr = ( (hWght*hL)*S_b(i,j,k) + (hWght*hR + hR*hL)*S_b(i,j+1,k) ) * iDenom
      else
        Ttl = T_t(i,j,k) ; Tbl = T_b(i,j,k) ; Ttr = T_t(i,j+1,k) ; Tbr = T_b(i,j+1,k)
        Stl = S_t(i,j,k) ; Sbl = S_b(i,j,k) ; Str = S_t(i,j+1,k) ; Sbr = S_b(i,j+1,k)
      endif

      do m=2,4
        w_left = wt_t(m) ; w_right = wt_b(m)
        dz_y(m,i) = (w_left*(e(i,j,K) - e(i,j,K+1))) + (w_right*(e(i,j+1,K) - e(i,j+1,K+1)))

        ! Salinity and temperature points are linearly interpolated in
        ! the horizontal. The subscript (1) refers to the top value in
        ! the vertical profile while subscript (5) refers to the bottom
        ! value in the vertical profile.
        pos = i*15+(m-2)*5
        T15(pos+1) = (w_left*Ttl) + (w_right*Ttr)
        T15(pos+5) = (w_left*Tbl) + (w_right*Tbr)

        S15(pos+1) = (w_left*Stl) + (w_right*Str)
        S15(pos+5) = (w_left*Sbl) + (w_right*Sbr)

        p15(pos+1) = -GxRho * ((w_left*(e(i,j,K)-z0pres(i,j))) + (w_right*(e(i,j+1,K)-z0pres(i,j+1))))

        ! Pressure
        do n=2,5
          p15(pos+n) = p15(pos+n-1) + GxRho*0.25*dz_y(m,i)
        enddo

        ! Salinity and temperature (linear interpolation in the vertical)
        do n=2,4
          S15(pos+n) = wt_t(n) * S15(pos+1) + wt_b(n) * S15(pos+5)
          T15(pos+n) = wt_t(n) * T15(pos+1) + wt_b(n) * T15(pos+5)
        enddo
        if (use_varT) T215(pos+1:pos+5) = (w_left*tv%varT(i,j,k)) + (w_right*tv%varT(i,j+1,k))
        if (use_covarTS) TS15(pos+1:pos+5) = (w_left*tv%covarTS(i,j,k)) + (w_right*tv%covarTS(i,j+1,k))
        if (use_varS) S215(pos+1:pos+5) = (w_left*tv%varS(i,j,k)) + (w_right*tv%varS(i,j+1,k))
      enddo
    enddo

    if (use_stanley_eos) then
      call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                             T215(15*HI%isc+1:), TS15(15*HI%isc+1:), S215(15*HI%isc+1:), &
                             r15(15*HI%isc+1:), EOS, EOSdom_h15, rho_ref=rho_ref)
    else
      if (use_rho_ref) then
        call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                               r15(15*HI%isc+1:), EOS, EOSdom_h15, rho_ref=rho_ref)
      else
        call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                               r15(15*HI%isc+1:), EOS, EOSdom_h15)
      endif
    endif

    do i=HI%isc,HI%iec
      intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)

      ! Use Boole's rule to estimate the pressure anomaly change.
      if (use_rho_ref) then
        do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = (G_e*dz_y(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + &
                                           32.0*(r15(pos+2)+r15(pos+4)) + &
                                           12.0*r15(pos+3)) ))
        enddo
      else
        do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = (G_e*dz_y(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + &
                                           32.0*(r15(pos+2)+r15(pos+4)) + &
                                           12.0*r15(pos+3)) - rho_ref ))
        enddo
      endif
      ! Use Boole's rule to integrate the values.
      inty_dpa(i,J) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    enddo
  enddo ; endif

end subroutine int_density_dz_generic_plm

end submodule MOM_density_integrals_s

!> \namespace mom_density_integrals
!!
