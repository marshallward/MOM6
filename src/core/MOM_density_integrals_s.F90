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
  real :: GxRho      ! The product of the gravitational acceleration and reference density [R L2 Z-1 T-2 ~> Pa m-1]
  real :: z0pres(HI%isd:HI%ied,HI%jsd:HI%jed) ! The height at which the pressure is zero [Z ~> m]
  real :: massWeightToggle          ! A non-dimensional toggle factor for near-bottom mass weighting (0 or 1) [nondim]
  real :: TopWeightToggle           ! A non-dimensional toggle factor for near-surface mass weighting (0 or 1) [nondim]
  real :: massWeightNVonlyToggle    ! A non-dimensional toggle factor for only using mass weighting
                                    ! if at least one side vanished (0 or 1) [nondim]
  real :: h_nonvanished             ! nonvanished height [Z ~> m]
  logical :: use_rho_ref ! Pass rho_ref to the equation of state for more accurate calculation
                         ! of density anomalies.
  logical :: use_varT, use_varS, use_covarTS ! Logicals for SGS variances fields
  integer :: Isq, Ieq, Jsq, Jeq, i, j
  integer :: TILE_SIZE_X, TILE_SIZE_Y

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

  ! 1. Compute vertical integrals
  TILE_SIZE_X = Ieq+1-Isq+1 ; TILE_SIZE_Y = Jeq+1-Jsq+1
  call generic_plm_update_dpa(TILE_SIZE_X, TILE_SIZE_Y, k, tv, T_t, T_b, S_t, S_b, e, &
                              z0pres, dpa, intz_dpa, GxRho, G_e, &
                              use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS, &
                              rho_ref, EOS, HI, GV)

  ! 2. Compute horizontal integrals in the x direction
  if (present(intx_dpa)) then
    TILE_SIZE_X = Ieq-Isq+1 ; TILE_SIZE_Y = HI%jec-HI%jsc+1
    call generic_plm_update_intx_dpa(TILE_SIZE_X, TILE_SIZE_Y, k, tv, T_t, T_b, S_t, S_b, e, &
                                     bathyT, z0pres, dpa, intx_dpa, GxRho, G_e, dz_subroundoff, &
                                     massWeightToggle, TopWeightToggle, massWeightNVonlyToggle, h_nonvanished, &
                                     use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS, &
                                     rho_ref, EOS, HI, GV)
  endif

  ! 3. Compute horizontal integrals in the y direction
  if (present(inty_dpa)) then
    TILE_SIZE_X = HI%iec-HI%isc+1 ; TILE_SIZE_Y = Jeq-Jsq+1
    call generic_plm_update_inty_dpa(TILE_SIZE_X, TILE_SIZE_Y, k, tv, T_t, T_b, S_t, S_b, e, &
                                     bathyT, z0pres, dpa, inty_dpa, GxRho, G_e, dz_subroundoff, &
                                     massWeightToggle, TopWeightToggle, massWeightNVonlyToggle, h_nonvanished, &
                                     use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS, &
                                     rho_ref, EOS, HI, GV)
  endif

end subroutine int_density_dz_generic_plm

subroutine generic_plm_update_dpa(TILE_SIZE_X, TILE_SIZE_Y, k, tv, T_t, T_b, S_t, S_b, e, &
                                  z0pres, dpa, intz_dpa, GxRho, G_e, &
                                  use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS, &
                                  rho_ref, EOS, HI, GV)
  integer,              intent(in)  :: TILE_SIZE_X, TILE_SIZE_Y
  integer,              intent(in)  :: k
  type(hor_index_type), intent(in)  :: HI
  type(verticalGrid_type), intent(in) :: GV
  type(thermo_var_ptrs), intent(in) :: tv
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), intent(in) :: T_t, T_b, S_t, S_b
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)+1), intent(in) :: e
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), intent(in) :: z0pres
  real, dimension(SZI_(HI),SZJ_(HI)), intent(inout) :: dpa
  real, dimension(SZI_(HI),SZJ_(HI)), optional, intent(inout) :: intz_dpa
  real,                 intent(in)  :: GxRho, G_e, rho_ref
  logical,              intent(in)  :: use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS
  type(EOS_type),       intent(in)  :: EOS

  real :: T5(5*TILE_SIZE_X,TILE_SIZE_Y)
  real :: S5(5*TILE_SIZE_X,TILE_SIZE_Y)
  real :: T25(5*TILE_SIZE_X,TILE_SIZE_Y)
  real :: TS5(5*TILE_SIZE_X,TILE_SIZE_Y)
  real :: S25(5*TILE_SIZE_X,TILE_SIZE_Y)
  real :: p5(5*TILE_SIZE_X,TILE_SIZE_Y)
  real :: r5(5*TILE_SIZE_X,TILE_SIZE_Y)
  real :: u5(5*TILE_SIZE_X,TILE_SIZE_Y)
  real :: wt_t(5), wt_b(5)
  real :: rho_anom, dz
  real, parameter :: C1_90 = 1.0/90.0
  integer, dimension(2,2) :: EOSdom_h5
  integer :: Isq, Ieq, Jsq, Jeq, i, j, n, jstart, jend, istart, iend, ii, jj

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  T25(:,:) = 0.
  TS5(:,:) = 0.
  S25(:,:) = 0.

  do n = 1, 5
    wt_t(n) = 0.25 * real(5-n)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  do jstart=Jsq,Jeq+1,TILE_SIZE_Y ; do istart=Isq,Ieq+1,TILE_SIZE_X
    jend = min(Jeq+1, jstart+TILE_SIZE_Y-1)
    iend = min(Ieq+1, istart+TILE_SIZE_X-1)

    do j=jstart,jend ; do i=istart,iend
      ii = i-istart+1 ; jj = j-jstart+1
      dz = e(i,j,K) - e(i,j,K+1)
      do n=1,5
        p5((ii-1)*5+n,jj) = -GxRho*((e(i,j,K) - z0pres(i,j)) - 0.25*real(n-1)*dz)
        S5((ii-1)*5+n,jj) = wt_t(n) * S_t(i,j,k) + wt_b(n) * S_b(i,j,k)
        T5((ii-1)*5+n,jj) = wt_t(n) * T_t(i,j,k) + wt_b(n) * T_b(i,j,k)
      enddo
      if (use_varT)    T25((ii-1)*5+1:(ii-1)*5+5,jj) = tv%varT(i,j,k)
      if (use_covarTS) TS5((ii-1)*5+1:(ii-1)*5+5,jj) = tv%covarTS(i,j,k)
      if (use_varS)    S25((ii-1)*5+1:(ii-1)*5+5,jj) = tv%varS(i,j,k)
    enddo ; enddo

    EOSdom_h5(1,1) = 1 ; EOSdom_h5(1,2) = 5*(iend-istart+1)
    EOSdom_h5(2,1) = 1 ; EOSdom_h5(2,2) = jend-jstart+1

    if (use_stanley_eos) then
      call calculate_density(T5, S5, p5, T25, TS5, S25, r5, EOS, EOSdom_h5, rho_ref=rho_ref)
    else
      if (use_rho_ref) then
        call calculate_density(T5, S5, p5, r5, EOS, EOSdom_h5, rho_ref=rho_ref)
      else
        call calculate_density(T5, S5, p5, r5, EOS, EOSdom_h5)
        do j=jstart,jend ; do i=istart,iend
          ii = i-istart+1 ; jj = j-jstart+1
          u5((ii-1)*5+1:(ii-1)*5+5,jj) = r5((ii-1)*5+1:(ii-1)*5+5,jj) - rho_ref
        enddo ; enddo
      endif
    endif

    if (use_rho_ref) then
      do j=jstart,jend ; do i=istart,iend
        ii = i-istart+1 ; jj = j-jstart+1
        dz = e(i,j,K) - e(i,j,K+1)
        rho_anom = C1_90*(7.0*(r5((ii-1)*5+1,jj)+r5((ii-1)*5+5,jj)) + 32.0*(r5((ii-1)*5+2,jj)+r5((ii-1)*5+4,jj)) + 12.0*r5((ii-1)*5+3,jj))
        dpa(i,j) = G_e*dz*rho_anom
        if (present(intz_dpa)) then
          intz_dpa(i,j) = 0.5*G_e*dz**2 * &
                  (rho_anom - C1_90*(16.0*(r5((ii-1)*5+4,jj)-r5((ii-1)*5+2,jj)) + 7.0*(r5((ii-1)*5+5,jj)-r5((ii-1)*5+1,jj))) )
        endif
      enddo ; enddo
    else
      do j=jstart,jend ; do i=istart,iend
        ii = i-istart+1 ; jj = j-jstart+1
        dz = e(i,j,K) - e(i,j,K+1)
        rho_anom = C1_90*(7.0*(r5((ii-1)*5+1,jj)+r5((ii-1)*5+5,jj)) + 32.0*(r5((ii-1)*5+2,jj)+r5((ii-1)*5+4,jj)) + 12.0*r5((ii-1)*5+3,jj)) &
                   - rho_ref
        dpa(i,j) = G_e*dz*rho_anom
        if (present(intz_dpa)) then
          intz_dpa(i,j) = 0.5*G_e*dz**2 * &
                  (rho_anom - C1_90*(16.0*(u5((ii-1)*5+4,jj)-u5((ii-1)*5+2,jj)) + 7.0*(u5((ii-1)*5+5,jj)-u5((ii-1)*5+1,jj))) )
        endif
      enddo ; enddo
    endif

  enddo ; enddo

end subroutine generic_plm_update_dpa

subroutine generic_plm_update_intx_dpa(TILE_SIZE_X, TILE_SIZE_Y, k, tv, T_t, T_b, S_t, S_b, e, &
                                       bathyT, z0pres, dpa, intx_dpa, GxRho, G_e, dz_subroundoff, &
                                       massWeightToggle, TopWeightToggle, massWeightNVonlyToggle, h_nonvanished, &
                                       use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS, &
                                       rho_ref, EOS, HI, GV)
  integer,              intent(in)  :: TILE_SIZE_X, TILE_SIZE_Y
  integer,              intent(in)  :: k
  type(hor_index_type), intent(in)  :: HI
  type(verticalGrid_type), intent(in) :: GV
  type(thermo_var_ptrs), intent(in) :: tv
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), intent(in) :: T_t, T_b, S_t, S_b
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)+1), intent(in) :: e
  real, dimension(SZI_(HI),SZJ_(HI)), intent(in) :: bathyT
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), intent(in) :: z0pres
  real, dimension(SZI_(HI),SZJ_(HI)), intent(in) :: dpa
  real, dimension(SZIB_(HI),SZJ_(HI)), intent(inout) :: intx_dpa
  real,                 intent(in)  :: GxRho, G_e, dz_subroundoff, rho_ref
  real,                 intent(in)  :: massWeightToggle, TopWeightToggle, massWeightNVonlyToggle, h_nonvanished
  logical,              intent(in)  :: use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS
  type(EOS_type),       intent(in)  :: EOS

  real :: T15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: S15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: T215(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: TS15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: S215(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: p15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: r15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: dz_x(5,TILE_SIZE_X,TILE_SIZE_Y)
  real :: wt_t(5), wt_b(5)
  real :: intz(5)
  real, parameter :: C1_90 = 1.0/90.0
  real :: w_left, w_right, hWght, hWghtTop, iDenom, hL, hR
  real :: Ttl, Tbl, Ttr, Tbr, Stl, Sbl, Str, Sbr
  integer, dimension(2,2) :: EOSdom_q15
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n, pos, jstart, jend, istart, iend, ii, jj

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  T215(:,:) = 0.
  TS15(:,:) = 0.
  S215(:,:) = 0.

  do n = 1, 5
    wt_t(n) = 0.25 * real(5-n)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  do jstart=HI%jsc,HI%jec,TILE_SIZE_Y ; do istart=Isq,Ieq,TILE_SIZE_X
    jend = min(HI%jec, jstart+TILE_SIZE_Y-1) ; iend = min(Ieq, istart+TILE_SIZE_X-1)

    do j=jstart,jend ; do I=istart,iend
      ii = i-istart+1 ; jj = j-jstart+1
      hWght = massWeightToggle * &
              max(0., -bathyT(i,j)-e(i+1,j,K), -bathyT(i+1,j)-e(i,j,K))
      hWghtTop = TopWeightToggle * &
              max(0., e(i+1,j,K+1)-e(i,j,1), e(i,j,K+1)-e(i+1,j,1))
      hWght = max(hWght, hWghtTop)
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
        dz_x(m,ii,jj) = (w_left*(e(i,j,K) - e(i,j,K+1))) + (w_right*(e(i+1,j,K) - e(i+1,j,K+1)))

        pos = (ii-1)*15+(m-2)*5
        T15(pos+1,jj) = (w_left*Ttl) + (w_right*Ttr)
        T15(pos+5,jj) = (w_left*Tbl) + (w_right*Tbr)

        S15(pos+1,jj) = (w_left*Stl) + (w_right*Str)
        S15(pos+5,jj) = (w_left*Sbl) + (w_right*Sbr)

        p15(pos+1,jj) = -GxRho * ((w_left*(e(i,j,K)-z0pres(i,j))) + (w_right*(e(i+1,j,K)-z0pres(i+1,j))))

        do n=2,5
          p15(pos+n,jj) = p15(pos+n-1,jj) + GxRho*0.25*dz_x(m,ii,jj)
        enddo

        do n=2,4
          S15(pos+n,jj) = wt_t(n) * S15(pos+1,jj) + wt_b(n) * S15(pos+5,jj)
          T15(pos+n,jj) = wt_t(n) * T15(pos+1,jj) + wt_b(n) * T15(pos+5,jj)
        enddo
        if (use_varT) T215(pos+1:pos+5,jj) = (w_left*tv%varT(i,j,k)) + (w_right*tv%varT(i+1,j,k))
        if (use_covarTS) TS15(pos+1:pos+5,jj) = (w_left*tv%covarTS(i,j,k)) + (w_right*tv%covarTS(i+1,j,k))
        if (use_varS) S215(pos+1:pos+5,jj) = (w_left*tv%varS(i,j,k)) + (w_right*tv%varS(i+1,j,k))
      enddo
    enddo ; enddo

    EOSdom_q15(1,1) = 1 ; EOSdom_q15(1,2) = 15*(iend-istart+1)
    EOSdom_q15(2,1) = 1 ; EOSdom_q15(2,2) = jend-jstart+1

    if (use_stanley_eos) then
      call calculate_density(T15, S15, p15, T215, TS15, S215, r15, EOS, EOSdom_q15, rho_ref=rho_ref)
    else
      if (use_rho_ref) then
        call calculate_density(T15, S15, p15, r15, EOS, EOSdom_q15, rho_ref=rho_ref)
      else
        call calculate_density(T15, S15, p15, r15, EOS, EOSdom_q15)
      endif
    endif

    do j=jstart,jend ; do I=istart,iend
      ii = i-istart+1 ; jj = j-jstart+1
      intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)

      if (use_rho_ref) then
        do m = 2,4
          pos = (ii-1)*15+(m-2)*5
          intz(m) = (G_e*dz_x(m,ii,jj)*( C1_90*(7.0*(r15(pos+1,jj)+r15(pos+5,jj)) + 32.0*(r15(pos+2,jj)+r15(pos+4,jj)) + &
                            12.0*r15(pos+3,jj)) ))
        enddo
      else
        do m = 2,4
          pos = (ii-1)*15+(m-2)*5
          intz(m) = (G_e*dz_x(m,ii,jj)*( C1_90*(7.0*(r15(pos+1,jj)+r15(pos+5,jj)) + 32.0*(r15(pos+2,jj)+r15(pos+4,jj)) + &
                            12.0*r15(pos+3,jj)) - rho_ref ))
        enddo
      endif
      intx_dpa(I,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    enddo ; enddo

  enddo ; enddo

end subroutine generic_plm_update_intx_dpa

subroutine generic_plm_update_inty_dpa(TILE_SIZE_X, TILE_SIZE_Y, k, tv, T_t, T_b, S_t, S_b, e, &
                                       bathyT, z0pres, dpa, inty_dpa, GxRho, G_e, dz_subroundoff, &
                                       massWeightToggle, TopWeightToggle, massWeightNVonlyToggle, h_nonvanished, &
                                       use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS, &
                                       rho_ref, EOS, HI, GV)
  integer,              intent(in)  :: TILE_SIZE_X, TILE_SIZE_Y
  integer,              intent(in)  :: k
  type(hor_index_type), intent(in)  :: HI
  type(verticalGrid_type), intent(in) :: GV
  type(thermo_var_ptrs), intent(in) :: tv
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), intent(in) :: T_t, T_b, S_t, S_b
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)+1), intent(in) :: e
  real, dimension(SZI_(HI),SZJ_(HI)), intent(in) :: bathyT
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), intent(in) :: z0pres
  real, dimension(SZI_(HI),SZJ_(HI)), intent(in) :: dpa
  real, dimension(SZI_(HI),SZJB_(HI)), intent(inout) :: inty_dpa
  real,                 intent(in)  :: GxRho, G_e, dz_subroundoff, rho_ref
  real,                 intent(in)  :: massWeightToggle, TopWeightToggle, massWeightNVonlyToggle, h_nonvanished
  logical,              intent(in)  :: use_rho_ref, use_stanley_eos, use_varT, use_covarTS, use_varS
  type(EOS_type),       intent(in)  :: EOS

  real :: T15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: S15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: T215(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: TS15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: S215(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: p15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: r15(15*TILE_SIZE_X,TILE_SIZE_Y)
  real :: dz_y(5,TILE_SIZE_X,TILE_SIZE_Y)
  real :: wt_t(5), wt_b(5)
  real :: intz(5)
  real, parameter :: C1_90 = 1.0/90.0
  real :: w_left, w_right, hWght, hWghtTop, iDenom, hL, hR
  real :: Ttl, Tbl, Ttr, Tbr, Stl, Sbl, Str, Sbr
  integer, dimension(2,2) :: EOSdom_h15
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n, pos, jstart, jend, istart, iend, ii, jj

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  T215(:,:) = 0.
  TS15(:,:) = 0.
  S215(:,:) = 0.

  do n = 1, 5
    wt_t(n) = 0.25 * real(5-n)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  do jstart=Jsq,Jeq,TILE_SIZE_Y ; do istart=HI%isc,HI%iec,TILE_SIZE_X
    jend = min(Jeq, jstart+TILE_SIZE_Y-1) ; iend = min(HI%iec, istart+TILE_SIZE_X-1)

    do j=jstart,jend ; do i=istart,iend
      ii = i-istart+1 ; jj = j-jstart+1
      hWght = massWeightToggle * &
              max(0., -bathyT(i,j)-e(i,j+1,K), -bathyT(i,j+1)-e(i,j,K))
      hWghtTop = TopWeightToggle * &
              max(0., e(i,j+1,K+1)-e(i,j,1), e(i,j,K+1)-e(i,j+1,1))
      hWght = max(hWght, hWghtTop)
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
        dz_y(m,ii,jj) = (w_left*(e(i,j,K) - e(i,j,K+1))) + (w_right*(e(i,j+1,K) - e(i,j+1,K+1)))

        pos = (ii-1)*15+(m-2)*5
        T15(pos+1,jj) = (w_left*Ttl) + (w_right*Ttr)
        T15(pos+5,jj) = (w_left*Tbl) + (w_right*Tbr)

        S15(pos+1,jj) = (w_left*Stl) + (w_right*Str)
        S15(pos+5,jj) = (w_left*Sbl) + (w_right*Sbr)

        p15(pos+1,jj) = -GxRho * ((w_left*(e(i,j,K)-z0pres(i,j))) + (w_right*(e(i,j+1,K)-z0pres(i,j+1))))

        do n=2,5
          p15(pos+n,jj) = p15(pos+n-1,jj) + GxRho*0.25*dz_y(m,ii,jj)
        enddo

        do n=2,4
          S15(pos+n,jj) = wt_t(n) * S15(pos+1,jj) + wt_b(n) * S15(pos+5,jj)
          T15(pos+n,jj) = wt_t(n) * T15(pos+1,jj) + wt_b(n) * T15(pos+5,jj)
        enddo
        if (use_varT) T215(pos+1:pos+5,jj) = (w_left*tv%varT(i,j,k)) + (w_right*tv%varT(i,j+1,k))
        if (use_covarTS) TS15(pos+1:pos+5,jj) = (w_left*tv%covarTS(i,j,k)) + (w_right*tv%covarTS(i,j+1,k))
        if (use_varS) S215(pos+1:pos+5,jj) = (w_left*tv%varS(i,j,k)) + (w_right*tv%varS(i,j+1,k))
      enddo
    enddo ; enddo

    EOSdom_h15(1,1) = 1 ; EOSdom_h15(1,2) = 15*(iend-istart+1)
    EOSdom_h15(2,1) = 1 ; EOSdom_h15(2,2) = jend-jstart+1

    if (use_stanley_eos) then
      call calculate_density(T15, S15, p15, T215, TS15, S215, r15, EOS, EOSdom_h15, rho_ref=rho_ref)
    else
      if (use_rho_ref) then
        call calculate_density(T15, S15, p15, r15, EOS, EOSdom_h15, rho_ref=rho_ref)
      else
        call calculate_density(T15, S15, p15, r15, EOS, EOSdom_h15)
      endif
    endif

    do j=jstart,jend ; do i=istart,iend
      ii = i-istart+1 ; jj = j-jstart+1
      intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)

      if (use_rho_ref) then
        do m = 2,4
          pos = (ii-1)*15+(m-2)*5
          intz(m) = (G_e*dz_y(m,ii,jj)*( C1_90*(7.0*(r15(pos+1,jj)+r15(pos+5,jj)) + &
                                           32.0*(r15(pos+2,jj)+r15(pos+4,jj)) + &
                                           12.0*r15(pos+3,jj)) ))
        enddo
      else
        do m = 2,4
          pos = (ii-1)*15+(m-2)*5
          intz(m) = (G_e*dz_y(m,ii,jj)*( C1_90*(7.0*(r15(pos+1,jj)+r15(pos+5,jj)) + &
                                           32.0*(r15(pos+2,jj)+r15(pos+4,jj)) + &
                                           12.0*r15(pos+3,jj)) - rho_ref ))
        enddo
      endif
      inty_dpa(i,J) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    enddo ; enddo

  enddo ; enddo

end subroutine generic_plm_update_inty_dpa

end submodule MOM_density_integrals_s

!> \namespace mom_density_integrals
!!
