! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Isopycnal height diffusion (or Gent McWilliams diffusion)
submodule (MOM_thickness_diffuse) MOM_thickness_diffuse_s
#include <MOM_memory.h>
contains

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

  ! Local variables
  real :: e(SZI_(G),SZJ_(G),SZK_(GV)+1) ! heights of interfaces, relative to mean
                                         ! sea level [Z ~> m], positive up.
  real :: uhD(SZIB_(G),SZJ_(G),SZK_(GV)) ! Diffusive u*h fluxes [L2 H T-1 ~> m3 s-1 or kg s-1]
  real :: vhD(SZI_(G),SZJB_(G),SZK_(GV)) ! Diffusive v*h fluxes [L2 H T-1 ~> m3 s-1 or kg s-1]

  real :: Sfn_unlim_u_3D(SZIB_(G), SZJ_(G),SZK_(GV)+1) ! Volume streamfunction for u-points [Z L2 T-1 ~> m3 s-1]
  real :: Sfn_unlim_v_3D(SZI_(G), SZJB_(G),SZK_(GV)+1)  ! Volume streamfunction for v-points [Z L2 T-1 ~> m3 s-1]

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: &
    KH_u, &       ! Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
    int_slope_u   ! A nondimensional ratio from 0 to 1 that gives the relative
                  ! weighting of the interface slopes to that calculated also
                  ! using density gradients at u points.  The physically correct
                  ! slopes occur at 0, while 1 is used for numerical closures [nondim].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: &
    KH_v, &       ! Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
    int_slope_v   ! A nondimensional ratio from 0 to 1 that gives the relative
                  ! weighting of the interface slopes to that calculated also
                  ! using density gradients at v points.  The physically correct
                  ! slopes occur at 0, while 1 is used for numerical closures [nondim].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    KH_t          ! diagnosed diffusivity at tracer points [L2 T-1 ~> m2 s-1]

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    KH_u_CFL      ! The maximum stable isopycnal height diffusivity at u grid points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    KH_v_CFL      ! The maximum stable isopycnal height diffusivity at v grid points [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot          ! The sum of the total layer thicknesses [H ~> m or kg m-2]
  real :: Khth_Loc_u(SZIB_(G),SZJ_(G)) ! The isopycnal height diffusivity at u points [L2 T-1 ~> m2 s-1]
  real :: Khth_Loc_v(SZI_(G),SZJB_(G)) ! The isopycnal height diffusivity at v points [L2 T-1 ~> m2 s-1]
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real, dimension(:,:), pointer :: cg1 => null() !< Wave speed [L T-1 ~> m s-1]
  real :: hu(SZI_(G),SZJ_(G))       ! A thickness-based mask at u points, used for diagnostics [nondim]
  real :: hv(SZI_(G),SZJ_(G))       ! A thickness-based mask at v points, used for diagnostics [nondim]
  real :: KH_u_lay(SZI_(G),SZJ_(G)) ! Diagnostic of isopycnal height diffusivities at u-points averaged
                                    ! to layer centers [L2 T-1 ~> m2 s-1]
  real :: KH_v_lay(SZI_(G),SZJ_(G)) ! Diagnostic of isopycnal height diffusivities at v-points averaged
                                    ! to layer centers [L2 T-1 ~> m2 s-1]
  logical :: use_VarMix, Resoln_scaled, Depth_scaled, use_stored_slopes, khth_use_vert_struct, use_Visbeck
  logical :: use_QG_Leith
  integer :: i, j, k, is, ie, js, je, nz

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_thickness_diffuse: "//&
         "Module must be initialized before it is used.")

  if ((.not.CS%thickness_diffuse) &
      .or. .not. (CS%Khth > 0.0 .or. CS%read_khth &
      .or. VarMix%use_variable_mixing)) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  h_neglect = GV%H_subroundoff

  if (allocated(MEKE%GM_src)) then
    do j=js,je ; do i=is,ie ; MEKE%GM_src(i,j) = 0. ; enddo ; enddo
  endif

  use_VarMix = .false. ; Resoln_scaled = .false. ; use_stored_slopes = .false.
  khth_use_vert_struct = .false. ; use_Visbeck = .false. ; use_QG_Leith = .false.
  Depth_scaled = .false.

  if (VarMix%use_variable_mixing) then
    use_VarMix = VarMix%use_variable_mixing .and. (CS%KHTH_Slope_Cff > 0.)
    Resoln_scaled = VarMix%Resoln_scaled_KhTh
    Depth_scaled = VarMix%Depth_scaled_KhTh
    use_stored_slopes = VarMix%use_stored_slopes
    khth_use_vert_struct = allocated(VarMix%khth_struct)
    use_Visbeck = VarMix%use_Visbeck
    use_QG_Leith = VarMix%use_QG_Leith_GM
    if (allocated(VarMix%cg1)) cg1 => VarMix%cg1
  else
    cg1 => null()
  endif

  if (is_root_pe()) then
    write(*,'(A)') "[thickness_diffuse] Branch flags:"
    write(*,'(2X,A,L1)') "use_VarMix           = ", use_VarMix
    write(*,'(2X,A,L1)') "Resoln_scaled        = ", Resoln_scaled
    write(*,'(2X,A,L1)') "Depth_scaled         = ", Depth_scaled
    write(*,'(2X,A,L1)') "use_stored_slopes    = ", use_stored_slopes
    write(*,'(2X,A,L1)') "khth_use_vert_struct = ", khth_use_vert_struct
    write(*,'(2X,A,L1)') "use_Visbeck          = ", use_Visbeck
    write(*,'(2X,A,L1)') "use_QG_Leith         = ", use_QG_Leith
    write(*,'(2X,A,L1)') "CS%read_khth         = ", CS%read_khth
    write(*,'(2X,A,L1)') "MEKE%Kh allocated    = ", allocated(MEKE%Kh)
    if (allocated(MEKE%Kh)) &
      write(*,'(2X,A,L1)') "CS%MEKE_GEOMETRIC    = ", CS%MEKE_GEOMETRIC
    write(*,'(2X,A,L1)') "full_depth_khth_min  = ", CS%full_depth_khth_min
    write(*,'(2X,A,L1)') "use_GME_thickness    = ", CS%use_GME_thickness_diffuse
    write(*,'(2X,A,L1)') "Khth_Max > 0         = ", (CS%Khth_Max > 0.0)
    write(*,'(2X,A,L1)') "max_Khth_CFL > 0     = ", (CS%max_Khth_CFL > 0.0)
  endif

  !$omp target update to(CS, MEKE)
  !$omp target enter data map(alloc: KH_u_CFL, KH_v_CFL, Khth_Loc_u, Khth_Loc_v, int_slope_u, int_slope_v, e, KH_u, KH_v)

  do concurrent (j=js:je, I=is-1:ie)
    KH_u_CFL(I,j) = (0.25*CS%max_Khth_CFL) /  &
      (dt * ((G%IdxCu(I,j)*G%IdxCu(I,j)) + (G%IdyCu(I,j)*G%IdyCu(I,j))))
  enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-1,je ; do i=is,ie
    KH_v_CFL(i,J) = (0.25*CS%max_Khth_CFL) / &
      (dt * ((G%IdxCv(i,J)*G%IdxCv(i,J)) + (G%IdyCv(i,J)*G%IdyCv(i,J))))
  enddo ; enddo

  ! Calculates interface heights, e, in [Z ~> m].
  if (CS%use_meso_sfn_ANN) then
    ! The ANN streamfunction needs a wider halo on e.
    call find_eta(h, tv, G, GV, US, e, halo_size=3)
  else
    call find_eta(h, tv, G, GV, US, e, halo_size=1)
  endif

  ! Set the diffusivities.
  if (.not. CS%read_khth) then
    do concurrent (j=js:je, I=is-1:ie)
      Khth_loc_u(I,j) = CS%Khth
    enddo
  else ! use 2d KHTH that was read in from file
    do concurrent (j=js:je, I=is-1:ie)
      Khth_loc_u(I,j) = 0.5 * (CS%khth2d(i,j) + CS%khth2d(i+1,j))
    enddo
  endif

  if (use_VarMix) then
    if (use_Visbeck) then
      !$omp target update from( VarMix%L2u, VarMix%SN_u)
      do concurrent (j=js:je, I=is-1:ie)
        Khth_loc_u(I,j) = Khth_loc_u(I,j) + &
          CS%KHTH_Slope_Cff*VarMix%L2u(I,j) * VarMix%SN_u(I,j)
      enddo
    endif
  endif

  if (allocated(MEKE%Kh)) then
    if (CS%MEKE_GEOMETRIC) then
      !$omp target update from( VarMix%SN_u)
      do concurrent (j=js:je, I=is-1:ie)
        Khth_loc_u(I,j) = Khth_loc_u(I,j) + G%OBCmaskCu(I,j) * CS%MEKE_GEOMETRIC_alpha * &
                          0.5*(MEKE%MEKE(i,j)+MEKE%MEKE(i+1,j)) / &
                          (VarMix%SN_u(I,j) + CS%MEKE_GEOMETRIC_epsilon)
      enddo
    else
      do concurrent (j=js:je, I=is-1:ie)
        Khth_loc_u(I,j) = Khth_loc_u(I,j) + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i+1,j))
      enddo
    endif
  endif

  if (Resoln_scaled) then
    !$omp target update from( VarMix%Res_fn_u )
    do concurrent (j=js:je, I=is-1:ie)
      Khth_loc_u(I,j) = Khth_loc_u(I,j) * VarMix%Res_fn_u(I,j)
    enddo
  endif

  if (Depth_scaled) then
    !$omp target update from( VarMix%Depth_fn_u )
    do concurrent (j=js:je, I=is-1:ie)
      Khth_loc_u(I,j) = Khth_loc_u(I,j) * VarMix%Depth_fn_u(I,j)
    enddo
  endif

  if (CS%Khth_Max > 0) then
    do concurrent (j=js:je, I=is-1:ie)
      Khth_loc_u(I,j) = max(CS%Khth_Min, min(Khth_loc_u(I,j), CS%Khth_Max))
    enddo
  else
    do concurrent (j=js:je, I=is-1:ie)
      Khth_loc_u(I,j) = max(CS%Khth_Min, Khth_loc_u(I,j))
    enddo
  endif
  do concurrent(j=js:je, I=is-1:ie)
    KH_u(I,j,1) = min(KH_u_CFL(I,j), Khth_loc_u(I,j))
  enddo

  if (khth_use_vert_struct) then
    if (CS%full_depth_khth_min) then
      do concurrent (K=2:nz+1, j=js:je, I=is-1:ie)
        KH_u(I,j,K) = KH_u(I,j,1) * 0.5 * ( VarMix%khth_struct(i,j,k-1) + VarMix%khth_struct(i+1,j,k-1) )
        KH_u(I,j,K) = max(KH_u(I,j,K), CS%Khth_Min)
      enddo
    else
      do concurrent (K=2:nz+1, j=js:je, I=is-1:ie)
        KH_u(I,j,K) = KH_u(I,j,1) * 0.5 * ( VarMix%khth_struct(i,j,k-1) + VarMix%khth_struct(i+1,j,k-1) )
      enddo
    endif
  else
    do concurrent (K=2:nz+1, j=js:je, I=is-1:ie)
      KH_u(I,j,K) = KH_u(I,j,1)
    enddo
  endif

  if (use_VarMix) then
    if (use_QG_Leith) then
      do concurrent (k=1:nz, j=js:je, I=is-1:ie)
        KH_u(I,j,k) = VarMix%KH_u_QG(I,j,k)
      enddo
    endif
  endif

  if (CS%use_GME_thickness_diffuse) then
    do concurrent (k=1:nz+1, j=js:je, I=is-1:ie)
      CS%KH_u_GME(I,j,k) = KH_u(I,j,k)
    enddo
  endif

  if (.not. CS%read_khth) then
    do concurrent (J=js-1:je, i=is:ie)
      Khth_loc_v(i,J) = CS%Khth
    enddo
  else ! read KHTH from file
    do concurrent (J=js-1:je, i=is:ie)
      Khth_loc_v(i,J) = 0.5 * (CS%khth2d(i,j) + CS%khth2d(i,j+1))
    enddo
  endif

  if (use_VarMix) then
    if (use_Visbeck) then
      !$omp target update from( VarMix%L2v, VarMix%SN_v )
      do concurrent (J=js-1:je, i=is:ie)
        Khth_loc_v(i,J) = Khth_loc_v(i,J) + CS%KHTH_Slope_Cff*VarMix%L2v(i,J)*VarMix%SN_v(i,J)
      enddo
    endif
  endif
  if (allocated(MEKE%Kh)) then
    if (CS%MEKE_GEOMETRIC) then
      !$omp target update from( VarMix%SN_v )
      do concurrent (J=js-1:je, i=is:ie)
        Khth_loc_v(i,J) = Khth_loc_v(i,J) + G%OBCmaskCv(i,J) * CS%MEKE_GEOMETRIC_alpha * &
                        0.5*(MEKE%MEKE(i,j)+MEKE%MEKE(i,j+1)) / &
                        (VarMix%SN_v(i,J) + CS%MEKE_GEOMETRIC_epsilon)
      enddo
    else
      do concurrent (J=js-1:je, i=is:ie)
        Khth_loc_v(i,J) = Khth_loc_v(i,J) + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i,j+1))
      enddo
    endif
  endif

  if (Resoln_scaled) then
    !$omp target update from( VarMix%Res_fn_v )
    do concurrent (J=js-1:je, i=is:ie)
      Khth_loc_v(i,J) = Khth_loc_v(i,J) * VarMix%Res_fn_v(i,J)
    enddo
  endif

  if (Depth_scaled) then
    !$omp target update from( VarMix%Depth_fn_v )
    do concurrent (J=js-1:je, i=is:ie)
      Khth_loc_v(i,J) = Khth_loc_v(i,J) * VarMix%Depth_fn_v(i,J)
    enddo
  endif

  if (CS%Khth_Max > 0) then
    do concurrent (J=js-1:je, i=is:ie)
      Khth_loc_v(i,J) = max(CS%Khth_Min, min(Khth_loc_v(i,J), CS%Khth_Max))
    enddo
  else
    do concurrent (J=js-1:je, i=is:ie)
      Khth_loc_v(i,J) = max(CS%Khth_Min, Khth_loc_v(i,J))
    enddo
  endif

  if (CS%max_Khth_CFL > 0.0) then
    do concurrent (J=js-1:je, i=is:ie)
      KH_v(i,J,1) = min(KH_v_CFL(i,J), Khth_loc_v(i,J))
    enddo
  endif

  if (khth_use_vert_struct) then
      if (CS%full_depth_khth_min) then
      do concurrent (K=2:nz+1, J=js-1:je, i=is:ie)
        KH_v(i,J,K) = KH_v(i,J,1) * 0.5 * ( VarMix%khth_struct(i,j,k-1) + VarMix%khth_struct(i,j+1,k-1) )
        KH_v(i,J,K) = max(KH_v(i,J,K), CS%Khth_Min)
      enddo
    else
      do concurrent (K=2:nz+1, J=js-1:je, i=is:ie)
        KH_v(i,J,K) = KH_v(i,J,1) * 0.5 * ( VarMix%khth_struct(i,j,k-1) + VarMix%khth_struct(i,j+1,k-1) )
      enddo
    endif
  else
    do concurrent (K=2:nz+1, J=js-1:je, i=is:ie)
      KH_v(i,J,K) = KH_v(i,J,1)
    enddo
  endif

  if (use_VarMix) then
    if (use_QG_Leith) then
      do concurrent (k=1:nz, J=js-1:je, i=is:ie)
        KH_v(i,J,k) = VarMix%KH_v_QG(i,J,k)
      enddo
    endif
  endif

  if (CS%use_GME_thickness_diffuse) then
    do concurrent (k=1:nz+1, J=js-1:je, i=is:ie)
      CS%KH_v_GME(i,J,k) = KH_v(i,J,k)
    enddo
  endif

  if (allocated(MEKE%Kh)) then
    if (CS%MEKE_GEOMETRIC) then
      !$omp target update from( VarMix%SN_u, VarMix%SN_v )
      if (CS%MEKE_GEOM_answer_date < 20190101) then
        do concurrent (j=js:je, i=is:ie)
          ! This does not give bitwise rotational symmetry.
          MEKE%Kh(i,j) = CS%MEKE_GEOMETRIC_alpha * MEKE%MEKE(i,j) / &
                         (0.25*(VarMix%SN_u(I,j)+VarMix%SN_u(I-1,j) + &
                                VarMix%SN_v(i,J)+VarMix%SN_v(i,J-1)) + &
                          CS%MEKE_GEOMETRIC_epsilon)
        enddo
      else
        do concurrent (j=js:je, i=is:ie)
          ! With the additional parentheses this gives bitwise rotational symmetry.
          MEKE%Kh(i,j) = CS%MEKE_GEOMETRIC_alpha * MEKE%MEKE(i,j) / &
                         (0.25*((VarMix%SN_u(I,j)+VarMix%SN_u(I-1,j)) + &
                                (VarMix%SN_v(i,J)+VarMix%SN_v(i,J-1))) + &
                          CS%MEKE_GEOMETRIC_epsilon)
        enddo
      endif
    endif
  endif

  do concurrent (K=1:nz+1, j=js:je, I=is-1:ie) ; int_slope_u(I,j,K) = 0.0 ; enddo
  do concurrent (K=1:nz+1, J=js-1:je, i=is:ie) ; int_slope_v(i,J,K) = 0.0 ; enddo

  !$omp target update from(KH_u_CFL, KH_v_CFL, Khth_Loc_u, Khth_Loc_v, int_slope_u, int_slope_v, e, KH_u, KH_v)

  if (is_root_pe()) write(*,'(A)') "[thickness_diffuse] Post-parallel branches:"
  if (CS%detangle_interfaces) then
    call add_detangling_Kh(h, e, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, tv, dt, G, GV, US, &
                           CS, int_slope_u, int_slope_v)
  endif

  if ((CS%Kh_eta_bg > 0.0) .or. (CS%Kh_eta_vel > 0.0)) then
    call add_interface_Kh(G, GV, US, CS, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, int_slope_u, int_slope_v)
  endif

  if (CS%use_meso_sfn_ANN) then
    call meso_sfn_ANN_compute(h, e, Sfn_unlim_u_3D, Sfn_unlim_v_3D, G, GV, US, tv, &
                              CS%meso_sfn_ANN_CS, dt, u, v)
  endif

  if (CS%debug) then
    call uvchksum("Kh_[uv]", Kh_u, Kh_v, G%HI, haloshift=0, &
                  unscale=(US%L_to_m**2)*US%s_to_T, scalar_pair=.true.)
    call uvchksum("Kh_[uv]_CFL", Kh_u_CFL, Kh_v_CFL, G%HI, haloshift=0, &
                  unscale=(US%L_to_m**2)*US%s_to_T, scalar_pair=.true.)
    if (Resoln_scaled) then
      call uvchksum("Res_fn_[uv]", VarMix%Res_fn_u, VarMix%Res_fn_v, G%HI, haloshift=0, &
                    unscale=1.0, scalar_pair=.true.)
    endif
    call uvchksum("int_slope_[uv]", int_slope_u, int_slope_v, G%HI, haloshift=0)
    call hchksum(h, "thickness_diffuse_1 h", G%HI, haloshift=1, unscale=GV%H_to_m)
    call hchksum(e, "thickness_diffuse_1 e", G%HI, haloshift=1, unscale=US%Z_to_m)
    if (use_stored_slopes) then
      !$omp target update from(VarMix%slope_x, VarMix%slope_y)
      call uvchksum("VarMix%slope_[xy]", VarMix%slope_x, VarMix%slope_y, &
                    G%HI, haloshift=0, unscale=US%Z_to_L)
    endif
    if (associated(tv%eqn_of_state)) then
      call hchksum(tv%T, "thickness_diffuse T", G%HI, haloshift=1, unscale=US%C_to_degC)
      call hchksum(tv%S, "thickness_diffuse S", G%HI, haloshift=1, unscale=US%S_to_ppt)
    endif
  endif

  ! Calculate uhD, vhD from h, e, KH_u, KH_v, tv%T/S
  if (STOCH%skeb_use_gm) then
    if (use_stored_slopes) then
      !$omp target update from(VarMix%slope_x, VarMix%slope_y)
      call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, CS, &
                                  int_slope_u, int_slope_v, VarMix%slope_x, VarMix%slope_y, &
                                  STOCH=STOCH, VarMix=VarMix, &
                                  Sfn_unlim_u_3D=Sfn_unlim_u_3D, Sfn_unlim_v_3D=Sfn_unlim_v_3D)
    else
      call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, CS, &
                                  int_slope_u, int_slope_v, STOCH=STOCH, VarMix=VarMix, &
                                  Sfn_unlim_u_3D=Sfn_unlim_u_3D, Sfn_unlim_v_3D=Sfn_unlim_v_3D)
    endif
  else
    if (use_stored_slopes) then
      !$omp target update from(VarMix%slope_x, VarMix%slope_y)
      call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, CS, &
                                  int_slope_u, int_slope_v, VarMix%slope_x, VarMix%slope_y, &
                                  Sfn_unlim_u_3D=Sfn_unlim_u_3D, Sfn_unlim_v_3D=Sfn_unlim_v_3D)
    else
      call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, cg1, dt, G, GV, US, MEKE, CS, &
                                  int_slope_u, int_slope_v, &
                                  Sfn_unlim_u_3D=Sfn_unlim_u_3D, Sfn_unlim_v_3D=Sfn_unlim_v_3D)
    endif
  endif

  if (VarMix%use_variable_mixing) then
    if (allocated(MEKE%Rd_dx_h) .and. allocated(VarMix%Rd_dx_h)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%Rd_dx_h(i,j) = VarMix%Rd_dx_h(i,j)
      enddo ; enddo
    endif
  endif

  ! offer diagnostic fields for averaging
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_uhGM > 0)   call post_data(CS%id_uhGM, uhD, CS%diag)
    if (CS%id_vhGM > 0)   call post_data(CS%id_vhGM, vhD, CS%diag)
    if (CS%id_GMwork > 0) call post_data(CS%id_GMwork, CS%GMwork, CS%diag)
    if (CS%id_KH_u > 0)   call post_data(CS%id_KH_u, KH_u, CS%diag)
    if (CS%id_KH_v > 0)   call post_data(CS%id_KH_v, KH_v, CS%diag)
    if (CS%id_KH_u1 > 0)  call post_data(CS%id_KH_u1, KH_u(:,:,1), CS%diag)
    if (CS%id_KH_v1 > 0)  call post_data(CS%id_KH_v1, KH_v(:,:,1), CS%diag)

    ! Diagnose diffusivity at T-cell point.  Do a simple average, rather than a
    ! thickness-weighted average, so that KH_t is depth-independent when KH_u and KH_v
    ! are depth independent.  If a thickness-weighted average were used, the variations
    ! of thickness could give a spurious depth dependence to the diagnosed KH_t.
    if (CS%id_KH_t > 0 .or. CS%id_KH_t1 > 0 .or. CS%Use_KH_in_MEKE) then
      do k=1,nz
        ! thicknesses across u and v faces, converted to 0/1 mask
        ! layer average of the interface diffusivities KH_u and KH_v
        do j=js,je ; do I=is-1,ie
          ! This expression uses harmonic mean thicknesses:
          ! hu(I,j)       = 2.0*h(i,j,k)*h(i+1,j,k) / (h(i,j,k)+h(i+1,j,k)+h_neglect)
          ! This expression is a 0/1 mask based on depths where there are thick layers:
          hu(I,j) = 0.0 ; if (h(i,j,k)*h(i+1,j,k) /= 0.0) hu(I,j) = 1.0
          KH_u_lay(I,j) = 0.5*(KH_u(I,j,k)+KH_u(I,j,k+1))
        enddo ; enddo
        do J=js-1,je ; do i=is,ie
          ! This expression uses harmonic mean thicknesses:
          ! hv(i,J)       = 2.0*h(i,j,k)*h(i,j+1,k)/(h(i,j,k)+h(i,j+1,k)+h_neglect)
          ! This expression is a 0/1 mask based on depths where there are thick layers:
          hv(i,J) = 0.0 ; if (h(i,j,k)*h(i,j+1,k) /= 0.0) hv(i,J) = 1.0
          KH_v_lay(i,J) = 0.5*(KH_v(i,J,k)+KH_v(i,J,k+1))
        enddo ; enddo
        ! diagnose diffusivity at T-points
        do j=js,je ; do i=is,ie
          Kh_t(i,j,k) = (((hu(I-1,j)*KH_u_lay(i-1,j)) + (hu(I,j)*KH_u_lay(I,j))) + &
                         ((hv(i,J-1)*KH_v_lay(i,J-1)) + (hv(i,J)*KH_v_lay(i,J)))) / &
                        ((hu(I-1,j)+hu(I,j)) + (hv(i,J-1)+hv(i,J)) + 1.0e-20)
          ! Use this denominator instead if hu and hv are actual thicknesses rather than a 0/1 mask:
          !              ((hu(I-1,j)+hu(I,j)) + (hv(i,J-1)+hv(i,J)) + h_neglect)
        enddo ; enddo
      enddo

      if (CS%Use_KH_in_MEKE) then
        MEKE%Kh_diff(:,:) = 0.0
        htot(:,:) = 0.0
        do k=1,nz
          do j=js,je ; do i=is,ie
            MEKE%Kh_diff(i,j) = MEKE%Kh_diff(i,j) + Kh_t(i,j,k) * h(i,j,k)
            htot(i,j) = htot(i,j) + h(i,j,k)
          enddo ; enddo
        enddo

        do j=js,je ; do i=is,ie
          MEKE%Kh_diff(i,j) = MEKE%Kh_diff(i,j) / MAX(CS%MEKE_min_depth_diff, htot(i,j))
        enddo ; enddo
      endif

      if (CS%id_KH_t  > 0) call post_data(CS%id_KH_t,  KH_t,        CS%diag)
      if (CS%id_KH_t1 > 0) call post_data(CS%id_KH_t1, KH_t(:,:,1), CS%diag)
    endif

  endif

  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      uhtr(I,j,k) = uhtr(I,j,k) + uhD(I,j,k) * dt
      if (associated(CDp%uhGM)) CDp%uhGM(I,j,k) = uhD(I,j,k)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      vhtr(i,J,k) = vhtr(i,J,k) + vhD(i,J,k) * dt
      if (associated(CDp%vhGM)) CDp%vhGM(i,J,k) = vhD(i,J,k)
    enddo ; enddo
    do j=js,je ; do i=is,ie
      h(i,j,k) = h(i,j,k) - dt * G%IareaT(i,j) * &
          ((uhD(I,j,k) - uhD(I-1,j,k)) + (vhD(i,J,k) - vhD(i,J-1,k)))
      if (h(i,j,k) < GV%Angstrom_H) h(i,j,k) = GV%Angstrom_H
    enddo ; enddo
  enddo

  !$omp target exit data map(release: KH_u_CFL, KH_v_CFL, Khth_Loc_u, Khth_Loc_v, int_slope_u, int_slope_v, e, KH_u, KH_v)

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  ! This needs to happen after the H update and before the next post_data.
  call diag_update_remap_grids(CS%diag)

  if (CS%debug) then
    call uvchksum("thickness_diffuse [uv]hD", uhD, vhD, &
                  G%HI, haloshift=0, unscale=GV%H_to_m*US%L_to_m**2*US%s_to_T)
    call uvchksum("thickness_diffuse [uv]htr", uhtr, vhtr, &
                  G%HI, haloshift=0, unscale=US%L_to_m**2*GV%H_to_m)
    call hchksum(h, "thickness_diffuse h", G%HI, haloshift=0, unscale=GV%H_to_m)
  endif

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

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    T, &          ! The temperature [C ~> degC], with the values in
                  ! in massless layers filled vertically by diffusion.
    S, &          ! The filled salinity [S ~> ppt], with the values in
                  ! in massless layers filled vertically by diffusion.
    h_avail, &    ! The mass available for diffusion out of each face, divided
                  ! by dt [H L2 T-1 ~> m3 s-1 or kg s-1].
    h_frac        ! The fraction of the mass in the column above the bottom
                  ! interface of a layer that is within a layer [nondim]. 0<h_frac<=1
  real :: dz(SZI_(G),SZJ_(G),SZK_(GV)) ! Height change across layers [Z ~> m]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: &
    Slope_y_PE, &  ! 3D array of neutral slopes at v-points, set equal to Slope (below) [Z L-1 ~> nondim]
    hN2_y_PE       ! Harmonic mean of thicknesses around the interfaces times the buoyancy frequency
                   ! at v-points with unit conversion factors [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2],
                   ! used for calculating the potential energy release
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: &
    Slope_x_PE, &  ! 3D array of neutral slopes at u-points, set equal to Slope (below) [Z L-1 ~> nondim]
    hN2_x_PE       ! Harmonic mean of thicknesses around the interfaces times the buoyancy frequency
                   ! at u-points  with unit conversion factors [H L2 Z-2 T-2 ~> m s-2 or kg m-2 s-2],
                   ! used for calculating the potential energy release
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    pres, &       ! The pressure at an interface [R L2 T-2 ~> Pa].
    h_avail_rsum  ! The running sum of h_avail above an interface [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G)) :: &
    drho_dT_u, &  ! The derivative of density with temperature at u points [R C-1 ~> kg m-3 degC-1]
    drho_dS_u     ! The derivative of density with salinity at u points [R S-1 ~> kg m-3 ppt-1].
  real, dimension(SZIB_(G)) :: scrap ! An array to pass to calculate_density_second_derivs()
                  ! with various units that will be ignored [various]
  real, dimension(SZI_(G)) :: &
    drho_dT_v, &  ! The derivative of density with temperature at v points [R C-1 ~> kg m-3 degC-1]
    drho_dS_v, &  ! The derivative of density with salinity at v points [R S-1 ~> kg m-3 ppt-1].
    drho_dT_dT_h, & ! The second derivative of density with temperature at h points [R C-2 ~> kg m-3 degC-2]
    drho_dT_dT_hr ! The second derivative of density with temperature at h (+1) points [R C-2 ~> kg m-3 degC-2]
  real :: uhtot(SZIB_(G),SZJ_(G))  ! The vertical sum of uhD [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: vhtot(SZI_(G),SZJB_(G))  ! The vertical sum of vhD [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G)) :: &
    T_u, &        ! Temperature on the interface at the u-point [C ~> degC].
    S_u, &        ! Salinity on the interface at the u-point [S ~> ppt].
    pres_u        ! Pressure on the interface at the u-point [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G)) :: &
    T_v, &        ! Temperature on the interface at the v-point [C ~> degC].
    S_v, &        ! Salinity on the interface at the v-point [S ~> ppt].
    pres_v, &     ! Pressure on the interface at the v-point [R L2 T-2 ~> Pa].
    T_h, &        ! Temperature on the interface at the h-point [C ~> degC].
    S_h, &        ! Salinity on the interface at the h-point [S ~> ppt].
    pres_h, &     ! Pressure on the interface at the h-point [R L2 T-2 ~> Pa].
    T_hr, &       ! Temperature on the interface at the h (+1) point [C ~> degC].
    S_hr, &       ! Salinity on the interface at the h (+1) point [S ~> ppt].
    pres_hr       ! Pressure on the interface at the h (+1) point [R L2 T-2 ~> Pa].
  real :: Work_u(SZIB_(G),SZJ_(G)) ! The work done by the isopycnal height diffusion
                                   ! integrated over u-point water columns [R Z L4 T-3 ~> W]
  real :: Work_v(SZI_(G),SZJB_(G)) ! The work done by the isopycnal height diffusion
                                   ! integrated over v-point water columns [R Z L4 T-3 ~> W]
  real :: Work_h        ! The work averaged over an h-cell [R Z L2 T-3 ~> W m-2].
  real :: PE_release_h  ! The amount of potential energy released by GM averaged over an h-cell
                        ! [R Z L2 T-3 ~> W m-2].  The calculation equals rho0 * h * S^2 * N^2 * kappa_GM.
  real :: I4dt          ! 1 / 4 dt [T-1 ~> s-1].
  real :: drdiA, drdiB  ! Along layer zonal potential density  gradients in the layers above (A)
                        ! and below (B) the interface times the grid spacing [R ~> kg m-3].
  real :: drdjA, drdjB  ! Along layer meridional potential density  gradients in the layers above (A)
                        ! and below (B) the interface times the grid spacing [R ~> kg m-3].
  real :: drdkL, drdkR  ! Vertical density differences across an interface [R ~> kg m-3].
  real :: drdi_u(SZIB_(G),SZK_(GV)) ! Copy of drdi at u-points [R ~> kg m-3].
  real :: drdj_v(SZI_(G),SZK_(GV)) ! Copy of drdj at v-points [R ~> kg m-3].
  real :: drdkDe_u(SZIB_(G),SZK_(GV)+1) ! Lateral difference of product of drdk and e at u-points
                                        ! [Z R ~> kg m-2].
  real :: drdkDe_v(SZI_(G),SZK_(GV)+1)  ! Lateral difference of product of drdk and e at v-points
                                        ! [Z R ~> kg m-2].
  real :: hg2A, hg2B, hg2L, hg2R ! Squares of geometric mean thicknesses [H2 ~> m2 or kg2 m-4].
  real :: haA, haB, haL, haR     ! Arithmetic mean thicknesses [H ~> m or kg m-2].
  real :: dzg2A, dzg2B  ! Squares of geometric mean vertical layer extents [Z2 ~> m2].
  real :: dzaA, dzaB    ! Arithmetic mean vertical layer extents [Z ~> m].
  real :: dzaL, dzaR    ! Temporary vertical layer extents [Z ~> m]
  real :: wtA, wtB      ! Unnormalized weights of the slopes above and below [H3 ~> m3 or kg3 m-6]
  real :: wtL, wtR      ! Unnormalized weights of the slopes to the left and right [H3 Z ~> m4 or kg3 m-5]
  real :: drdx, drdy    ! Zonal and meridional density gradients [R L-1 ~> kg m-4].
  real :: drdz          ! Vertical density gradient [R Z-1 ~> kg m-4].
  real :: dz_harm       ! Harmonic mean layer vertical extent [Z ~> m].
  real :: c2_dz_u(SZIB_(G),SZK_(GV)+1) ! Wave speed squared divided by dz at u-points [L2 Z-1 T-2 ~> m s-2]
  real :: c2_dz_v(SZI_(G),SZK_(GV)+1)  ! Wave speed squared divided by dz at v-points [L2 Z-1 T-2 ~> m s-2]
  real :: dzN2_u(SZIB_(G),SZK_(GV)+1) ! Vertical extent times N2 at interfaces above u-points times
                        ! rescaling factors from vertical to horizontal distances [L2 Z-1 T-2 ~> m s-2]
  real :: dzN2_v(SZI_(G),SZK_(GV)+1)  ! Vertical extent times N2 at interfaces above v-points times
                        ! rescaling factors from vertical to horizontal distances [L2 Z-1 T-2 ~> m s-2]
  real :: Sfn_est       ! A preliminary estimate (before limiting) of the overturning
                        ! streamfunction [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: Sfn_unlim_u(SZIB_(G),SZK_(GV)+1) ! Volume streamfunction for u-points [Z L2 T-1 ~> m3 s-1]
  real :: Sfn_unlim_v(SZI_(G),SZK_(GV)+1)  ! Volume streamfunction for v-points [Z L2 T-1 ~> m3 s-1]
  real :: slope2_Ratio_u(SZIB_(G),SZK_(GV)+1) ! The ratio of the slope squared to slope_max squared [nondim]
  real :: slope2_Ratio_v(SZI_(G),SZK_(GV)+1)  ! The ratio of the slope squared to slope_max squared [nondim]
  real :: Sfn_in_h      ! The overturning streamfunction [H L2 T-1 ~> m3 s-1 or kg s-1] (note that
                        ! the units are different from other Sfn vars).
  real :: Sfn_safe      ! The streamfunction that goes linearly back to 0 at the surface
                        ! [H L2 T-1 ~> m3 s-1 or kg s-1].  This is a good value to use when the
                        ! slope is so large as to be meaningless, usually due to weak stratification.
  real :: Slope         ! The slope of density surfaces, calculated in a way that is always
                        ! between -1 and 1 after undoing dimensional scaling, [Z L-1 ~> nondim]
  real :: mag_grad2     ! The squared magnitude of the 3-d density gradient [R2 L-2 ~> kg2 m-8].
  real :: I_slope_max2  ! The inverse of slope_max squared [L2 Z-2 ~> nondim].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: hn_2          ! Half of h_neglect [H ~> m or kg m-2].
  real :: h_neglect2    ! h_neglect^2 [H2 ~> m2 or kg2 m-4].
  real :: dz_neglect    ! A thickness [Z ~> m], that is so small it is usually lost
                        ! in roundoff and can be neglected [Z ~> m].
  real :: dz_neglect2   ! dz_neglect^2 [Z2 ~> m2]
  real :: G_scale       ! The gravitational acceleration times a unit conversion
                        ! factor [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.
  logical :: find_work  ! If true, find the change in energy due to the fluxes.
  integer :: nk_linear  ! The number of layers over which the streamfunction goes to 0.
  real :: G_rho0        ! g/Rho0 [L2 R-1 Z-1 T-2 ~> m4 kg-1 s-2].
  real :: Rho_avg       ! The in situ density averaged to an interface [R ~> kg m-3]
  real :: N2_floor      ! A floor for N2 to avoid degeneracy in the elliptic solver
                        ! times unit conversion factors [L2 Z-2 T-2 ~> s-2]
  real :: N2_unlim      ! An unlimited estimate of the buoyancy frequency
                        ! times unit conversion factors [L2 Z-2 T-2 ~> s-2]
  real :: Z_to_H        ! A conversion factor from heights to thicknesses, perhaps based on
                        ! a spatially variable local density [H Z-1 ~> nondim or kg m-3]
  real :: diag_sfn_x(SZIB_(G),SZJ_(G),SZK_(GV)+1)       ! Diagnostic of the x-face streamfunction
                                                        ! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: diag_sfn_unlim_x(SZIB_(G),SZJ_(G),SZK_(GV)+1) ! Diagnostic of the x-face streamfunction before
                                                        ! applying limiters [Z L2 T-1 ~> m3 s-1]
  real :: diag_sfn_y(SZI_(G),SZJB_(G),SZK_(GV)+1)       ! Diagnostic of the y-face streamfunction
                                                        ! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: diag_sfn_unlim_y(SZI_(G),SZJB_(G),SZK_(GV)+1) ! Diagnostic of the y-face streamfunction before
                                                        ! applying limiters [Z L2 T-1 ~> m3 s-1]
                                                        ! applying limiters [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, allocatable :: skeb_gm_work(:,:)                ! Temp array to hold GM work for SKEB
  real, allocatable :: skeb_ebt_norm2(:,:)              ! Used to normalize EBT for SKEB

  logical :: present_slope_x, present_slope_y, calc_derivatives
  integer, dimension(2) :: EOSdom_u  ! The shifted I-computational domain to use for equation of
                                     ! state calculations at u-points.
  integer, dimension(2) :: EOSdom_v  ! The shifted i-computational domain to use for equation of
                                     ! state calculations at v-points.
  integer, dimension(2) :: EOSdom_h1 ! The shifted i-computational domain to use for equation of
                                     ! state calculations at h points with 1 extra halo point
  logical :: use_stanley, skeb_use_gm
  integer :: is, ie, js, je, nz, IsdB, halo
  integer :: i, j, k
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke ; IsdB = G%IsdB

  I4dt = 0.25 / dt
  I_slope_max2 = 1.0 / (CS%slope_max**2)

  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect**2 ; hn_2 = 0.5*h_neglect
  dz_neglect = GV%dZ_subroundoff ; dz_neglect2 = dz_neglect**2
  if (GV%Boussinesq) G_rho0 = GV%g_Earth / GV%Rho0
  N2_floor = CS%N2_floor

  use_EOS = associated(tv%eqn_of_state)
  present_slope_x = PRESENT(slope_x)
  present_slope_y = PRESENT(slope_y)

  use_stanley = CS%use_stanley_gm

  skeb_use_gm = .false.
  if (present(STOCH)) skeb_use_gm = STOCH%skeb_use_gm
  if (skeb_use_gm) then
    allocate(skeb_gm_work(is:ie,js:je), source=0.)
    allocate(skeb_ebt_norm2(is:ie,js:je), source=0.)
  endif

  nk_linear = max(GV%nkml, 1)

  do concurrent (k=1:nz+1, j=G%jsd:G%jed, i=G%isdB:G%iedB)
    Slope_x_PE(i,j,k) = 0.0
    hN2_x_PE(i,j,k) = 0.0
  enddo
  do concurrent (k=1:nz+1, j=G%jsdB:G%jedB, i=G%isd:G%ied)
    Slope_y_PE(i,j,k) = 0.0
    hN2_y_PE(i,j,k) = 0.0
  enddo

  find_work = allocated(MEKE%GM_src)
  find_work = (allocated(CS%GMwork) .or. find_work)
  find_work = (skeb_use_gm .or. find_work)

  if (use_EOS) then
    halo = 1 ! Default halo to fill is 1
    !$omp target enter data map(to: h, tv%T, tv%S)
    !$omp target enter data map(alloc: T, S)
    call vert_fill_TS(h, tv%T, tv%S, CS%kappa_smooth*dt, T, S, G, GV, US, halo, larger_h_denom=.true.)
    !$omp target exit data map(from: T, S)
    !$omp target exit data map(release: h, tv%T, tv%S)
  endif

  ! Rescale the thicknesses, perhaps using the specific volume.
  !$omp target enter data map(to: tv, tv%SpV_avg) map(alloc: dz)
  call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=1, do_offload=.true.)
  !$omp target exit data map(release: tv, tv%SpV_avg) map(from: dz)

  if (CS%use_FGNV_streamfn .and. .not. associated(cg1)) call MOM_error(FATAL, &
       "cg1 must be associated when using FGNV streamfunction.")

  ! Find the maximum and minimum permitted streamfunction.
  do concurrent (j=js-1:je+1, i=is-1:ie+1)
    h_avail_rsum(i,j,1) = 0.0
    pres(i,j,1) = 0.0
    if (associated(tv%p_surf)) then ; pres(i,j,1) = tv%p_surf(i,j) ; endif

    h_avail(i,j,1) = max(I4dt*G%areaT(i,j)*(h(i,j,1)-GV%Angstrom_H),0.0)
    h_avail_rsum(i,j,2) = h_avail(i,j,1)
    h_frac(i,j,1) = 1.0
    pres(i,j,2) = pres(i,j,1) + (GV%g_Earth*GV%H_to_RZ) * h(i,j,1)
  enddo
  do concurrent (j=js-1:je+1)
    do k=2,nz ; do concurrent (i=is-1:ie+1)
      h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
      h_avail_rsum(i,j,k+1) = h_avail_rsum(i,j,k) + h_avail(i,j,k)
      h_frac(i,j,k) = 0.0 ; if (h_avail(i,j,k) > 0.0) &
        h_frac(i,j,k) = h_avail(i,j,k) / h_avail_rsum(i,j,k+1)
      pres(i,j,K+1) = pres(i,j,K) + (GV%g_Earth*GV%H_to_RZ) * h(i,j,k)
    enddo ; enddo
  enddo
  do concurrent (j=js:je, i=is-1:ie)
    uhtot(I,j) = 0.0 ; Work_u(I,j) = 0.0
  enddo
  do concurrent (J=js-1:je, i=is:ie)
    vhtot(i,J) = 0.0 ; Work_v(i,J) = 0.0
  enddo

  if (CS%id_sfn_x > 0) then ; diag_sfn_x(:,:,1) = 0.0 ; diag_sfn_x(:,:,nz+1) = 0.0 ; endif
  if (CS%id_sfn_y > 0) then ; diag_sfn_y(:,:,1) = 0.0 ; diag_sfn_y(:,:,nz+1) = 0.0 ; endif
  if (CS%id_sfn_unlim_x > 0) then ; diag_sfn_unlim_x(:,:,1) = 0.0 ; diag_sfn_unlim_x(:,:,nz+1) = 0.0 ; endif
  if (CS%id_sfn_unlim_y > 0) then ; diag_sfn_unlim_y(:,:,1) = 0.0 ; diag_sfn_unlim_y(:,:,nz+1) = 0.0 ; endif

  EOSdom_u(1) = (is-1) - (G%IsdB-1) ; EOSdom_u(2) = ie - (G%IsdB-1)
  EOSdom_v(:) = EOS_domain(G%HI)
  EOSdom_h1(:) = EOS_domain(G%HI, halo=1)

  !$OMP parallel do default(none) shared(nz,is,ie,js,je,find_work,use_EOS,G,GV,US,pres,T,S, &
  !$OMP                                  nk_linear,IsdB,tv,h,h_neglect,e,dz,dz_neglect,dz_neglect2, &
  !$OMP                                  h_neglect2,hn_2,I_slope_max2,int_slope_u,KH_u,uhtot, &
  !$OMP                                  h_frac,h_avail_rsum,uhD,h_avail,Work_u,CS,slope_x,cg1, &
  !$OMP                                  diag_sfn_x,diag_sfn_unlim_x,N2_floor,EOSdom_u,EOSdom_h1, &
  !$OMP                                  Sfn_unlim_u_3D, &
  !$OMP                                  use_stanley,present_slope_x,G_rho0,Slope_x_PE,hN2_x_PE) &
  !$OMP                          private(drdiA,drdiB,drdkL,drdkR,pres_u,T_u,S_u,G_scale, &
  !$OMP                                  drho_dT_u,drho_dS_u,hg2A,hg2B,hg2L,hg2R,haA, &
  !$OMP                                  drho_dT_dT_h,scrap,pres_h,T_h,S_h,N2_unlim,  &
  !$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,  &
  !$OMP                                  dzg2A,dzg2B,dzaA,dzaB,dz_harm,Z_to_H, &
  !$OMP                                  drdx,mag_grad2,Slope,slope2_Ratio_u,dzN2_u,  &
  !$OMP                                  Sfn_unlim_u,Rho_avg,drdi_u,drdkDe_u,c2_dz_u, &
  !$OMP                                  Sfn_safe,Sfn_est,Sfn_in_h,calc_derivatives)
  do j=js,je
    do I=is-1,ie ; dzN2_u(I,1) = 0. ; dzN2_u(I,nz+1) = 0. ; enddo
    do K=nz,2,-1
      if (find_work .and. .not.(use_EOS)) then
        drdiA = 0.0 ; drdiB = 0.0
        drdkL = GV%Rlay(k) - GV%Rlay(k-1) ; drdkR = drdkL
      endif

      calc_derivatives = use_EOS .and. (k >= nk_linear) .and. &
         (find_work .or. .not. present_slope_x .or. CS%use_FGNV_streamfn .or. use_stanley)

      ! Calculate the zonal fluxes and gradients.
      if (calc_derivatives) then
        do I=is-1,ie
          pres_u(I) = 0.5*(pres(i,j,K) + pres(i+1,j,K))
          T_u(I) = 0.25*((T(i,j,k) + T(i+1,j,k)) + (T(i,j,k-1) + T(i+1,j,k-1)))
          S_u(I) = 0.25*((S(i,j,k) + S(i+1,j,k)) + (S(i,j,k-1) + S(i+1,j,k-1)))
        enddo
        call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u, drho_dS_u, &
                                      tv%eqn_of_state, EOSdom_u)
      endif
      if (use_stanley) then
        do i=is-1,ie+1
          pres_h(i) = pres(i,j,K)
          T_h(i) = 0.5*(T(i,j,k) + T(i,j,k-1))
          S_h(i) = 0.5*(S(i,j,k) + S(i,j,k-1))
        enddo

        ! The second line below would correspond to arguments
        !            drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP, &
        call calculate_density_second_derivs(T_h, S_h, pres_h, &
                     scrap, scrap, drho_dT_dT_h, scrap, scrap, &
                     tv%eqn_of_state, EOSdom_h1)
      endif

      do I=is-1,ie
        if (calc_derivatives) then
          ! Estimate the horizontal density gradients along layers.
          drdiA = drho_dT_u(I) * (T(i+1,j,k-1)-T(i,j,k-1)) + &
                  drho_dS_u(I) * (S(i+1,j,k-1)-S(i,j,k-1))
          drdiB = drho_dT_u(I) * (T(i+1,j,k)-T(i,j,k)) + &
                  drho_dS_u(I) * (S(i+1,j,k)-S(i,j,k))

          ! Estimate the vertical density gradients times the grid spacing.
          drdkL = (drho_dT_u(I) * (T(i,j,k)-T(i,j,k-1)) + &
                   drho_dS_u(I) * (S(i,j,k)-S(i,j,k-1)))
          drdkR = (drho_dT_u(I) * (T(i+1,j,k)-T(i+1,j,k-1)) + &
                   drho_dS_u(I) * (S(i+1,j,k)-S(i+1,j,k-1)))
          drdkDe_u(I,K) = (drdkR * e(i+1,j,K)) - (drdkL * e(i,j,K))
        elseif (find_work) then ! This is used in pure stacked SW mode
          drdkDe_u(I,K) = (drdkR * e(i+1,j,K)) - (drdkL * e(i,j,K))
        endif
        if (use_stanley) then
          ! Correction to the horizontal density gradient due to nonlinearity in
          ! the EOS rectifying SGS temperature anomalies
          drdiA = drdiA + 0.5 * ((drho_dT_dT_h(i+1) * tv%varT(i+1,j,k-1)) - &
                                (drho_dT_dT_h(i) * tv%varT(i,j,k-1)) )
          drdiB = drdiB + 0.5 * ((drho_dT_dT_h(i+1) * tv%varT(i+1,j,k)) - &
                                (drho_dT_dT_h(i) * tv%varT(i,j,k)) )
        endif
        if (find_work) drdi_u(I,k) = drdiB

        if (k > nk_linear) then
          if (use_EOS) then
            if (CS%use_FGNV_streamfn .or. find_work .or. .not.present_slope_x) then
              hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
              hg2R = h(i+1,j,k-1)*h(i+1,j,k) + h_neglect2
              haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
              haR = 0.5*(h(i+1,j,k-1) + h(i+1,j,k)) + h_neglect
              if (GV%Boussinesq) then
                dzaL = haL * GV%H_to_Z ; dzaR = haR * GV%H_to_Z
              elseif (GV%semi_Boussinesq) then
                dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
                dzaR = 0.5*(e(i+1,j,K-1) - e(i+1,j,K+1)) + dz_neglect
              else
                dzaL = 0.5*(dz(i,j,k-1) + dz(i,j,k)) + dz_neglect
                dzaR = 0.5*(dz(i+1,j,k-1) + dz(i+1,j,k)) + dz_neglect
              endif
              ! Use the harmonic mean thicknesses to weight the horizontal gradients.
              ! These unnormalized weights have been rearranged to minimize divisions.
              wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

              drdz = ((wtL * drdkL) + (wtR * drdkR)) / ((dzaL*wtL) + (dzaR*wtR))
              ! The expression for drdz above is mathematically equivalent to:
              !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
              !          ((hg2L/haL) + (hg2R/haR))
              hg2A = h(i,j,k-1)*h(i+1,j,k-1) + h_neglect2
              hg2B = h(i,j,k)*h(i+1,j,k) + h_neglect2
              haA = 0.5*(h(i,j,k-1) + h(i+1,j,k-1)) + h_neglect
              haB = 0.5*(h(i,j,k) + h(i+1,j,k)) + h_neglect

              if (GV%Boussinesq) then
                N2_unlim = drdz*G_rho0
              else
                N2_unlim = (GV%g_Earth*GV%RZ_to_H) * &
                           (((wtL * drdkL) + (wtR * drdkR)) / ((haL*wtL) + (haR*wtR)))
              endif

              dzg2A = dz(i,j,k-1)*dz(i+1,j,k-1) + dz_neglect2
              dzg2B = dz(i,j,k)*dz(i+1,j,k) + dz_neglect2
              dzaA = 0.5*(dz(i,j,k-1) + dz(i+1,j,k-1)) + dz_neglect
              dzaB = 0.5*(dz(i,j,k) + dz(i+1,j,k)) + dz_neglect
              ! dzN2_u is used with the FGNV streamfunction formulation
              dzN2_u(I,K) = (0.5 * ( dzg2A / dzaA + dzg2B / dzaB )) * max(N2_unlim, N2_floor)
              if (find_work .and. CS%GM_src_alt) &
                hN2_x_PE(I,j,k) = (0.5 * ( hg2A / haA + hg2B / haB )) * max(N2_unlim, N2_floor)
            endif

            if (present_slope_x) then
              Slope = slope_x(I,j,k)
              slope2_Ratio_u(I,K) = Slope**2 * I_slope_max2
            else
              ! Use the harmonic mean thicknesses to weight the horizontal gradients.
              ! These unnormalized weights have been rearranged to minimize divisions.
              wtA = hg2A*haB ; wtB = hg2B*haA
              ! This is the gradient of density along geopotentials.
              drdx = ((wtA * drdiA + wtB * drdiB) / (wtA + wtB) - &
                      drdz * (e(i,j,K)-e(i+1,j,K))) * G%IdxCu(I,j)

              ! This estimate of slope is accurate for small slopes, but bounded
              ! to be between -1 and 1.
              mag_grad2 = (US%Z_to_L*drdx)**2 + drdz**2
              if (mag_grad2 > 0.0) then
                Slope = drdx / sqrt(mag_grad2)
                slope2_Ratio_u(I,K) = Slope**2 * I_slope_max2
              else ! Just in case mag_grad2 = 0 ever.
                Slope = 0.0
                slope2_Ratio_u(I,K) = 1.0e20  ! Force the use of the safe streamfunction.
              endif
            endif

            ! Adjust real slope by weights that bias towards slope of interfaces
            ! that ignore density gradients along layers.
            Slope = (1.0 - int_slope_u(I,j,K)) * Slope + &
                    int_slope_u(I,j,K) * ((e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j))
            slope2_Ratio_u(I,K) = (1.0 - int_slope_u(I,j,K)) * slope2_Ratio_u(I,K)

            if (CS%MEKE_src_slope_bug) then
              Slope_x_PE(I,j,k) = MIN(Slope, CS%slope_max)
            else
              Slope_x_PE(I,j,k) = Slope
              if (Slope > CS%slope_max) Slope_x_PE(I,j,k) = CS%slope_max
              if (Slope < -CS%slope_max) Slope_x_PE(I,j,k) = -CS%slope_max
            endif
            if (CS%id_slope_x > 0) CS%diagSlopeX(I,j,k) = Slope

            ! Estimate the streamfunction at each interface [H L2 T-1 ~> m3 s-1 or kg s-1].
            Sfn_unlim_u(I,K) = -(KH_u(I,j,K)*G%dy_Cu(I,j))*Slope

            if (CS%use_meso_sfn_ANN) then
              Sfn_unlim_u(I,K) = Sfn_unlim_u(I,K) + Sfn_unlim_u_3D(I,j,K)
            endif

            ! Avoid moving dense water upslope from below the level of
            ! the bottom on the receiving side.
            if (Sfn_unlim_u(I,K) > 0.0) then ! The flow below this interface is positive.
              if (e(i,j,K) < e(i+1,j,nz+1)) then
                Sfn_unlim_u(I,K) = 0.0 ! This is not uhtot, because it may compensate for
                                ! deeper flow in very unusual cases.
              elseif (e(i+1,j,nz+1) > e(i,j,K+1)) then
                ! Scale the transport with the fraction of the donor layer above
                ! the bottom on the receiving side.
                Sfn_unlim_u(I,K) = Sfn_unlim_u(I,K) * ((e(i,j,K) - e(i+1,j,nz+1)) / &
                                         ((e(i,j,K) - e(i,j,K+1)) + dz_neglect))
              endif
            else
              if (e(i+1,j,K) < e(i,j,nz+1)) then ; Sfn_unlim_u(I,K) = 0.0
              elseif (e(i,j,nz+1) > e(i+1,j,K+1)) then
                Sfn_unlim_u(I,K) = Sfn_unlim_u(I,K) * ((e(i+1,j,K) - e(i,j,nz+1)) / &
                                       ((e(i+1,j,K) - e(i+1,j,K+1)) + dz_neglect))
              endif
            endif

          else ! .not. use_EOS
            if (present_slope_x) then
              Slope = slope_x(I,j,k)
            else
              Slope = (e(i+1,j,K)-e(i,j,K)) * G%IdxCu_OBCmask(I,j)
            endif
            if (CS%id_slope_x > 0) CS%diagSlopeX(I,j,k) = Slope
            Sfn_unlim_u(I,K) = -(KH_u(I,j,K)*G%dy_Cu(I,j))*Slope
            dzN2_u(I,K) = GV%g_prime(K)

            if (CS%use_meso_sfn_ANN) then
              Sfn_unlim_u(I,K) = Sfn_unlim_u(I,K) + Sfn_unlim_u_3D(I,j,K)

              ! Avoid moving dense water upslope from below the level of
              ! the bottom on the receiving side.
              if (Sfn_unlim_u(I,K) > 0.0) then ! The flow below this interface is positive.
                if (e(i,j,K) < e(i+1,j,nz+1)) then
                  Sfn_unlim_u(I,K) = 0.0 ! This is not uhtot, because it may compensate for
                                  ! deeper flow in very unusual cases.
                elseif (e(i+1,j,nz+1) > e(i,j,K+1)) then
                  ! Scale the transport with the fraction of the donor layer above
                  ! the bottom on the receiving side.
                  Sfn_unlim_u(I,K) = Sfn_unlim_u(I,K) * ((e(i,j,K) - e(i+1,j,nz+1)) / &
                                           ((e(i,j,K) - e(i,j,K+1)) + dz_neglect))
                endif
              else
                if (e(i+1,j,K) < e(i,j,nz+1)) then ; Sfn_unlim_u(I,K) = 0.0
                elseif (e(i,j,nz+1) > e(i+1,j,K+1)) then
                  Sfn_unlim_u(I,K) = Sfn_unlim_u(I,K) * ((e(i+1,j,K) - e(i,j,nz+1)) / &
                                         ((e(i+1,j,K) - e(i+1,j,K+1)) + dz_neglect))
                endif
              endif
            endif

          endif ! if (use_EOS)
        else ! if (k > nk_linear)
          dzN2_u(I,K) = N2_floor * dz_neglect
          Sfn_unlim_u(I,K) = 0.
        endif ! if (k > nk_linear)
        if (CS%id_sfn_unlim_x>0) diag_sfn_unlim_x(I,j,K) = Sfn_unlim_u(I,K)
      enddo ! i-loop
    enddo ! k-loop

    if (CS%use_FGNV_streamfn) then
      do k=1,nz ; do I=is-1,ie ; if (G%OBCmaskCu(I,j)>0.) then
        dz_harm = max( dz_neglect, &
              2. * dz(i,j,k) * dz(i+1,j,k) / ( ( dz(i,j,k) + dz(i+1,j,k) ) + dz_neglect ) )
        c2_dz_u(I,k) = CS%FGNV_scale * ( 0.5*( cg1(i,j) + cg1(i+1,j) ) )**2 / dz_harm
      endif ; enddo ; enddo

      ! Solve an elliptic equation for the streamfunction following Ferrari et al., 2010.
      do I=is-1,ie
        if (G%OBCmaskCu(I,j)>0.) then
          do K=2,nz
            Sfn_unlim_u(I,K) = (1. + CS%FGNV_scale) * Sfn_unlim_u(I,K)
          enddo
          call streamfn_solver(nz, c2_dz_u(I,:), dzN2_u(I,:), Sfn_unlim_u(I,:))
        else
          do K=2,nz
            Sfn_unlim_u(I,K) = 0.
          enddo
        endif
      enddo
    endif

    do K=nz,2,-1
      do I=is-1,ie

        if (allocated(tv%SpV_avg) .and. (find_work .or. (k > nk_linear)) ) then
          Rho_avg = ( ((h(i,j,k) + h(i,j,k-1)) + (h(i+1,j,k) + h(i+1,j,k-1))) + 4.0*hn_2 ) / &
                ( (((h(i,j,k)+hn_2) * tv%SpV_avg(i,j,k))   + ((h(i,j,k-1)+hn_2) * tv%SpV_avg(i,j,k-1))) + &
                  (((h(i+1,j,k)+hn_2)*tv%SpV_avg(i+1,j,k)) + ((h(i+1,j,k-1)+hn_2)*tv%SpV_avg(i+1,j,k-1))) )
          ! Use an average density to convert the volume streamfunction estimate into a mass streamfunction.
          Z_to_H = GV%RZ_to_H*Rho_avg
        else
          Z_to_H = GV%Z_to_H
        endif

        if (k > nk_linear) then
          if (use_EOS) then

            if (uhtot(I,j) <= 0.0) then
              ! The transport that must balance the transport below is positive.
              Sfn_safe = uhtot(I,j) * (1.0 - h_frac(i,j,k))
            else !  (uhtot(I,j) > 0.0)
              Sfn_safe = uhtot(I,j) * (1.0 - h_frac(i+1,j,k))
            endif

            ! Determine the actual streamfunction at each interface.
            Sfn_est = (Z_to_H*Sfn_unlim_u(I,K) + slope2_Ratio_u(I,K)*Sfn_safe) / (1.0 + slope2_Ratio_u(I,K))
          else  ! When use_EOS is false, the layers are constant density.
            Sfn_est = Z_to_H*Sfn_unlim_u(I,K)
          endif

          ! Make sure that there is enough mass above to allow the streamfunction
          ! to satisfy the boundary condition of 0 at the surface.
          Sfn_in_H = min(max(Sfn_est, -h_avail_rsum(i,j,K)), h_avail_rsum(i+1,j,K))

          ! The actual transport is limited by the mass available in the two
          ! neighboring grid cells.
          uhD(I,j,k) = max(min((Sfn_in_H - uhtot(I,j)), h_avail(i,j,k)), &
                           -h_avail(i+1,j,k))

          if (CS%id_sfn_x>0) diag_sfn_x(I,j,K) = diag_sfn_x(I,j,K+1) + uhD(I,j,k)
!         sfn_x(I,j,K) = max(min(Sfn_in_h, uhtot(I,j)+h_avail(i,j,k)), &
!                            uhtot(I,j)-h_avail(i+1,j,K))
!         sfn_slope_x(I,j,K) = max(uhtot(I,j)-h_avail(i+1,j,k), &
!                                  min(uhtot(I,j)+h_avail(i,j,k), &
!               min(h_avail_rsum(i+1,j,K), max(-h_avail_rsum(i,j,K), &
!               (KH_u(I,j,K)*G%dy_Cu(I,j)) * ((e(i,j,K)-e(i+1,j,K))*G%IdxCu(I,j)) )) ))
        else ! k <= nk_linear
          ! Balance the deeper flow with a return flow uniformly distributed
          ! though the remaining near-surface layers.  This is the same as
          ! using Sfn_safe above.  There is no need to apply the limiters in
          ! this case.
          if (uhtot(I,j) <= 0.0) then
            uhD(I,j,k) = -uhtot(I,j) * h_frac(i,j,k)
          else !  (uhtot(I,j) > 0.0)
            uhD(I,j,k) = -uhtot(I,j) * h_frac(i+1,j,k)
          endif

          if (CS%id_sfn_x>0) diag_sfn_x(I,j,K) = diag_sfn_x(I,j,K+1) + uhD(I,j,k)
!         if (sfn_slope_x(I,j,K+1) <= 0.0) then
!           sfn_slope_x(I,j,K) = sfn_slope_x(I,j,K+1) * (1.0 - h_frac(i,j,k))
!         else
!           sfn_slope_x(I,j,K) = sfn_slope_x(I,j,K+1) * (1.0 - h_frac(i+1,j,k))
!         endif

        endif

        uhtot(I,j) = uhtot(I,j) + uhD(I,j,k)

        if (find_work) then
          !   This is the energy tendency based on the original profiles, and does
          ! not include any nonlinear terms due to a finite time step (which would
          ! involve interactions between the fluxes through the different faces.
          !   A second order centered estimate is used for the density transferred
          ! between water columns.

          if (allocated(tv%SpV_avg)) then
            G_scale = GV%H_to_RZ * GV%g_Earth / Rho_avg
          else
            G_scale = GV%g_Earth * GV%H_to_Z
          endif

          Work_u(I,j) = Work_u(I,j) + G_scale * &
            ( uhtot(I,j) * drdkDe_u(I,K) - &
              (uhD(I,j,k) * drdi_u(I,k)) * 0.25 * &
              ((e(i,j,K) + e(i,j,K+1)) + (e(i+1,j,K) + e(i+1,j,K+1))) )
        endif

      enddo
    enddo ! end of k-loop
  enddo ! end of j-loop

  ! Calculate the meridional fluxes and gradients.

  !$OMP parallel do default(none) shared(nz,is,ie,js,je,find_work,use_EOS,G,GV,US,pres,T,S,dz, &
  !$OMP                                  nk_linear,IsdB,tv,h,h_neglect,e,dz_neglect,dz_neglect2, &
  !$OMP                                  h_neglect2,int_slope_v,KH_v,vhtot,h_frac,h_avail_rsum, &
  !$OMP                                  I_slope_max2,vhD,h_avail,Work_v,CS,slope_y,cg1,hn_2,&
  !$OMP                                  diag_sfn_y,diag_sfn_unlim_y,N2_floor,EOSdom_v,use_stanley,&
  !$OMP                                  Sfn_unlim_v_3D, &
  !$OMP                                  present_slope_y,G_rho0,Slope_y_PE,hN2_y_PE)  &
  !$OMP                          private(drdjA,drdjB,drdkL,drdkR,pres_v,T_v,S_v,S_h,S_hr,    &
  !$OMP                                  drho_dT_v,drho_dS_v,hg2A,hg2B,hg2L,hg2R,haA,G_scale, &
  !$OMP                                  drho_dT_dT_h,drho_dT_dT_hr,scrap,pres_h,T_h,T_hr,   &
  !$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,pres_hr, &
  !$OMP                                  dzg2A,dzg2B,dzaA,dzaB,dz_harm,Z_to_H, &
  !$OMP                                  drdy,mag_grad2,Slope,slope2_Ratio_v,dzN2_v,N2_unlim, &
  !$OMP                                  Sfn_unlim_v,Rho_avg,drdj_v,drdkDe_v,c2_dz_v, &
  !$OMP                                  Sfn_safe,Sfn_est,Sfn_in_h,calc_derivatives)
  do J=js-1,je
    do K=nz,2,-1
      if (find_work .and. .not.(use_EOS)) then
        drdjA = 0.0 ; drdjB = 0.0
        drdkL = GV%Rlay(k) - GV%Rlay(k-1) ; drdkR = drdkL
      endif

      calc_derivatives = use_EOS .and. (k >= nk_linear) .and. &
         (find_work .or. .not. present_slope_y .or. CS%use_FGNV_streamfn .or. use_stanley)

      if (calc_derivatives) then
        do i=is,ie
          pres_v(i) = 0.5*(pres(i,j,K) + pres(i,j+1,K))
          T_v(i) = 0.25*((T(i,j,k) + T(i,j+1,k)) + (T(i,j,k-1) + T(i,j+1,k-1)))
          S_v(i) = 0.25*((S(i,j,k) + S(i,j+1,k)) + (S(i,j,k-1) + S(i,j+1,k-1)))
        enddo
        call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v, drho_dS_v, &
                                      tv%eqn_of_state, EOSdom_v)
      endif
      if (use_stanley) then
        do i=is,ie
          pres_h(i) = pres(i,j,K)
          T_h(i) = 0.5*(T(i,j,k) + T(i,j,k-1))
          S_h(i) = 0.5*(S(i,j,k) + S(i,j,k-1))

          pres_hr(i) = pres(i,j+1,K)
          T_hr(i) = 0.5*(T(i,j+1,k) + T(i,j+1,k-1))
          S_hr(i) = 0.5*(S(i,j+1,k) + S(i,j+1,k-1))
        enddo

        ! The second line below would correspond to arguments
        !            drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, drho_dT_dP, &
        call calculate_density_second_derivs(T_h, S_h, pres_h, &
                     scrap, scrap, drho_dT_dT_h, scrap, scrap, &
                     tv%eqn_of_state, EOSdom_v)
        call calculate_density_second_derivs(T_hr, S_hr, pres_hr, &
                     scrap, scrap, drho_dT_dT_hr, scrap, scrap, &
                     tv%eqn_of_state, EOSdom_v)
      endif
      do i=is,ie
        if (calc_derivatives) then
          ! Estimate the horizontal density gradients along layers.
          drdjA = drho_dT_v(i) * (T(i,j+1,k-1)-T(i,j,k-1)) + &
                  drho_dS_v(i) * (S(i,j+1,k-1)-S(i,j,k-1))
          drdjB = drho_dT_v(i) * (T(i,j+1,k)-T(i,j,k)) + &
                  drho_dS_v(i) * (S(i,j+1,k)-S(i,j,k))

          ! Estimate the vertical density gradients times the grid spacing.
          drdkL = (drho_dT_v(i) * (T(i,j,k)-T(i,j,k-1)) + &
                   drho_dS_v(i) * (S(i,j,k)-S(i,j,k-1)))
          drdkR = (drho_dT_v(i) * (T(i,j+1,k)-T(i,j+1,k-1)) + &
                   drho_dS_v(i) * (S(i,j+1,k)-S(i,j+1,k-1)))
          drdkDe_v(i,K) =  (drdkR * e(i,j+1,K)) - (drdkL * e(i,j,K))
        elseif (find_work) then ! This is used in pure stacked SW mode
          drdkDe_v(i,K) =  (drdkR * e(i,j+1,K)) - (drdkL * e(i,j,K))
        endif
        if (use_stanley) then
          ! Correction to the horizontal density gradient due to nonlinearity in
          ! the EOS rectifying SGS temperature anomalies
          drdjA = drdjA + 0.5 * ((drho_dT_dT_hr(i) * tv%varT(i,j+1,k-1)) - &
                                (drho_dT_dT_h(i) * tv%varT(i,j,k-1)) )
          drdjB = drdjB + 0.5 * ((drho_dT_dT_hr(i) * tv%varT(i,j+1,k)) - &
                                (drho_dT_dT_h(i) * tv%varT(i,j,k)) )
        endif

        if (find_work) drdj_v(i,k) = drdjB

        if (k > nk_linear) then
          if (use_EOS) then
            if (CS%use_FGNV_streamfn .or. find_work .or. .not. present_slope_y) then
              hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
              hg2R = h(i,j+1,k-1)*h(i,j+1,k) + h_neglect2
              haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
              haR = 0.5*(h(i,j+1,k-1) + h(i,j+1,k)) + h_neglect

              if (GV%Boussinesq) then
                dzaL = haL * GV%H_to_Z ; dzaR = haR * GV%H_to_Z
              elseif (GV%semi_Boussinesq) then
                dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
                dzaR = 0.5*(e(i,j+1,K-1) - e(i,j+1,K+1)) + dz_neglect
              else
                dzaL = 0.5*(dz(i,j,k-1) + dz(i,j,k)) + dz_neglect
                dzaR = 0.5*(dz(i,j+1,k-1) + dz(i,j+1,k)) + dz_neglect
              endif
              ! Use the harmonic mean thicknesses to weight the horizontal gradients.
              ! These unnormalized weights have been rearranged to minimize divisions.
              wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

              drdz = ((wtL * drdkL) + (wtR * drdkR)) / ((dzaL*wtL) + (dzaR*wtR))
              ! The expression for drdz above is mathematically equivalent to:
              !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
              !          ((hg2L/haL) + (hg2R/haR))
              hg2A = h(i,j,k-1)*h(i,j+1,k-1) + h_neglect2
              hg2B = h(i,j,k)*h(i,j+1,k) + h_neglect2
              haA = 0.5*(h(i,j,k-1) + h(i,j+1,k-1)) + h_neglect
              haB = 0.5*(h(i,j,k) + h(i,j+1,k)) + h_neglect

              if (GV%Boussinesq) then
                N2_unlim = drdz*G_rho0
              else
                N2_unlim = (GV%g_Earth*GV%RZ_to_H) * &
                           (((wtL * drdkL) + (wtR * drdkR)) / ((haL*wtL) + (haR*wtR)))
              endif

              dzg2A = dz(i,j,k-1)*dz(i,j+1,k-1) + dz_neglect2
              dzg2B = dz(i,j,k)*dz(i,j+1,k) + dz_neglect2
              dzaA = 0.5*(dz(i,j,k-1) + dz(i,j+1,k-1)) + dz_neglect
              dzaB = 0.5*(dz(i,j,k) + dz(i,j+1,k)) + dz_neglect

              ! dzN2_v is used with the FGNV streamfunction formulation
              dzN2_v(i,K) = (0.5*( dzg2A / dzaA + dzg2B / dzaB )) * max(N2_unlim, N2_floor)
              if (find_work .and. CS%GM_src_alt) &
                hN2_y_PE(i,J,k) = (0.5*( hg2A / haA + hg2B / haB )) * max(N2_unlim, N2_floor)
            endif
            if (present_slope_y) then
              Slope = slope_y(i,J,k)
              slope2_Ratio_v(i,K) = Slope**2 * I_slope_max2
            else
              ! Use the harmonic mean thicknesses to weight the horizontal gradients.
              ! These unnormalized weights have been rearranged to minimize divisions.
              wtA = hg2A*haB ; wtB = hg2B*haA
              ! This is the gradient of density along geopotentials.
              drdy = ((wtA * drdjA + wtB * drdjB) / (wtA + wtB) - &
                      drdz * (e(i,j,K)-e(i,j+1,K))) * G%IdyCv(i,J)

              ! This estimate of slope is accurate for small slopes, but bounded
              ! to be between -1 and 1.
              mag_grad2 = (US%Z_to_L*drdy)**2 + drdz**2
              if (mag_grad2 > 0.0) then
                Slope = drdy / sqrt(mag_grad2)
                slope2_Ratio_v(i,K) = Slope**2 * I_slope_max2
              else ! Just in case mag_grad2 = 0 ever.
                Slope = 0.0
                slope2_Ratio_v(i,K) = 1.0e20  ! Force the use of the safe streamfunction.
              endif
            endif

            ! Adjust real slope by weights that bias towards slope of interfaces
            ! that ignore density gradients along layers.
            Slope = (1.0 - int_slope_v(i,J,K)) * Slope + &
                    int_slope_v(i,J,K) * ((e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J))
            slope2_Ratio_v(i,K) = (1.0 - int_slope_v(i,J,K)) * slope2_Ratio_v(i,K)

            if (CS%MEKE_src_slope_bug) then
              Slope_y_PE(i,J,k) = MIN(Slope, CS%slope_max)
            else
              Slope_y_PE(i,J,k) = Slope
              if (Slope > CS%slope_max) Slope_y_PE(i,J,k) = CS%slope_max
              if (Slope < -CS%slope_max) Slope_y_PE(i,J,k) = -CS%slope_max
            endif
            if (CS%id_slope_y > 0) CS%diagSlopeY(I,j,k) = Slope

            Sfn_unlim_v(i,K) = -((KH_v(i,J,K)*G%dx_Cv(i,J))*Slope)

            if (CS%use_meso_sfn_ANN) then
              Sfn_unlim_v(i,K) = Sfn_unlim_v(i,K) + Sfn_unlim_v_3D(i,J,k)
            endif

            ! Avoid moving dense water upslope from below the level of
            ! the bottom on the receiving side.
            if (Sfn_unlim_v(i,K) > 0.0) then ! The flow below this interface is positive.
              if (e(i,j,K) < e(i,j+1,nz+1)) then
                Sfn_unlim_v(i,K) = 0.0 ! This is not vhtot, because it may compensate for
                                ! deeper flow in very unusual cases.
              elseif (e(i,j+1,nz+1) > e(i,j,K+1)) then
                ! Scale the transport with the fraction of the donor layer above
                ! the bottom on the receiving side.
                Sfn_unlim_v(i,K) = Sfn_unlim_v(i,K) * ((e(i,j,K) - e(i,j+1,nz+1)) / &
                                         ((e(i,j,K) - e(i,j,K+1)) + dz_neglect))
              endif
            else
              if (e(i,j+1,K) < e(i,j,nz+1)) then ; Sfn_unlim_v(i,K) = 0.0
              elseif (e(i,j,nz+1) > e(i,j+1,K+1)) then
                Sfn_unlim_v(i,K) = Sfn_unlim_v(i,K) * ((e(i,j+1,K) - e(i,j,nz+1)) / &
                                       ((e(i,j+1,K) - e(i,j+1,K+1)) + dz_neglect))
              endif
            endif

          else ! .not. use_EOS
            if (present_slope_y) then
              Slope = slope_y(i,J,k)
            else
              Slope = (e(i,j+1,K)-e(i,j,K)) * G%IdyCv_OBCmask(i,J)
            endif
            if (CS%id_slope_y > 0) CS%diagSlopeY(I,j,k) = Slope
            Sfn_unlim_v(i,K) = -((KH_v(i,J,K)*G%dx_Cv(i,J))*Slope)
            dzN2_v(i,K) = GV%g_prime(K)

            if (CS%use_meso_sfn_ANN) then
              Sfn_unlim_v(i,K) = Sfn_unlim_v(i,K) + Sfn_unlim_v_3D(i,J,k)

              ! Avoid moving dense water upslope from below the level of
              ! the bottom on the receiving side.
              if (Sfn_unlim_v(i,K) > 0.0) then ! The flow below this interface is positive.
                if (e(i,j,K) < e(i,j+1,nz+1)) then
                  Sfn_unlim_v(i,K) = 0.0 ! This is not vhtot, because it may compensate for
                                  ! deeper flow in very unusual cases.
                elseif (e(i,j+1,nz+1) > e(i,j,K+1)) then
                  ! Scale the transport with the fraction of the donor layer above
                  ! the bottom on the receiving side.
                  Sfn_unlim_v(i,K) = Sfn_unlim_v(i,K) * ((e(i,j,K) - e(i,j+1,nz+1)) / &
                                           ((e(i,j,K) - e(i,j,K+1)) + dz_neglect))
                endif
              else
                if (e(i,j+1,K) < e(i,j,nz+1)) then ; Sfn_unlim_v(i,K) = 0.0
                elseif (e(i,j,nz+1) > e(i,j+1,K+1)) then
                  Sfn_unlim_v(i,K) = Sfn_unlim_v(i,K) * ((e(i,j+1,K) - e(i,j,nz+1)) / &
                                         ((e(i,j+1,K) - e(i,j+1,K+1)) + dz_neglect))
                endif
              endif
            endif

          endif ! if (use_EOS)
        else ! if (k > nk_linear)
          dzN2_v(i,K) = N2_floor * dz_neglect
          Sfn_unlim_v(i,K) = 0.
        endif ! if (k > nk_linear)
        if (CS%id_sfn_unlim_y>0) diag_sfn_unlim_y(i,J,K) = Sfn_unlim_v(i,K)
      enddo ! i-loop
    enddo ! k-loop

    if (CS%use_FGNV_streamfn) then
      do k=1,nz ; do i=is,ie ; if (G%OBCmaskCv(i,J)>0.) then
        dz_harm = max( dz_neglect, &
              2. * dz(i,j,k) * dz(i,j+1,k) / ( ( dz(i,j,k) + dz(i,j+1,k) ) + dz_neglect ) )
        c2_dz_v(i,k) = CS%FGNV_scale * ( 0.5*( cg1(i,j) + cg1(i,j+1) ) )**2 / dz_harm
      endif ; enddo ; enddo

      ! Solve an elliptic equation for the streamfunction following Ferrari et al., 2010.
      do i=is,ie
        if (G%OBCmaskCv(i,J)>0.) then
          do K=2,nz
            Sfn_unlim_v(i,K) = (1. + CS%FGNV_scale) * Sfn_unlim_v(i,K)
          enddo
          call streamfn_solver(nz, c2_dz_v(i,:), dzN2_v(i,:), Sfn_unlim_v(i,:))
        else
          do K=2,nz
            Sfn_unlim_v(i,K) = 0.
          enddo
        endif
      enddo
    endif

    do K=nz,2,-1
      do i=is,ie
        if (allocated(tv%SpV_avg) .and. (find_work .or. (k > nk_linear)) ) then
          Rho_avg = ( ((h(i,j,k) + h(i,j,k-1)) + (h(i,j+1,k) + h(i,j+1,k-1))) + 4.0*hn_2 ) / &
              ( (((h(i,j,k)+hn_2) * tv%SpV_avg(i,j,k))   + ((h(i,j,k-1)+hn_2) * tv%SpV_avg(i,j,k-1))) + &
                (((h(i,j+1,k)+hn_2)*tv%SpV_avg(i,j+1,k)) + ((h(i,j+1,k-1)+hn_2)*tv%SpV_avg(i,j+1,k-1))) )
          ! Use an average density to convert the volume streamfunction estimate into a mass streamfunction.
          Z_to_H = GV%RZ_to_H*Rho_avg
        else
          Z_to_H = GV%Z_to_H
        endif

        if (k > nk_linear) then
          if (use_EOS) then

            if (vhtot(i,J) <= 0.0) then
              ! The transport that must balance the transport below is positive.
              Sfn_safe = vhtot(i,J) * (1.0 - h_frac(i,j,k))
            else !  (vhtot(I,j) > 0.0)
              Sfn_safe = vhtot(i,J) * (1.0 - h_frac(i,j+1,k))
            endif

            ! Find the actual streamfunction at each interface.
            Sfn_est = (Z_to_H*Sfn_unlim_v(i,K) + slope2_Ratio_v(i,K)*Sfn_safe) / (1.0 + slope2_Ratio_v(i,K))
          else  ! When use_EOS is false, the layers are constant density.
            Sfn_est = Z_to_H*Sfn_unlim_v(i,K)
          endif

          ! Make sure that there is enough mass above to allow the streamfunction
          ! to satisfy the boundary condition of 0 at the surface.
          Sfn_in_H = min(max(Sfn_est, -h_avail_rsum(i,j,K)), h_avail_rsum(i,j+1,K))

          ! The actual transport is limited by the mass available in the two
          ! neighboring grid cells.
          vhD(i,J,k) = max(min((Sfn_in_H - vhtot(i,J)), h_avail(i,j,k)), -h_avail(i,j+1,k))

          if (CS%id_sfn_y>0) diag_sfn_y(i,J,K) = diag_sfn_y(i,J,K+1) + vhD(i,J,k)
!         sfn_y(i,J,K) = max(min(Sfn_in_h, vhtot(i,J)+h_avail(i,j,k)), &
!                            vhtot(i,J)-h_avail(i,j+1,k))
!         sfn_slope_y(i,J,K) = max(vhtot(i,J)-h_avail(i,j+1,k), &
!                                  min(vhtot(i,J)+h_avail(i,j,k), &
!               min(h_avail_rsum(i,j+1,K), max(-h_avail_rsum(i,j,K), &
!               (KH_v(i,J,K)*G%dx_Cv(i,J)) * ((e(i,j,K)-e(i,j+1,K))*G%IdyCv(i,J)) )) ))
        else ! k <= nk_linear
          ! Balance the deeper flow with a return flow uniformly distributed
          ! though the remaining near-surface layers.  This is the same as
          ! using Sfn_safe above.  There is no need to apply the limiters in
          ! this case.
          if (vhtot(i,J) <= 0.0) then
            vhD(i,J,k) = -vhtot(i,J) * h_frac(i,j,k)
          else !  (vhtot(i,J) > 0.0)
            vhD(i,J,k) = -vhtot(i,J) * h_frac(i,j+1,k)
          endif

          if (CS%id_sfn_y>0) diag_sfn_y(i,J,K) = diag_sfn_y(i,J,K+1) + vhD(i,J,k)
!         if (sfn_slope_y(i,J,K+1) <= 0.0) then
!           sfn_slope_y(i,J,K) = sfn_slope_y(i,J,K+1) * (1.0 - h_frac(i,j,k))
!         else
!           sfn_slope_y(i,J,K) = sfn_slope_y(i,J,K+1) * (1.0 - h_frac(i,j+1,k))
!         endif
        endif

        vhtot(i,J) = vhtot(i,J)  + vhD(i,J,k)

        if (find_work) then
          !   This is the energy tendency based on the original profiles, and does
          ! not include any nonlinear terms due to a finite time step (which would
          ! involve interactions between the fluxes through the different faces.
          !   A second order centered estimate is used for the density transferred
          ! between water columns.

          if (allocated(tv%SpV_avg)) then
            G_scale = GV%H_to_RZ * GV%g_Earth / Rho_avg
          else
            G_scale = GV%g_Earth * GV%H_to_Z
          endif

          Work_v(i,J) = Work_v(i,J) + G_scale * &
            ( vhtot(i,J) * drdkDe_v(i,K) - &
             (vhD(i,J,k) * drdj_v(i,k)) * 0.25 * &
             ((e(i,j,K) + e(i,j,K+1)) + (e(i,j+1,K) + e(i,j+1,K+1))) )
        endif

      enddo
    enddo ! end of k-loop
  enddo ! end of j-loop

  ! In layer 1, enforce the boundary conditions that Sfn(z=0) = 0.0
  if (.not.find_work .or. .not.(use_EOS)) then
    do j=js,je ; do I=is-1,ie ; uhD(I,j,1) = -uhtot(I,j) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhD(i,J,1) = -vhtot(i,J) ; enddo ; enddo
  else
    EOSdom_u(1) = (is-1) - (G%IsdB-1) ; EOSdom_u(2) = ie - (G%IsdB-1)
    !$OMP parallel do default(shared) private(pres_u,T_u,S_u,drho_dT_u,drho_dS_u,drdiB,G_scale)
    do j=js,je
      if (use_EOS) then
        do I=is-1,ie
          pres_u(I) = 0.5*(pres(i,j,1) + pres(i+1,j,1))
          T_u(I) = 0.5*(T(i,j,1) + T(i+1,j,1))
          S_u(I) = 0.5*(S(i,j,1) + S(i+1,j,1))
        enddo
        call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u, drho_dS_u, &
                                      tv%eqn_of_state, EOSdom_u )
      endif
      do I=is-1,ie
        uhD(I,j,1) = -uhtot(I,j)

        G_scale = GV%g_Earth * GV%H_to_Z
        if (use_EOS) then
          drdiB = drho_dT_u(I) * (T(i+1,j,1)-T(i,j,1)) + &
                  drho_dS_u(I) * (S(i+1,j,1)-S(i,j,1))
          if (allocated(tv%SpV_avg)) then
            G_scale = GV%H_to_RZ * GV%g_Earth * &
                ( ( ((h(i,j,1)+hn_2) * tv%SpV_avg(i,j,1)) + ((h(i+1,j,1)+hn_2) * tv%SpV_avg(i+1,j,1)) ) / &
                  ( (h(i,j,1) + h(i+1,j,1)) + 2.0*hn_2 ) )
          endif
        endif
        if (CS%use_GM_work_bug) then
          Work_u(I,j) = Work_u(I,j) + G_scale * &
              ( (uhD(I,j,1) * drdiB) * 0.25 * &
                ((e(i,j,1) + e(i,j,2)) + (e(i+1,j,1) + e(i+1,j,2))) )
        else
          Work_u(I,j) = Work_u(I,j) - G_scale * &
              ( (uhD(I,j,1) * drdiB) * 0.25 * &
                ((e(i,j,1) + e(i,j,2)) + (e(i+1,j,1) + e(i+1,j,2))) )
        endif
      enddo
    enddo

    EOSdom_v(:) = EOS_domain(G%HI)
    !$OMP parallel do default(shared) private(pres_v,T_v,S_v,drho_dT_v,drho_dS_v,drdjB,G_scale)
    do J=js-1,je
      if (use_EOS) then
        do i=is,ie
          pres_v(i) = 0.5*(pres(i,j,1) + pres(i,j+1,1))
          T_v(i) = 0.5*(T(i,j,1) + T(i,j+1,1))
          S_v(i) = 0.5*(S(i,j,1) + S(i,j+1,1))
        enddo
        call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v, drho_dS_v, &
                                      tv%eqn_of_state, EOSdom_v)
      endif
      do i=is,ie
        vhD(i,J,1) = -vhtot(i,J)

        G_scale = GV%g_Earth * GV%H_to_Z
        if (use_EOS) then
          drdjB = drho_dT_v(i) * (T(i,j+1,1)-T(i,j,1)) + &
                  drho_dS_v(i) * (S(i,j+1,1)-S(i,j,1))
          if (allocated(tv%SpV_avg)) then
            G_scale = GV%H_to_RZ * GV%g_Earth * &
                ( ( ((h(i,j,1)+hn_2) * tv%SpV_avg(i,j,1)) + ((h(i,j+1,1)+hn_2) * tv%SpV_avg(i,j+1,1)) ) / &
                  ( (h(i,j,1) + h(i,j+1,1)) + 2.0*hn_2 ) )
          endif
        endif
        Work_v(i,J) = Work_v(i,J) - G_scale * &
            ( (vhD(i,J,1) * drdjB) * 0.25 * &
              ((e(i,j,1) + e(i,j,2)) + (e(i,j+1,1) + e(i,j+1,2))) )
      enddo
    enddo
  endif

  if (find_work) then ; do j=js,je ; do i=is,ie
    ! Note that the units of Work_v and Work_u are [R Z L4 T-3 ~> W], while Work_h is in [R Z L2 T-3 ~> W m-2].
    Work_h = 0.5 * G%IareaT(i,j) * &
      ((Work_u(I-1,j) + Work_u(I,j)) + (Work_v(i,J-1) + Work_v(i,J)))
    if (allocated(CS%GMwork)) CS%GMwork(i,j) = Work_h
    if (.not. CS%GM_src_alt) then ; if (allocated(MEKE%GM_src)) then
      MEKE%GM_src(i,j) = MEKE%GM_src(i,j) + Work_h
    endif ; endif
    if (skeb_use_gm) then
      skeb_gm_work(i,j)   = STOCH%skeb_gm_coef * Work_h
      skeb_ebt_norm2(i,j) = 0.0
      do k=1,nz
        skeb_ebt_norm2(i,j) = skeb_ebt_norm2(i,j) + h(i,j,k) * VarMix%ebt_struct(i,j,k)**2
      enddo
      skeb_ebt_norm2(i,j) = GV%H_to_RZ * (skeb_ebt_norm2(i,j) + h_neglect)
    endif
  enddo ; enddo ; endif

  if (skeb_use_gm) then
    ! This block spreads the GM work down through the column using the ebt vertical structure, squared.
    ! Note the sign convention.
    do k=1,nz ; do j=js,je ; do i=is,ie
      STOCH%skeb_diss(i,j,k) = STOCH%skeb_diss(i,j,k) - skeb_gm_work(i,j) * &
                               VarMix%ebt_struct(i,j,k)**2 / skeb_ebt_norm2(i,j)
    enddo ; enddo ; enddo
  endif

  if (find_work .and. CS%GM_src_alt) then ; if (allocated(MEKE%GM_src)) then
    if (CS%MEKE_src_answer_date >= 20240601) then
      do j=js,je ; do i=is,ie ; do k=nz,1,-1
        PE_release_h = -0.25 * GV%H_to_RZ * &
                         ( ((KH_u(I,j,k)*(Slope_x_PE(I,j,k)**2) * hN2_x_PE(I,j,k)) + &
                            (Kh_u(I-1,j,k)*(Slope_x_PE(I-1,j,k)**2) * hN2_x_PE(I-1,j,k))) + &
                           ((Kh_v(i,J,k)*(Slope_y_PE(i,J,k)**2) * hN2_y_PE(i,J,k)) + &
                            (Kh_v(i,J-1,k)*(Slope_y_PE(i,J-1,k)**2) * hN2_y_PE(i,J-1,k))) )
        MEKE%GM_src(i,j) = MEKE%GM_src(i,j) + PE_release_h
      enddo ; enddo ; enddo
    else
      do j=js,je ; do i=is,ie ; do k=nz,1,-1
        PE_release_h = -0.25 * GV%H_to_RZ * &
                           ((KH_u(I,j,k)*(Slope_x_PE(I,j,k)**2) * hN2_x_PE(I,j,k)) + &
                            (Kh_u(I-1,j,k)*(Slope_x_PE(I-1,j,k)**2) * hN2_x_PE(I-1,j,k)) + &
                            (Kh_v(i,J,k)*(Slope_y_PE(i,J,k)**2) * hN2_y_PE(i,J,k)) + &
                            (Kh_v(i,J-1,k)*(Slope_y_PE(i,J-1,k)**2) * hN2_y_PE(i,J-1,k)))
        MEKE%GM_src(i,j) = MEKE%GM_src(i,j) + PE_release_h
      enddo ; enddo ; enddo
    endif

    if (CS%debug) then
      call hchksum(MEKE%GM_src, 'MEKE%GM_src', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      call uvchksum("KH_[uv]", Kh_u, Kh_v, G%HI, unscale=US%L_to_m**2*US%s_to_T, &
                    scalar_pair=.true.)
      call uvchksum("Slope_[xy]_PE", Slope_x_PE, Slope_y_PE, G%HI, unscale=US%Z_to_L)
      call uvchksum("hN2_[xy]_PE", hN2_x_PE, hN2_y_PE, G%HI, unscale=GV%H_to_mks*US%L_to_Z**2*US%s_to_T**2, &
                    scalar_pair=.true.)
    endif
  endif ; endif

  if (CS%id_slope_x > 0) call post_data(CS%id_slope_x, CS%diagSlopeX, CS%diag)
  if (CS%id_slope_y > 0) call post_data(CS%id_slope_y, CS%diagSlopeY, CS%diag)
  if (CS%id_sfn_x > 0) call post_data(CS%id_sfn_x, diag_sfn_x, CS%diag)
  if (CS%id_sfn_y > 0) call post_data(CS%id_sfn_y, diag_sfn_y, CS%diag)
  if (CS%id_sfn_unlim_x > 0) call post_data(CS%id_sfn_unlim_x, diag_sfn_unlim_x, CS%diag)
  if (CS%id_sfn_unlim_y > 0) call post_data(CS%id_sfn_unlim_y, diag_sfn_unlim_y, CS%diag)

end subroutine thickness_diffuse_full

end submodule MOM_thickness_diffuse_s
