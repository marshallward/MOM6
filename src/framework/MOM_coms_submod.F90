! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

#include "do_concurrent_compat.h"

!> Interfaces to non-domain-oriented communication subroutines, including the
!! MOM6 reproducing sums facility
submodule (MOM_coms) MOM_coms_submod

use, intrinsic :: iso_fortran_env, only : int64
use MOM_coms_infra,    only : PE_here, root_PE, num_PEs, set_rootPE, Set_PElist, Get_PElist
use MOM_coms_infra,    only : broadcast, field_chksum, MOM_infra_init, MOM_infra_end
use MOM_coms_infra,    only : sum_across_PEs, max_across_PEs, min_across_PEs
use MOM_coms_infra,    only : all_across_PEs, any_across_PEs
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_coms_infra,    only : sync_PEs

implicit none

integer(kind=int64), parameter :: prec = (2_int64)**46 !< The precision of each integer.
real, parameter :: r_prec=2.0**46  !< A real version of prec [nondim].
real, parameter :: I_prec=1.0/(2.0**46) !< The inverse of prec [nondim].
integer, parameter :: max_count_prec=2**(63-46)-1
                              !< The number of values that can be added together
                              !! with the current value of prec before there will
                              !! be roundoff problems.

integer, parameter :: ni=6    !< The number of long integers to use to represent
                              !< a real number.
real, parameter, dimension(ni) :: &
  pr = (/ r_prec**2, r_prec, 1.0, 1.0/r_prec, 1.0/r_prec**2, 1.0/r_prec**3 /)
    !< An array of the real precision of each of the integers in arbitrary units [a]
real, parameter, dimension(ni) :: &
  I_pr = (/ 1.0/r_prec**2, 1.0/r_prec, 1.0, r_prec, r_prec**2, r_prec**3 /)
    !< An array of the inverse of the real precision of each of the integers in arbitrary units [a-1]
real, parameter :: max_efp_float = pr(1) * (2.**63 - 1.)

!$omp declare target(pr, I_pr)

contains

!**!!> Decompose one real into its 6 signed EFP bin contributions. PURE so it can be
!**!!! called from do concurrent (stdpar offloads it automatically) as well as from
!**!!! the host. The arithmetic is exact (pr/I_pr are powers of two), so the bins are
!**!!! bit-identical CPU<->GPU. Scalar outputs avoid a per-thread local array, so no
!**!!! forced inlining is required. NaN/overflow are reported via flags rather than
!**!!! the module-level error logicals, so the routine is side-effect free.
!**!pure subroutine efp_decompose(r, e1, e2, e3, e4, e5, e6, rmag, is_nan, is_ovf)
!**!  real,                intent(in)  :: r      !< The real number being decomposed [a]
!**!  integer(kind=int64), intent(out) :: e1     !< Signed contribution to EFP bin 1
!**!  integer(kind=int64), intent(out) :: e2     !< Signed contribution to EFP bin 2
!**!  integer(kind=int64), intent(out) :: e3     !< Signed contribution to EFP bin 3
!**!  integer(kind=int64), intent(out) :: e4     !< Signed contribution to EFP bin 4
!**!  integer(kind=int64), intent(out) :: e5     !< Signed contribution to EFP bin 5
!**!  integer(kind=int64), intent(out) :: e6     !< Signed contribution to EFP bin 6
!**!  real,                intent(out) :: rmag   !< abs(r), or 0 if r is NaN/Inf [a]
!**!  integer,             intent(out) :: is_nan !< 1 if r is a NaN or Inf, else 0
!**!  integer,             intent(out) :: is_ovf !< 1 if abs(r) has no EFP representation, else 0
!**!
!**!  real :: rs  ! The remaining value to add, in arbitrary units [a]
!**!  integer(kind=int64) :: ival
!**!  integer :: sgn
!**!
!**!  e1 = 0_int64 ; e2 = 0_int64 ; e3 = 0_int64
!**!  e4 = 0_int64 ; e5 = 0_int64 ; e6 = 0_int64
!**!  rmag = 0.0 ; is_nan = 0 ; is_ovf = 0
!**!
!**!  if ((r >= 1e30) .eqv. (r < 1e30)) then ; is_nan = 1 ; return ; endif
!**!  sgn = 1 ; if (r < 0.0) sgn = -1
!**!  rs = abs(r) ; rmag = rs
!**!
!**!  ! Abort if the number has no EFP representation
!**!  if (rs > max_efp_float) then ; is_ovf = 1 ; return ; endif
!**!
!**!  ival = int(rs*I_pr(1), kind=int64) ; rs = rs - ival*pr(1) ; e1 = sgn*ival
!**!  ival = int(rs*I_pr(2), kind=int64) ; rs = rs - ival*pr(2) ; e2 = sgn*ival
!**!  ival = int(rs*I_pr(3), kind=int64) ; rs = rs - ival*pr(3) ; e3 = sgn*ival
!**!  ival = int(rs*I_pr(4), kind=int64) ; rs = rs - ival*pr(4) ; e4 = sgn*ival
!**!  ival = int(rs*I_pr(5), kind=int64) ; rs = rs - ival*pr(5) ; e5 = sgn*ival
!**!  ival = int(rs*I_pr(6), kind=int64) ; rs = rs - ival*pr(6) ; e6 = sgn*ival
!**!end subroutine efp_decompose

!> Decompose one real into its 6 signed EFP bin contributions. PURE so it can be
!! called from do concurrent (stdpar offloads it automatically) as well as from
!! the host. The arithmetic is exact (pr/I_pr are powers of two), so the bins are
!! bit-identical CPU<->GPU. Scalar outputs avoid a per-thread local array, so no
!! forced inlining is required. NaN/overflow are reported via flags rather than
!! the module-level error logicals, so the routine is side-effect free.
pure subroutine efp_decompose(r, e, rmag, is_nan, is_ovf)
  !$omp declare target
  real,                intent(in)  :: r      !< The real number being decomposed [a]
  integer(kind=int64), intent(out) :: e(ni) !< Signed contribution to EFP bin 1
  real,                intent(out) :: rmag   !< abs(r), or 0 if r is NaN/Inf [a]
  integer,             intent(out) :: is_nan !< 1 if r is a NaN or Inf, else 0
  integer,             intent(out) :: is_ovf !< 1 if abs(r) has no EFP representation, else 0

  real :: rs  ! The remaining value to add, in arbitrary units [a]
  integer(kind=int64) :: ival
  integer :: sgn
  integer :: n

  e(:) = 0
  rmag = 0.0 ; is_nan = 0 ; is_ovf = 0

  if ((r >= 1e30) .eqv. (r < 1e30)) then ; is_nan = 1 ; return ; endif
  sgn = 1 ; if (r < 0.0) sgn = -1
  rs = abs(r) ; rmag = rs

  ! Abort if the number has no EFP representation
  if (rs > max_efp_float) then ; is_ovf = 1 ; return ; endif

  do n=1,ni
    ival = int(rs*I_pr(n), kind=int64)
    rs = rs - ival*pr(n)
    e(n) = sgn*ival
  enddo
end subroutine efp_decompose

!> Increment an EFP number with a real number without doing any carrying of
!! of overflows and using only minimal error checking.
module subroutine increment_ints_faster(int_sum, r, max_mag_term)
  integer(kind=int64), dimension(ni), intent(inout) :: int_sum  !< The array of EFP integers being incremented
  real,                           intent(in)    :: r        !< The real number being added in arbitrary units [a]
  real,                           intent(inout) :: max_mag_term !< A running maximum magnitude of the r's
                                                            !! in arbitrary units [a]

  ! This subroutine increments a number with another, both using the integer
  ! representation in real_to_ints, but without doing any carrying of overflow.
  ! The per-element decomposition is shared with the GPU kernel via efp_decompose.
  !integer(kind=int64) :: e1, e2, e3, e4, e5, e6
  integer(kind=int64) :: e(ni)
  real    :: rmag
  integer :: is_nan, is_ovf

  call efp_decompose(r, e, rmag, is_nan, is_ovf)

  if (is_nan /= 0) then ; NaN_error = .true. ; return ; endif
  if (rmag > abs(max_mag_term)) max_mag_term = r

  ! Abort if the number has no EFP representation
  if (is_ovf /= 0) then ; overflow_error = .true. ; return ; endif

  int_sum(:) = int_sum(:) + e(:)

end subroutine increment_ints_faster

!> Per-slice EFP bin accumulation over a 2d window, GPU-friendly. Equivalent to
!! looping increment_ints_faster over (is:ie, js:je), modulo: max_mag_term is
!! updated with the magnitude (>=0) rather than the signed last-winner. Only
!! consumer of max_mag_term is abs() in the overflow guard, so values are
!! unchanged; only the sign in one FATAL message can differ.
module subroutine increment_ints_2d(array, is, ie, js, je, descale, ints_sum, max_mag_term)
  real, dimension(:,:),               intent(in)    :: array  !< Input slice in arbitrary units [a]
  integer,                            intent(in)    :: is !< Starting i-index of the window (1-based)
  integer,                            intent(in)    :: ie !< Ending i-index of the window (1-based)
  integer,                            intent(in)    :: js !< Starting j-index of the window (1-based)
  integer,                            intent(in)    :: je !< Ending j-index of the window (1-based)
  real,                               intent(in)    :: descale !< unscale factor or 1.0 [a A-1 ~> 1]
  integer(kind=int64), dimension(ni), intent(inout) :: ints_sum !< EFP bins, incremented in place
  real,                               intent(inout) :: max_mag_term !< Running max magnitude [a]

  integer :: i, j
  integer(kind=int64) :: e(ni)
  integer(kind=int64) :: s(ni)
  real :: r, rmag, mmag
  integer :: inan, iovf, lnan, lovf

  mmag = abs(max_mag_term)
  inan = 0 ; iovf = 0

  ! The per-element decomposition is shared with the host path via efp_decompose;
  ! the six int64 bins reduce as scalars (libnvomp array-section reductions are
  ! avoided on purpose). array is copied in by stdpar if not already resident.
  do concurrent (j=js:je, i=is:ie) &
      DO_LOCALITY(local(r, e, rmag, lnan, lovf)) &
      DO_LOCALITY(reduce(+: ints_sum) reduce(max: mmag, inan, iovf))
    r = descale*array(i,j)

    call efp_decompose(r, e, rmag, lnan, lovf)

    inan = max(inan, lnan)
    iovf = max(iovf, lovf)

    if (rmag > mmag) mmag = rmag

    ints_sum(:) = ints_sum(:) + e(:)
  enddo

  max_mag_term = mmag
  if (inan /= 0) NaN_error = .true.
  if (iovf /= 0) overflow_error = .true.
end subroutine increment_ints_2d

end submodule MOM_coms_submod
