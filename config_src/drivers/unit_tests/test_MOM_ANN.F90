! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Driver for the artificial neural network unit tests.
program test_MOM_ANN

use, intrinsic :: iso_fortran_env, only : error_unit
use MOM_error_handler, only : set_skip_mpi
use MOM_ANN_tests, only : collect_tests
use testdrive, only : run_testsuite

implicit none

integer :: status ! Number of failed tests

call set_skip_mpi(.true.) ! This unit tests is not expecting MPI to be used

status = 0
call run_testsuite(collect_tests, error_unit, status, parallel=.false.)
if (status > 0) stop 1

end program test_MOM_ANN
