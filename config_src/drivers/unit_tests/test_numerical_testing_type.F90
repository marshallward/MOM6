! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Driver for the numerical testing type unit tests.
program test_numerical_testing_type

use, intrinsic :: iso_fortran_env, only : error_unit
use testdrive, only : run_testsuite
use numerical_testing_type_tests, only : collect_tests

implicit none

integer :: status ! Number of failed tests

status = 0
call run_testsuite(collect_tests, error_unit, status, parallel=.false.)
if (status > 0) stop 1

end program test_numerical_testing_type
