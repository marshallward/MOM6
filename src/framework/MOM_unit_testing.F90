! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

module MOM_unit_testing

use posix, only : chmod
use posix, only : sigsetjmp
use posix, only : sigjmp_buf

use MOM_coms, only : num_PEs, sync_PEs
use MOM_error_handler, only : is_root_pe
use MOM_error_handler, only : disable_fatal_errors
use MOM_error_handler, only : enable_fatal_errors
use MOM_error_handler, only : query_skip_mpi
use MOM_error_handler, only : assert

implicit none ; private

public :: string
public :: create_test_file
public :: delete_test_file
public :: TestSuite


!> String container type
type :: string
  character(len=:), allocatable :: s
    !< Internal character array of string
end type string


!> String constructor
interface string
  module procedure init_string_char
  module procedure init_string_int
end interface string


!> A generalized instance of a unit test.
type, abstract :: UnitTest
  private
  procedure(), nopass, pointer :: cleanup => null()
    !< Cleanup function to be run after proc
  character(len=:), allocatable :: name
    !< Unit test name (usually set to name of proc)
contains
  procedure :: run => run_unit_test_base
    !< Run the unit test
end type UnitTest


!> A unit test defined by a procedure.
type, extends(UnitTest) :: ProcedureUnitTest
  private
  procedure(), nopass, pointer :: proc => null()
    !< Unit test function/subroutine
  logical :: is_fatal
    !< True if proc() is expected to fail
contains
  procedure :: run => run_procedure_unit_test
    !< Run the unit test function, proc
end type ProcedureUnitTest


!> A unit test defined by a scalar value comparison.
type, extends(UnitTest) :: ScalarUnitTest
  private
  real :: result_value
    !< Scalar value being tested [arbitrary]
  real :: reference
    !< Reference scalar value [arbitrary]
  real :: tolerance
    !< Absolute comparison tolerance [arbitrary]
  logical :: has_tolerance = .false.
    !< True if scalar comparison uses an absolute tolerance
contains
  procedure :: run => run_scalar_unit_test
    !< Run the scalar value comparison test
end type ScalarUnitTest


!> Unit test constructor
interface UnitTest
  module procedure create_procedure_unit_test_basic
  module procedure create_procedure_unit_test_full
  module procedure create_scalar_unit_test
end interface UnitTest


!> Collection of unit tests
type :: TestSuite
  private
  type(UnitTestNode), pointer :: head => null()
    !< Head of the unit test linked list
  type(UnitTestNode), pointer :: tail => null()
    !< Tail of the unit test linked list (pre-allocated and unconfigured)

  ! Public API
  procedure(), nopass, pointer, public :: cleanup => null()
    !< Default cleanup function for unit tests in suite
contains
  private
  procedure :: add_basic => add_unit_test_basic
    !< Add a unit test without a cleanup function
  procedure :: add_full => add_unit_test_full
    !< Add a unit test with an explicit cleanup function
  generic, public :: add => add_basic, add_full
    !< Add a unit test to the test suite
  procedure, public :: add_scalar => add_unit_test_scalar
    !< Add a scalar value comparison test to the test suite
  procedure, public :: run => run_test_suite
    !< Run all unit tests in the suite
end type TestSuite


!> TestSuite constructor
interface TestSuite
  module procedure create_test_suite
end interface TestSuite


!> UnitTest node of TestSuite's linked list
type :: UnitTestNode
  private
  class(UnitTest), pointer :: test => null()
    !< Node contents
  type(UnitTestNode), pointer :: next => null()
    !< Pointer to next node in list
end type UnitTestNode

contains

!> Return a new procedure unit test without a cleanup function.
function create_procedure_unit_test_basic(proc, name, fatal) result(test)
  procedure() :: proc
    !< Subroutine which defines the unit test
  character(len=*), intent(in) :: name
    !< Name of the unit test
  logical, intent(in), optional :: fatal
    !< True if the test is expected to raise a FATAL error
  type(ProcedureUnitTest) :: test

  procedure(), pointer :: cleanup
  cleanup => null()

  test = create_procedure_unit_test_full(proc, name, fatal, cleanup)
end function create_procedure_unit_test_basic


!> Return a new procedure unit test with an explicit cleanup function.
function create_procedure_unit_test_full(proc, name, fatal, cleanup) &
    result(test)
  procedure() :: proc
    !< Subroutine which defines the unit test
  character(len=*), intent(in) :: name
    !< Name of the unit test
  logical, optional :: fatal
    !< True if the test is expected to raise a FATAL error
  procedure() :: cleanup
    !< Cleanup subroutine, called after test
  type(ProcedureUnitTest) :: test

  test%proc => proc
  test%name = name
  test%is_fatal = .false.
  if (present(fatal)) test%is_fatal = fatal
  test%cleanup => cleanup
end function create_procedure_unit_test_full


!> Return a new scalar value comparison unit test.
function create_scalar_unit_test(result_value, reference, name, tolerance, &
    cleanup) result(test)
  real, intent(in) :: result_value
    !< Scalar value being tested [arbitrary]
  real, intent(in) :: reference
    !< Reference scalar value [arbitrary]
  character(len=*), intent(in) :: name
    !< Name of the test
  real, optional, intent(in) :: tolerance
    !< Absolute comparison tolerance [arbitrary]
  procedure() :: cleanup
    !< Cleanup subroutine, called after test
  type(ScalarUnitTest) :: test

  test%cleanup => cleanup
  test%name = name
  test%result_value = result_value
  test%reference = reference
  test%tolerance = 0.
  test%has_tolerance = present(tolerance)
  if (present(tolerance)) test%tolerance = tolerance
end function create_scalar_unit_test


!> Base unit test run method.  Concrete unit test types override this method.
subroutine run_unit_test_base(test)
  class(UnitTest), intent(in) :: test

  call assert(.false., "MOM_unit_testing: run_unit_test_base called for " &
      // test%name)
end subroutine run_unit_test_base


!> Launch a procedure unit test with a custom cleanup procedure.
subroutine run_procedure_unit_test(test)
  class(ProcedureUnitTest), intent(in) :: test

  type(sigjmp_buf) :: env
  integer :: rc

  if (.not. query_skip_mpi()) then
    call sync_PEs

    ! FIXME: Some FATAL tests under MPI are unable to recover after jumpback,
    !   so we disable these tests for now.
    if (test%is_fatal .and. num_PEs() > 1) return
  endif

  if (test%is_fatal) then
    rc = sigsetjmp(env, 1)
    if (rc == 0) then
      call disable_fatal_errors(env)
      call test%proc
    endif
    call enable_fatal_errors
  else
    call test%proc
  endif

  if (associated(test%cleanup)) call test%cleanup
end subroutine run_procedure_unit_test


!> Run a scalar value comparison test.
subroutine run_scalar_unit_test(test)
  class(ScalarUnitTest), intent(in) :: test

  real :: err
    !< Absolute difference between test value and reference value [arbitrary]

  if (test%has_tolerance) then
    err = abs(test%result_value - test%reference)
    print '(2x, a, " = ", ES22.15, ", ref = ", ES22.15, ", abs err = ", &
        ES9.2, ", tol = ", ES9.2)', &
        test%name, test%result_value, test%reference, err, test%tolerance
    call assert(err <= test%tolerance, trim(test%name) &
        // " exceeds scalar tolerance")
  else
    print '(2x, a, " = ", ES22.15, ", ref = ", ES22.15)', &
        test%name, test%result_value, test%reference
    call assert(test%result_value == test%reference, trim(test%name) &
        // " differs from scalar reference")
  endif

  if (associated(test%cleanup)) call test%cleanup
end subroutine run_scalar_unit_test


!> Return a new test suite
function create_test_suite() result(suite)
  type(TestSuite) :: suite

  ! Setup the head node, but do not populate it
  allocate(suite%head)
  suite%tail => suite%head
end function create_test_suite


subroutine add_unit_test_basic(suite, test, name, fatal)
  class(TestSuite), intent(inout) :: suite
  procedure() :: test
  character(len=*), intent(in) :: name
  logical, intent(in), optional :: fatal

  procedure(), pointer :: cleanup

  cleanup => null()
  if (associated(suite%cleanup)) cleanup => suite%cleanup

  call add_unit_test_full(suite, test, name, fatal, cleanup)
end subroutine add_unit_test_basic


subroutine add_unit_test_full(suite, test, name, fatal, cleanup)
  class(TestSuite), intent(inout) :: suite
  procedure() :: test
  character(len=*), intent(in) :: name
  procedure() :: cleanup
  logical, intent(in), optional :: fatal

  type(ProcedureUnitTest), pointer :: utest

  allocate(utest)
  utest = UnitTest(test, name, fatal, cleanup)
  call append_unit_test(suite, utest)
end subroutine add_unit_test_full


subroutine add_unit_test_scalar(suite, result_value, reference, name, tolerance)
  class(TestSuite), intent(inout) :: suite
  real, intent(in) :: result_value
    !< Scalar value being tested [arbitrary]
  real, intent(in) :: reference
    !< Reference scalar value [arbitrary]
  character(len=*), intent(in) :: name
    !< Name of the test
  real, optional, intent(in) :: tolerance
    !< Absolute comparison tolerance [arbitrary]

  type(ScalarUnitTest), pointer :: utest
  procedure(), pointer :: cleanup

  cleanup => null()
  if (associated(suite%cleanup)) cleanup => suite%cleanup

  allocate(utest)
  utest = UnitTest(result_value, reference, name, tolerance, cleanup)
  call append_unit_test(suite, utest)
end subroutine add_unit_test_scalar


!> Append a unit test to a test suite.
subroutine append_unit_test(suite, test)
  class(TestSuite), intent(inout) :: suite
  class(UnitTest), pointer, intent(in) :: test
    !< Unit test to append

  type(UnitTestNode), pointer :: node

  suite%tail%test => test
  allocate(node)
  suite%tail%next => node
  suite%tail => node
end subroutine append_unit_test


subroutine run_test_suite(suite)
  class(TestSuite), intent(in) :: suite

  type(UnitTestNode), pointer :: node

  node => suite%head
  do while(associated(node%test))
    ! TODO: Capture FMS stdout/stderr
    print '(/a)', "=== "//node%test%name

    call node%test%run

    node => node%next
  enddo
end subroutine run_test_suite


!> Initialize string with a character array.
function init_string_char(c) result(str)
  character(len=*), dimension(:), intent(in) :: c
    !< List of character arrays
  type(string), dimension(size(c)) :: str
    !< String output

  integer :: i

  do i = 1, size(c)
    str(i)%s = c(i)
  enddo
end function init_string_char


!> Convert an integer to a string
function init_string_int(n) result(str)
  integer, intent(in) :: n
    !< Integer input
  type(string) :: str
    !< String output

  ! TODO: Estimate this with integer arithmetic
  character(1 + floor(log10(real(abs(n)))) + (1 - sign(1, n))/2) :: chr

  write(chr, '(i0)') n
  str = string(chr)
end function init_string_int


!> Create a text file for unit testing
subroutine create_test_file(filename, lines, mode)
  character(len=*), intent(in) :: filename
    !< Name of file to be created
  type(string), intent(in), optional :: lines(:)
    !< list of strings to write to file
  integer, optional, intent(in) :: mode
    !< Permissions of new file

  integer :: param_unit
  integer :: i
  integer :: rc
  logical :: sync

  if (is_root_PE()) then
    open(newunit=param_unit, file=filename, status='replace')
    if (present(lines)) then
      do i = 1, size(lines)
        write(param_unit, '(a)') lines(i)%s
      enddo
    endif
    close(param_unit)
    if (present(mode)) rc = chmod(filename, mode)
  endif
  if (.not. query_skip_mpi()) call sync_PEs
end subroutine create_test_file


!> Delete a file created during testing
subroutine delete_test_file(filename)
  character(len=*), intent(in) :: filename
    !< Name of file to be deleted

  logical :: is_file, is_open
  integer :: io_unit

  if (is_root_PE()) then
    inquire(file=filename, exist=is_file, opened=is_open, number=io_unit)

    if (is_file) then
      if (.not. is_open) open(newunit=io_unit, file=filename)
      close(io_unit, status='delete')
    endif
  endif
  if (.not. query_skip_mpi()) call sync_PEs
end subroutine delete_test_file

end module MOM_unit_testing
