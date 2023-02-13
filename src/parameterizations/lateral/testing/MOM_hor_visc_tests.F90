module MOM_hor_visc_tests

use MOM_unit_testing, only : TestSuite

implicit none ; private

public :: run_hor_visc_tests

contains

subroutine test_hor_visc_init
end subroutine test_hor_visc_init

subroutine run_hor_visc_tests
  type(TestSuite) :: suite

  ! Build the test suite
  suite = TestSuite()

  call suite%add(test_hor_visc_init, "test_hor_visc_init")

  call suite%run()
end subroutine run_hor_visc_tests

end module MOM_hor_visc_tests
