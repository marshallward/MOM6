module MOM_hor_visc_tests

use MOM_diag_mediator, only : time_type
use MOM_diag_mediator, only : diag_ctrl

use MOM_domains, only : MOM_domains_init
use MOM_grid, only : ocean_grid_type
use MOM_grid, only : MOM_grid_init

use MOM_verticalGrid, only : verticalGrid_type

use MOM_unit_scaling, only : unit_scale_type

use MOM_file_parser, only : param_file_type
use MOM_file_parser, only : open_param_file
use MOM_file_parser, only : close_param_file

use MOM_hor_visc, only : hor_visc_init
use MOM_hor_visc, only : hor_visc_CS
use MOM_hor_visc, only : hor_visc_end

use MOM_unit_testing, only : TestSuite
use MOM_unit_testing, only : string
use MOM_unit_testing, only : create_test_file

implicit none ; private

public :: run_hor_visc_tests

! Durr
character(len=*), parameter :: param_filename = 'TEST_input'

contains

subroutine test_hor_visc_init
  type(string) :: lines(3)

  type(time_type) :: Time
  type(ocean_grid_type) :: G
  type(verticalGrid_type) :: GV
  type(unit_scale_type) :: US
  type(param_file_type) :: param
  type(diag_ctrl) :: diag
  type(hor_visc_CS) :: CS

  lines = [ &
      string('DT = 0.1'), &
      string('NIGLOBAL = 8'), &
      string('NJGLOBAL = 8') &
  ]
  call create_test_file('MOM_input', lines)

  call open_param_file('MOM_input', param)
  call MOM_domains_init(G%Domain, param)
  call MOM_grid_init(G, param)

  call hor_visc_init(Time, G, GV, US, param, diag, CS)
  call hor_visc_end(CS)

  call close_param_file(param)
end subroutine test_hor_visc_init

subroutine run_hor_visc_tests
  type(TestSuite) :: suite

  ! Build the test suite
  suite = TestSuite()

  call suite%add(test_hor_visc_init, "test_hor_visc_init")

  call suite%run()
end subroutine run_hor_visc_tests

end module MOM_hor_visc_tests
