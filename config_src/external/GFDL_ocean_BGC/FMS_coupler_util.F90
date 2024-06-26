module FMS_coupler_util

use coupler_types_mod, only : coupler_2d_bc_type

implicit none ; private

public :: extract_coupler_values, set_coupler_values

contains

!> Get element and index of a boundary condition
subroutine extract_coupler_values(BC_struc, BC_index, BC_element, array_out, ilb, jlb, &
                                  is, ie, js, je, conversion)
  integer,                   intent(in)  :: ilb !< Lower bounds
  integer,                   intent(in)  :: jlb !< Lower bounds
  real, dimension(ilb:,jlb:),intent(out) :: array_out !< The array being filled with the input values
  type(coupler_2d_bc_type),  intent(in)  :: BC_struc !< The type from which the data is being extracted
  integer,                   intent(in)  :: BC_index !< The boundary condition number being extracted
  integer,                   intent(in)  :: BC_element !< The element of the boundary condition being extracted
  integer,        optional,  intent(in)  :: is !< The i- limits of array_out to be filled
  integer,        optional,  intent(in)  :: ie !< The i- limits of array_out to be filled
  integer,        optional,  intent(in)  :: js !< The j- limits of array_out to be filled
  integer,        optional,  intent(in)  :: je !< The j- limits of array_out to be filled
  real,           optional,  intent(in)  :: conversion !< A number that every element is multiplied by

  array_out(:,:) = -1.
end subroutine extract_coupler_values

!> Set element and index of a boundary condition
subroutine set_coupler_values(array_in, BC_struc, BC_index, BC_element, ilb, jlb,&
                              is, ie, js, je, conversion)
  integer,                  intent(in)    :: ilb !< Lower bounds
  integer,                  intent(in)    :: jlb !< Lower bounds
  real, dimension(ilb:,jlb:), intent(in)  :: array_in !< The array containing the values to load into the BC
  type(coupler_2d_bc_type), intent(inout) :: BC_struc !< The type into which the data is being loaded
  integer,                  intent(in)    :: BC_index !< The boundary condition number being set
  integer,                  intent(in)    :: BC_element !< The element of the boundary condition being set
  integer,        optional, intent(in)    :: is !< The i- limits of array_out to be filled
  integer,        optional, intent(in)    :: ie !< The i- limits of array_out to be filled
  integer,        optional, intent(in)    :: js !< The j- limits of array_out to be filled
  integer,        optional, intent(in)    :: je !< The j- limits of array_out to be filled
  real,           optional, intent(in)    :: conversion !< A number that every element is multiplied by
end subroutine set_coupler_values

end module FMS_coupler_util
