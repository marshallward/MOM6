! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Test the string functions independently of model startup.
module MOM_string_functions_tests

use testdrive, only : check, error_type, new_unittest, unittest_type
use MOM_string_functions, only : left_int, left_ints, left_real, left_reals
use MOM_string_functions, only : ints_to_string, trim_trailing_commas, extractWord
use MOM_string_functions, only : extract_word, remove_spaces, extract_integer, extract_real

implicit none ; private

public :: collect_tests

contains

!> Collect the string-function tests.
subroutine collect_tests(tests)
  type(unittest_type), allocatable, intent(out) :: tests(:) !< Named tests

  tests = [ &
    new_unittest("format integers", test_format_integers), &
    new_unittest("format reals", test_format_reals), &
    new_unittest("integer strings", test_integer_strings), &
    new_unittest("extract words", test_extract_words), &
    new_unittest("remove spaces", test_remove_spaces), &
    new_unittest("extract integers", test_extract_integers), &
    new_unittest("extract reals", test_extract_reals) ]
end subroutine collect_tests

!> Check integer formatting.
subroutine test_format_integers(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  integer, parameter :: values(5) = [ -1, 1, 3, 3, 0 ] ! Test values

  call check(error, trim(left_int(-1)), "-1")
  if (allocated(error)) return
  call check(error, trim(left_ints(values(:))), "-1, 1, 3, 3, 0")
end subroutine test_format_integers

!> Check real formatting and separators.
subroutine test_format_reals(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  real, parameter :: values(8) = [ 0., 1., -2., 1.3, 3.e-11, 3.e-11, 3.e-11, -5.1e12 ] !< Test values [A]

  call check(error, trim(left_real(0.)), "0.0")
  if (allocated(error)) return
  call check(error, trim(left_reals(values(:))), "0.0, 1.0, -2.0, 1.3, 3*3.0E-11, -5.1E+12")
  if (allocated(error)) return
  call check(error, trim(left_reals(values(:), sep=" ")), "0.0 1.0 -2.0 1.3 3*3.0E-11 -5.1E+12")
  if (allocated(error)) return
  call check(error, trim(left_reals(values(:), sep=",")), "0.0,1.0,-2.0,1.3,3*3.0E-11,-5.1E+12")
end subroutine test_format_reals

!> Check fixed-width integer strings.
subroutine test_integer_strings(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  integer, parameter :: values(5) = [ -1, 1, 3, 3, 0 ] ! Test values

  call check(error, trim(ints_to_string(values(:), 5)), "_-0001_0001_0003_0003_0000")
  if (allocated(error)) return
  call check(error, trim(ints_to_string(values(2:), 2)), "_0001_0003")
  if (allocated(error)) return
  call check(error, trim(ints_to_string(values(:))), "_-0001_0001_0003")
end subroutine test_integer_strings

!> Check word extraction and trailing-comma removal.
subroutine test_extract_words(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any

  call check(error, trim(trim_trailing_commas("One, Two, Three, ")), "One, Two, Three")
  if (allocated(error)) return
  call check(error, trim(extractWord("One Two,Three", 1)), "One")
  if (allocated(error)) return
  call check(error, trim(extractWord("One Two,Three", 2)), "Two")
  if (allocated(error)) return
  call check(error, trim(extractWord("One Two,Three", 3)), "Three")
  if (allocated(error)) return
  call check(error, trim(extractWord("One Two,  Three", 3)), "Three")
  if (allocated(error)) return
  call check(error, trim(extractWord(" One Two,Three", 1)), "One")
  if (allocated(error)) return
  call check(error, trim(extract_word("One,Two,Three", ",", 3)), "Three")
  if (allocated(error)) return
  call check(error, trim(extract_word("One,Two,Three", ",", 4)), "")
end subroutine test_extract_words

!> Check removal of spaces, including empty input.
subroutine test_remove_spaces(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any

  call check(error, trim(remove_spaces("1 2 3")), "123")
  if (allocated(error)) return
  call check(error, trim(remove_spaces(" 1 2 3")), "123")
  if (allocated(error)) return
  call check(error, trim(remove_spaces("1 2 3 ")), "123")
  if (allocated(error)) return
  call check(error, trim(remove_spaces("123")), "123")
  if (allocated(error)) return
  call check(error, trim(remove_spaces(" ")), "")
  if (allocated(error)) return
  call check(error, trim(remove_spaces("")), "")
end subroutine test_remove_spaces

!> Check integer extraction and defaults.
subroutine test_extract_integers(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any

  call check(error, extract_integer("1", "", 1), 1)
  if (allocated(error)) return
  call check(error, extract_integer("1,2,3", ",", 1), 1)
  if (allocated(error)) return
  call check(error, extract_integer("1,2", ",", 2), 2)
  if (allocated(error)) return
  call check(error, extract_integer("1,2", ",", 3), 0)
  if (allocated(error)) return
  call check(error, extract_integer("1,2", ",", 4, 4), 4)
end subroutine test_extract_integers

!> Check real extraction and defaults.
subroutine test_extract_reals(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any

  call check(error, extract_real("1.", "", 1) == 1.)
  if (allocated(error)) return
  call check(error, extract_real("1.,2.,3.", ",", 1) == 1.)
  if (allocated(error)) return
  call check(error, extract_real("1.,2.", ",", 2) == 2.)
  if (allocated(error)) return
  call check(error, extract_real("1.,2.", ",", 3) == 0.)
  if (allocated(error)) return
  call check(error, extract_real("1.,2.", ",", 4, 4.) == 4.)
end subroutine test_extract_reals

end module MOM_string_functions_tests
