! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Tests for artificial neural network inference.
module MOM_ANN_tests

use MOM_ANN, only : ANN_CS, ANN_allocate, ANN_apply, ANN_apply_vector_orig
use MOM_ANN, only : ANN_apply_vector_oi, ANN_apply_array_sio, ANN_end, ANN_random
use MOM_ANN, only : set_layer, set_input_normalization, set_output_normalization
use testdrive, only : check, error_type, new_unittest, unittest_type

implicit none ; private

public :: collect_tests

contains

!> Collect the neural network tests.
subroutine collect_tests(tests)
  type(unittest_type), allocatable, intent(out) :: tests(:) !< Named tests

  tests = [ &
    new_unittest("scalar identity", test_scalar_identity), &
    new_unittest("summation", test_summation), &
    new_unittest("vector identity", test_vector_identity), &
    new_unittest("rectifier", test_rectifier), &
    new_unittest("hidden layer, bias and normalization", test_hidden_layer), &
    new_unittest("four-layer inference variants", test_four_layers), &
    new_unittest("random inference equivalence", test_random_inference) ]
end subroutine collect_tests

!> Check a scalar identity network.
subroutine test_scalar_identity(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(ANN_CS) :: ann ! Network
  real :: result(1) ! Network output [A]

  call ANN_allocate(ann, 2, [1, 1])
  call set_layer(ann, 1, reshape([1.], [1, 1]), [0.], .false.)
  call ANN_apply([1.], result, ann)
  call ANN_end(ann)
  call check(error, result(1) == 1., "Scalar identity")
end subroutine test_scalar_identity

!> Check a network that sums its inputs.
subroutine test_summation(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(ANN_CS) :: ann ! Network
  real :: result(1) ! Network output [A]

  call ANN_allocate(ann, 2, [4, 1])
  call set_layer(ann, 1, reshape([1., 1., 1., 1.], [1, 4]), [0.], .false.)
  call ANN_apply([-1., 0., 1., 2.], result, ann)
  call ANN_end(ann)
  call check(error, result(1) == 2., "Summation")
end subroutine test_summation

!> Check a vector identity network.
subroutine test_vector_identity(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(ANN_CS) :: ann ! Network
  real :: result(3) ! Network output [A]

  call ANN_allocate(ann, 2, [3, 3])
  call set_layer(ann, 1, reshape([1., 0., 0., 0., 1., 0., 0., 0., 1.], [3, 3]), [0., 0., 0.], .false.)
  call ANN_apply([-1., 0., 1.], result, ann)
  call ANN_end(ann)
  call check(error, all(result(:) == [-1., 0., 1.]), "Vector identity")
end subroutine test_vector_identity

!> Check the rectifier activation.
subroutine test_rectifier(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(ANN_CS) :: ann ! Network
  real :: result(3) ! Network output [A]

  call ANN_allocate(ann, 2, [3, 3])
  call set_layer(ann, 1, reshape([1., 0., 0., 0., 1., 0., 0., 0., 1.], [3, 3]), [0., 0., 0.], .true.)
  call ANN_apply([-1., 0., 1.], result, ann)
  call ANN_end(ann)
  call check(error, all(result(:) == [0., 0., 1.]), "Rectifier")
end subroutine test_rectifier

!> Set the rectifying hidden layer and its summation layer.
subroutine set_hidden_layers(ann, bias)
  type(ANN_CS), intent(inout) :: ann !< Network
  real, intent(in) :: bias !< Hidden-layer bias [A]

  call set_layer(ann, 1, reshape([1., 0., 0., 0., 0., 1., 0., 0., &
                                  0., 0., 1., 0., 0., 0., 0., 1.], [4, 4]), &
                 [bias, bias, bias, bias], .true.)
  call set_layer(ann, 2, reshape([1., 1., 1., 1.], [1, 4]), [3.*bias], .false.)
end subroutine set_hidden_layers

!> Check hidden-layer inference with bias and normalization.
subroutine test_hidden_layer(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(ANN_CS) :: ann ! Network
  real :: result(1), without_bias, with_bias, normalized ! Network outputs [A]

  call ANN_allocate(ann, 3, [4, 4, 1])
  call set_hidden_layers(ann, 0.)
  call ANN_apply_vector_orig([-1., 0., 1., 2.], result, ann)
  without_bias = result(1)
  call set_hidden_layers(ann, 1.)
  call ANN_apply_vector_orig([-2., -1., 0., 1.], result, ann)
  with_bias = result(1)
  call set_input_normalization(ann, means=[4., 4., 4., 4.], norms=[0.5, 0.5, 0.5, 0.5])
  call set_output_normalization(ann, norms=[2.], means=[-10.])
  call ANN_apply_vector_orig([0., 2., 4., 6.], result, ann)
  normalized = result(1)
  call ANN_end(ann)

  call check(error, without_bias == 3., "Rectifier and summation")
  if (allocated(error)) return
  call check(error, with_bias == 6., "Rectifier, summation and bias")
  if (allocated(error)) return
  call check(error, normalized == 2., "Rectifier, summation, bias and normalization")
end subroutine test_hidden_layer

!> Check four-layer inference using the vector and array implementations.
subroutine test_four_layers(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  type(ANN_CS) :: ann ! Network
  real :: original(1), optimized(1), array_result(2,1) ! Network outputs [A]

  call ANN_allocate(ann, 4, [4, 4, 1, 1])
  call set_hidden_layers(ann, 1.)
  call set_layer(ann, 3, reshape([1.], [1, 1]), [0.], .false.)
  call set_input_normalization(ann, means=[4., 4., 4., 4.], norms=[0.5, 0.5, 0.5, 0.5])
  call set_output_normalization(ann, norms=[2.], means=[-10.])
  call ANN_apply_vector_orig([0., 2., 4., 6.], original, ann)
  call ANN_apply_vector_oi([0., 2., 4., 6.], optimized, ann)
  call ANN_apply_array_sio(2, reshape([0., 1., 2., 3., 4., 5., 6., 7.], [2, 4]), array_result, ann)
  call ANN_end(ann)

  call check(error, original(1) == 2., "Original vector inference")
  if (allocated(error)) return
  call check(error, optimized(1) == 2., "Optimized vector inference")
  if (allocated(error)) return
  call check(error, all(array_result(:,1) == [2., 5.]), "Array inference")
end subroutine test_four_layers

!> Compare inference implementations for seeded random networks and inputs.
subroutine test_random_inference(error)
  type(error_type), allocatable, intent(out) :: error !< Failure, if any
  integer, parameter :: max_layers = 10, num_samples = 1000, num_points = 20
  type(ANN_CS) :: ann ! Network
  integer, allocatable :: saved_seed(:), seed(:) ! Random-number generator state
  integer :: seed_size, widths(max_layers), num_layers, iter, i
  real :: draws(max_layers+1) ! Random network dimensions [nondim]
  real, allocatable :: input(:), original(:), optimized(:), inputs(:,:), outputs(:,:) ! Inputs and outputs [A]
  logical :: equal_results ! Whether all implementations agree
  character(len=100) :: message ! Failure context

  call random_seed(size=seed_size)
  allocate(saved_seed(seed_size), seed(seed_size))
  call random_seed(get=saved_seed)
  do i = 1, seed_size
    seed(i) = 12345 + i
  enddo
  call random_seed(put=seed)

  equal_results = .true.
  message = ""
  do iter = 1, num_samples
    call random_number(draws)
    num_layers = 2 + floor(draws(max_layers+1) * (max_layers - 1))
    widths(1:num_layers) = 1 + floor(draws(1:num_layers) * 8)
    call ANN_random(ann, num_layers, widths(1:num_layers))
    allocate(input(widths(1)), original(widths(num_layers)), optimized(widths(num_layers)))
    allocate(inputs(num_points,widths(1)), outputs(num_points,widths(num_layers)))
    call random_number(input)
    input(:) = (2. * input(:)) - 1.
    do i = 1, num_points
      inputs(i,:) = input(:)
    enddo
    call ANN_apply_vector_orig(input, original, ann)
    call ANN_apply_vector_oi(input, optimized, ann)
    call ANN_apply_array_sio(num_points, inputs, outputs, ann)
    equal_results = all(optimized(:) == original(:))
    do i = 1, num_points
      equal_results = equal_results .and. all(outputs(i,:) == original(:))
    enddo
    call ANN_end(ann)
    deallocate(input, original, optimized, inputs, outputs)
    if (.not. equal_results) then
      write(message, '("Inference variants differ for random network ",i0)') iter
      exit
    endif
  enddo
  call random_seed(put=saved_seed)
  call check(error, equal_results, trim(message))
end subroutine test_random_inference

end module MOM_ANN_tests
