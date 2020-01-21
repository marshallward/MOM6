!> Module for supporting the transformation of fields between different grids.
module MOM_array_transform

implicit none; private

public rotate_array
public rotate_quarter, rotate_half  ! Hopefully won't need these...

interface rotate_array
  module procedure rotate_array_real_2d
  module procedure rotate_array_real_3d
  module procedure rotate_array_real_4d
  module procedure rotate_array_integer
end interface rotate_array

interface rotate_quarter
  module procedure rotate_quarter_real
  module procedure rotate_quarter_integer
end interface rotate_quarter

interface rotate_half
  module procedure rotate_half_real, rotate_half_integer
end interface rotate_half

contains

subroutine rotate_array_real_4d(A, A_in, turns)
  real, intent(out) :: A(:,:,:,:)
  real, intent(in) :: A_in(:,:,:,:)
  integer, intent(in) :: turns

  integer :: n

  do n = lbound(A_in, 4), ubound(A_in, 4)
    call rotate_array(A(:,:,:,n), A_in(:,:,:,n), turns)
  enddo
end subroutine rotate_array_real_4d


subroutine rotate_array_real_3d(A, A_in, turns)
  real, intent(out) :: A(:,:,:)
  real, intent(in) :: A_in(:,:,:)
  integer, intent(in) :: turns

  integer :: k

  do k = lbound(A_in, 3), ubound(A_in, 3)
    call rotate_array(A(:,:,k), A_in(:,:,k), turns)
  enddo
end subroutine rotate_array_real_3d


subroutine rotate_array_real_2d(A, A_in, turns)
  real, intent(out) :: A(:,:)
  real, intent(in) :: A_in(:,:)
  integer, intent(in) :: turns

  select case (modulo(turns, 4))
    case(0)
      A = A_in
    case(1)
      A = rotate_quarter(A_in)
    case(2)
      A = rotate_half(A_in)
    case(3)
      A = rotate_quarter(A_in, clockwise=.true.)
  end select
end subroutine rotate_array_real_2d


subroutine rotate_array_integer(A, A_in, turns)
  integer, intent(out) :: A(:,:)
  integer, intent(in) :: A_in(:,:)
  integer, intent(in) :: turns

  select case (modulo(turns, 4))
    case(0)
      A = A_in
    case(1)
      A = rotate_quarter(A_in)
    case(2)
      A = rotate_half(A_in)
    case(3)
      A = rotate_quarter(A_in, clockwise=.true.)
  end select
end subroutine rotate_array_integer


function rotate_quarter_real(A_in, clockwise) result(A)
  real, intent(in)  :: A_in(:,:)
  logical, intent(in), optional :: clockwise
  real :: A(size(A_in, 2), size(A_in, 1))

  logical :: ccw    ! True if turn is counterclockwise (positive)
  integer :: m, n   ! Input array shape

  ccw = .true.
  if (present(clockwise)) ccw = .not. clockwise

  m = size(A_in, 1)
  n = size(A_in, 2)

  ! +90deg: B(i,j) = A(n-j,i)
  !                = transpose, then row reverse
  ! -90deg: B(i,j) = A(j,m-i)
  !                = row reverse, then transpose

  if (.not. ccw) then
    A = transpose(A_in(m:1:-1, :))
  else
    A = transpose(A_in)
  endif

  if (ccw) &
    A(:,:) = A(n:1:-1, :)

end function rotate_quarter_real


function rotate_quarter_integer(A_in, clockwise) result(A)
  integer, intent(in)  :: A_in(:,:)
  logical, intent(in), optional :: clockwise
  integer :: A(size(A_in, 2), size(A_in, 1))

  logical :: ccw    ! True if turn is counterclockwise (positive)
  integer :: m, n   ! Input array shape

  ccw = .true.
  if (present(clockwise)) ccw = .not. clockwise

  m = size(A_in, 1)
  n = size(A_in, 2)

  ! +90deg: B(i,j) = A(n-j,i)
  !                = transpose, then row reverse
  ! -90deg: B(i,j) = A(j,m-i)
  !                = row reverse, then transpose

  if (.not. ccw) then
    A = transpose(A_in(m:1:-1, :))
  else
    A = transpose(A_in)
  endif

  if (ccw) &
    A(:,:) = A(n:1:-1, :)

end function rotate_quarter_integer


function rotate_half_real(A_in) result(A)
  real, intent(in)  :: A_in(:,:)
  real :: A(size(A_in, 1), size(A_in, 2))

  integer :: m, n   ! Input array shape

  m = size(A_in, 1)
  n = size(A_in, 2)

  ! 180deg: B(i,j) = A(m-i,n-j)
  !                = row reversal + column reversal

  A(:,:) = A_in(m:1:-1, :)
  A(:,:) = A(:, n:1:-1)

end function rotate_half_real


function rotate_half_integer(A_in) result(A)
  integer, intent(in)  :: A_in(:,:)
  integer :: A(size(A_in, 1), size(A_in, 2))

  integer :: m, n   ! Input array shape

  m = size(A_in, 1)
  n = size(A_in, 2)

  ! 180deg: B(i,j) = A(m-i,n-j)
  !                = row reversal + column reversal

  A(:,:) = A_in(m:1:-1, :)
  A(:,:) = A(:, n:1:-1)

end function rotate_half_integer

end module MOM_array_transform
