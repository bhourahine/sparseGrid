!> Some computational constants
module accuracy
  implicit none
  public

  !> double precision
  integer, parameter :: dp = kind(0.0d0)

end module accuracy

!> Simple exampls of Smolyak sparse grids
module Smolyak
  use accuracy
  implicit none

  public :: enumerateGrids

contains

  !> Enumerate over all tensor product grids of orders sum(o_i) <= nLevel in nDims
  subroutine enumerateGrids(nLevel, nDim)

    !> Level of sparse grid
    integer, intent(in) :: nLevel

    !> number of spatial dimensions
    integer, intent(in) :: nDim

    !> loop counter
    integer :: counter(nDim)

    !> current dimention in the loop
    integer :: iDim

    counter(:) = 0
    outer_loop : do

      ! Generate that particular product of grids
      call evaluateGrid(counter)

      ! loop the index over which counter this is
      iDim = 1

      ! Calculate the next indices:
      inner_loop : do

        counter(iDim) = counter(iDim) + 1

        ! If this is still a valid index, exit the inner loop and go for the next iteration
        if (sum(counter) <= nLEvel) exit inner_loop

        ! This index has overshot, so reset the counter and move on to the next index.
        counter(iDim) = 0

        iDim = iDim + 1

        ! If the next index would be outside of nDim, the whole loop is finished.
        if (iDim > nDim) exit outer_loop

      end do inner_loop

    end do outer_loop

  end subroutine enumerateGrids

  ! Evaluate a particular tensor product of grids
  subroutine evaluateGrid(orders)

    !> Order of the individual grids
    integer, intent(in) :: orders(:)

    !> Number of dimensions
    integer :: nDim

    !> Dimensional index
    integer :: iDim

    !> Loop counters
    integer, allocatable :: counter(:)

    !> Resulting point in space
    real(dp), allocatable :: x(:)

    nDim = size(orders)

    allocate(counter(nDim))
    allocate(x(nDim))

    counter(:) = 1
    outer_loop : do

      call gridPoint(x,counter,orders)
      write(45,*)x

      iDim = 1
      ! Calculate the next indices:
      inner_loop : do
        counter(iDim) = counter(iDim) + 1

        ! If this is still a valid index, exit the inner loop and go for the next iteration
        if (counter(iDim) <= nPoints(orders(iDim))) exit inner_loop

        ! This index has overshot, so reset the counter and move to next index.
        counter(iDim) = 1

        iDim = iDim + 1

        ! If the next index would be outside of nDim, the whole loop is finished.
        if (iDim > size(orders)) exit outer_loop

      end do inner_loop
    end do outer_loop

  end subroutine evaluateGrid

  ! returns a point in the [0..1]^n space
  subroutine gridPoint(point, localIndex, level)

    !> resulting point in space
    real(dp), intent(out) :: point(:)

    !> point on each of the 1D grids
    integer, intent(in) :: localIndex(:)

    !> 1D grid levels for each dimension
    integer, intent(in) :: level(:)

    integer :: ii, n

    n = size(level)
    do ii = 1, n
      point(ii) = location(localIndex(ii), level(ii))
    end do

  end subroutine gridPoint

  ! Internal function that generates a point on a 1D grid for that level of grid refinement
  function location(iPt,nLevel)

    !> which point, from 1 .. nPoints
    integer, intent(in) :: iPt

    !> Level of grid refinement
    integer, intent(in) :: nLevel

    !> resulting value of grid
    real(dp) :: location

    location = 1.0 / 2**nLevel + real(iPt-1,dp) / max(1,2**(nLevel-1))
    if (nLevel == 0) then
      location = location -1.0_dp
    end if

  end function location

  !> Evaluate the number of points at the given 1D grid level
  function nPoints(nLevel)

    !> grid level
    integer, intent(in) :: nLevel

    !> total number of points
    integer :: nPoints

    select case(nLevel)
    case(0)
      nPoints = 2
    case(1)
      nPoints = 1
    case default
      nPoints = 2**(nLevel-1)
    end select

  end function nPoints

end module Smolyak

!> Simple input driver
program testSmolyakGrid
  use accuracy
  use smolyak
  implicit none

  integer :: nDim, nLevel

  write(*,*)'Level of refinement, Dimensionality'
  read(*,*)nLEvel, nDim

  call enumerateGrids(nLevel, nDim)

end program testSmolyakGrid
