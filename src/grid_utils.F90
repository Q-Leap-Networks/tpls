!> Utilities related process grids and domain grids and matching these.
!!
!! @author Mike Jackson.
!! @version $Revision: 328 $
!! @copyright (c) 2015, The University of Edinburgh, all rights
!! reserved. 
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module grid_utils

  implicit none

contains

  !> Calculate factors of n.
  !! Return an array of integers whose first nf entries are the
  !! factors found where:
  !! - factors(i, 1) and factors(i, 2) contain a pair of factors.
  !! - factors(i, 1) <= factors(i, 2).
  !! If n < 0 then behaviour is undefined.
  subroutine get_factors(n, nf, factors)
    integer, intent(in)  :: n  !< Number.
    integer, intent(out) :: nf !< Number of factors of n.
    integer, dimension(:,:), allocatable, intent(out) :: factors !< Factors.
    integer, dimension(:,:), allocatable :: local_factors
    integer :: i

    allocate(local_factors(n, 2))
    nf = 0
    do i = 1, int((n ** 0.5) + 1)
      if (mod(n, i) == 0) then
        nf = nf + 1
        local_factors(nf, 1) = n / i
        local_factors(nf, 2) = i
      end if
    end do
    ! Allocate an array of exactly the required size.
    allocate(factors(nf, 2))
    do i = 1, nf
      factors(i, 1) = local_factors(i,1)
      factors(i, 2) = local_factors(i,2)
    end do
    deallocate(local_factors)

  end subroutine get_factors


  !> Given an array of pairs of integers x,y, return an array of:
  !! - All pairs of integers x,y such that:
  !!   - x divides m.
  !!   - y divides n.
  !! - All pairs of integers y,x such that:
  !!   - y divides m.
  !!   - x divides n.
  subroutine get_divisors(m, n, np, pairs, nd, divisors)
    integer, intent(in) :: m  !< m of pair (m,n)
    integer, intent(in) :: n  !< n of pair (m,n)
    integer, intent(in) :: np !< Number of pairs.
    integer, dimension(:,:),intent(in) :: pairs !< Pairs.
    integer, intent(out) :: nd !< Number of divisor pairs.
    integer, dimension(:,:), allocatable, intent(out) :: divisors !< Divisor pairs.
    integer, dimension(:,:), allocatable :: local_divisors 
    integer :: i

    allocate(local_divisors(np * 2, 2))
    nd = 0
    do i = 1, np
      if (mod(m, pairs(i, 1)) == 0 .and. mod(n, pairs(i, 2)) == 0) then
        nd = nd + 1
        local_divisors(nd, 1) = pairs(i, 1)
        local_divisors(nd, 2) = pairs(i, 2)
      end if
      if (mod(m, pairs(i, 2)) == 0 .and. mod(n, pairs(i, 1)) == 0) then
        nd = nd + 1
        local_divisors(nd, 1) = pairs(i, 2)
        local_divisors(nd, 2) = pairs(i, 1)
      end if
    end do
    ! Allocate an array of exactly the required size.
    allocate(divisors(nd, 2))
    do i = 1, nd
      divisors(i, 1) = local_divisors(i,1)
      divisors(i, 2) = local_divisors(i,2)
    end do
    deallocate(local_divisors)

  end subroutine get_divisors


  !> Given an array of pairs of integers x,y, return the pair, x,y,
  !! whose ratio is closest to the ratio of m to n.
  subroutine get_closest_pair(m, n, np, pairs, x, y)
    integer, intent(in) :: m  !< m of pair (m,n)
    integer, intent(in) :: n  !< n of pair (m,n)
    integer, intent(in) :: np !< Number of pairs.
    integer, dimension(:,:), intent(in) :: pairs !< Pairs.
    integer, intent(out) :: x !< x of pair (x,y)
    integer, intent(out) :: y !< y of pair (x,y)
    double precision :: mton, xtoy, closest
    integer          :: i

    closest = 1000000
    x = 0
    y = 0
    mton = dble(m) / n
    do i = 1, np
      xtoy = dble(pairs(i, 1)) / pairs(i, 2)
      if (abs(xtoy - mton) <= closest) then
        closest = abs(xtoy - mton)
        x = pairs(i, 1)
        y = pairs(i, 2)
      end if
    end do
  end subroutine get_closest_pair


  !> Get a suitable 2D process grid for a 2D domain grid.
  !! Given a number of processes, np, and a 2D domain grid dimensions m,n
  !! return all suitable process grids x,y such that:
  !! - x * y == np.
  !! - x divides m.
  !! - y divides n.
  subroutine get_process_grids(np, m, n, ng, grids)
    integer, intent(in)  :: np !< Number of processes.
    integer, intent(in)  :: m  !< Domain grid X dimension.
    integer, intent(in)  :: n  !< Domain grid Y dimension.
    integer, intent(out) :: ng !< Number of grids.
    integer, dimension(:,:), allocatable, intent(out) :: grids !< Process grids.
    integer :: nf
    integer, dimension(:,:), allocatable :: factors

    call get_factors(np, nf, factors)
    call get_divisors(m, n, nf, factors, ng, grids)
    deallocate(factors)
  end subroutine get_process_grids


  !> Get a suitable 2D process grid for a 2D domain grid.
  !! Given a number of processes, np, and 2D domain grid dimensions m,n,
  !! calculate all suitable process grids x,y such that:
  !! - x * y == np.
  !! - x divides m.
  !! - y divides n.
  !! Out of all candidate x,y satisfying these criteria, return the
  !! x,y whose ratio closest matches that of m,n. 
  !! If there is no suitable candidate then 0,0 is returned.
  subroutine get_process_grid(np, m, n, x, y)
    integer, intent(in)  :: np !< Number of processes.
    integer, intent(in)  :: m  !< Domain grid X dimension.
    integer, intent(in)  :: n  !< Domain grid Y dimension.
    integer, intent(out) :: x  !< Process grid X dimensions.
    integer, intent(out) :: y  !< Process grid Y dimensions.
    integer :: nd
    integer, dimension(:,:), allocatable :: grids

    call get_process_grids(np, m, n, nd, grids)
    call get_closest_pair(m, n, nd, grids, x, y)
    deallocate(grids)
  end subroutine get_process_grid


  !> Given x and y return the first n pairs a,b in increasing order
  !! such that x divides a and y divides b.
  subroutine get_multiples(x, y, n, multiples)
    integer, intent(in) :: x !< x of pair(x,y)
    integer, intent(in) :: y !< y of pair(x,y)
    integer, intent(in) :: n !< Number of pairs to return.
    integer, dimension(n,2), intent(inout) :: multiples !< n pairs.
    integer :: i

    do i = 1, n
      multiples(i, 1) = x * i
      multiples(i, 2) = y * i
    end do
  end subroutine get_multiples

end module grid_utils
