!> Command-line tool to print information on possible process grids
!! and domain grids. 
!!
!! Usage: running as follows displays available options.
!! \code
!! Usage: grids <command> <arguments>
!! Commands:
!!   pgrids <N>
!!     Show all possible 2D process grids given N processes.
!!   ptodgrids <N> <L> <M>
!!     Given an LxM domain grid, show all 2D process grids given N processes
!!     where X divides L-1 and Y divides M-1.
!!   ptodgrid <N> <L> <M>
!!     Given an LxM domain grid, show the best 2D process grid given N processes
!!     where X divides L-1 and Y divides M-1.
!!     The best decomposition is that where the ratio X:Y is closest to the ratio L-1:M-1.
!!   dgrids <X> <Y> [<N>]
!!     Given an XxY process grid, show N valid domain grids LxM
!!     where X divides L-1 and Y divides M-1.
!!     If N is not provided then a default of 10 is used.
!!     This is not a complete set of all valid domains just a naively calculated subset.
!! \endcode
!!
!! @author Mike Jackson.
!! @version $Revision: 328 $
!! @copyright (c) 2015, The University of Edinburgh, all rights
!! reserved. 
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

program grids_program

  use grid_utils

  implicit none

  integer       :: num_arguments, np, nd, l, m, x, y, i, max
  character(40) :: command, argument
  integer, dimension(:,:), allocatable :: grid

  num_arguments = command_argument_count()
  call check_arguments(1, num_arguments)
  call get_command_argument(1, command)
  select case(command)
    case('pgrids')
      call check_arguments(2, num_arguments)
      call get_command_argument(2, argument)
      read(argument,*) np
      call get_factors(np, nd, grid)
      write(*,'(A,I6,A)') 'Process grids for ', np, ' processes:'
      do i = 1, nd
        write(*,'(I6,I6)') grid(i, 1), grid(i, 2)
      end do
      do i = 1, nd
        write(*,'(I6,I6)') grid(i, 2), grid(i, 1)
      end do
      deallocate(grid)
    case('ptodgrids')
      call check_arguments(4, num_arguments)
      call get_command_argument(2, argument)
      read(argument,*) np
      call get_command_argument(3, argument)
      read(argument,*) l
      call get_command_argument(4, argument)
      read(argument,*) m
      call get_process_grids(np, l - 1, m - 1, nd, grid)
      write(*,'(A,I6,A,I6,A,I6,A)') 'Process grids for placing a ', &
        l,' by ', m, ' domain grid onto ', np, ' processes:'
      do i = 1, nd
        write(*,'(I6,I6)') grid(i, 1), grid(i, 2)
      end do
      deallocate(grid)
    case('ptodgrid')
      call check_arguments(4, num_arguments)
      call get_command_argument(2, argument)
      read(argument,*) np
      call get_command_argument(3, argument)
      read(argument,*) l
      call get_command_argument(4, argument)
      read(argument,*) m
      call get_process_grid(np, l - 1, m - 1, x, y)
      write(*,'(A,I6,A,I6,A,I6,A)') 'Best process grid for placing a ', &
        l,' by ', m, ' domain grid onto ', np, ' processes:'
      write(*,'(I6,I6)') x, y
    case('dgrids')
      call check_arguments(3, num_arguments)
      call get_command_argument(2, argument)
      read(argument,*) l
      call get_command_argument(3, argument)
      read(argument,*) m
      max = 10
      if (num_arguments == 4) then
        call get_command_argument(4, argument)
        read(argument,*) max
      end if
      allocate(grid(max, 2))
      call get_multiples(l, m, max, grid)
      write(*,'(A,I6,A,I6,A)') 'Valid domain grids for a ', l, ' by ', m, ' process grid:'
      do i = 1, max
        write(*,'(I6,I6,A3,I6)') &
          grid(i, 1) + 1, grid(i, 2) + 1, ' = ', (grid(i, 1) + 1) * (grid(i, 2) + 1)
      end do
      deallocate(grid)
    case default
      call print_help()
  end select

end program grids_program

!> If actual number of arguments less than that expected then
!! print help and stop.
subroutine check_arguments(expected, actual)
  integer, intent(in) :: expected
  integer, intent(in) :: actual
 if (actual < expected) then
  call print_help()
  stop
 end if
end subroutine check_arguments

!> Print help.
subroutine print_help()
  write(*,*) 'Usage: grids <command> <arguments>'
  write(*,*) 'Commands:'
  write(*,*) '  pgrids <N>'
  write(*,*) '    Show all possible 2D process grids given N processes.'
  write(*,*) '  ptodgrids <N> <L> <M>'
  write(*,*) '    Given an LxM domain grid, show all 2D process grids given N processes'
  write(*,*) '    where X divides L-1 and Y divides M-1.'
  write(*,*) '  ptodgrid <N> <L> <M>'
  write(*,*) '    Given an LxM domain grid, show the best 2D process grid given N processes'
  write(*,*) '    where X divides L-1 and Y divides M-1.'
  write(*,*) '    The best decomposition is that where the ratio X:Y is closest to the ratio L-1:M-1.'
  write(*,*) '  dgrids <X> <Y> [<N>]'
  write(*,*) '    Given an XxY process grid, show N valid domain grids LxM'
  write(*,*) '    where X divides L-1 and Y divides M-1.'
  write(*,*) '    If N is not provided then a default of 10 is used.'
  write(*,*) '    This is not a complete set of all valid domains just a naively calculated subset.'
end subroutine print_help
