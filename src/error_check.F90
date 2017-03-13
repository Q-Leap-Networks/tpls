!> Error checking subroutines.
!!
!! @author Antonia Collis, Mike Jackson.
!! @version $Revision: 252 $
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt, The University of
!! Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module tpls_error_check

  implicit none

contains

  !> Print an error code and message and abort.
  !! Error printed to stderr in the format:
  !! \code
  !! message
  !! error=err
  !! \endcode
  subroutine abort(err, message)
    implicit none
    integer,          intent(in) :: err
    character(len=*), intent(in) :: message
    logical                      :: open_stderr
    integer                      :: stderr, stdout
    parameter(stderr=0,stdout=6)

    ! Flush so error messages don't get out of synch with
    ! outputs if both are redirected to same destination.
    flush stdout
    flush stderr
    open_stderr=.true.
    ! Open the stderr file, but not if the message is blank
    if(len_trim(message)==0)open_stderr=.false.
    if(open_stderr)then
      write(stderr,'(A)') trim(message)
      write(stderr,'(A6,I3)') "error=",err
    end if
    stop
  end subroutine abort

end module tpls_error_check
