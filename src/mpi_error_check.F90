!> TPLS_MPI_error_check subroutines.
!!
!! @author Toni Collis
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module tpls_mpi_error_check
    !-------------------------------------------------------------------!
    ! module T P L S _ M P I _ E R R O R _ C H E C K                    !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Toni Collis, TPLS V1.1, 6th May 2014                              !
    !===================================================================!

  use petsc
  
  use tpls_mpi

  implicit none

  public abort

contains

  subroutine abort(err, message, my_id)
    !-------------------------------------------------------------------!
    ! subroutine   A B O R T                                            !
    !                                                                   !
    ! 1) Writes 'message' to 'stderr'                                   !
    ! 2) Closes all open units                                          !
    ! 3) Aborts MPI                                                     !
    ! 4) Exits                                                          !
    !                                                                   !
    ! Note: if "message" is a null string then no stderr file is opened !
    !-------------------------------------------------------------------!
    ! Arguments:                                                        !
    !                                                                   !
    ! 1) message : input : error message written on abort               !
    !-------------------------------------------------------------------!
    ! Necessary conditions:                                             !
    !                                                                   !
    ! 1) File unit=0 must not have been used                            !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Toni Collis, V1.1, 6th May 2014                                   !
    !===================================================================!


    implicit none

    character(len=*), intent(in) :: message
    integer,          intent(in) :: err
    integer,          intent(in) :: my_id
    integer                      :: unit
    logical                      :: open_stderr
    integer                      :: stderr, stdout, mpi_err
    parameter(stderr=0,stdout=6)

    ! Flush stdout so error messages don't get out of synch with
    ! outputs if both are redirected to same destination.
    flush stdout

    open_stderr=.true.
    ! Open the stderr file, but not if the message is blank
    if(len_trim(message)==0)open_stderr=.false.
    if(open_stderr)then
      ! ** Write "message" to "stderr"
       if(is_master(my_id)) then
          write(stderr,'(a)') trim(message)
          write(stderr,'(A6,I3)') "error=",err
       end if
       write(stderr,'(A13,I6,A9)') "TPLS process ", my_id, " aborting"
    end if

    ! Close all units
    !TC to do - check if there is a specific list of units to be used in TPLS
    do unit=9,9
       close(unit)
    end do

    !Abort execution
    !Terminate MPI execution environment
    call MPI_Abort(PETSC_COMM_WORLD,1,mpi_err)
    !Terminate 
    stop

  end subroutine abort

end module tpls_mpi_error_check
