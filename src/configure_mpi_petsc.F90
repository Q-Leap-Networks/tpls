module tpls_configure_mpi_petsc
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune,
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

  use petsc

  use pressure_solver, only : compute_matrix, compute_rhs, copy_pressure_in_to_petsc, calculate_indices
  use tpls_mpi
  use tpls_mpi_error_check

  implicit none

  public :: initialise
  public :: setup_MPI
  public :: get_size_and_rank
  public :: compute_initial_strides

#include "petsc/finclude/petscdef.h"

contains

  subroutine initialise(ierr)
    !-------------------------------------------------------------------!
    ! subroutine   I N I T I A L I S E                                  !
    !                                                                   !
    ! 1) Initialises petsc                                              !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Arguments:                                                        !
    !                                                                   !
    ! 1) ierr - error code produce by petsc                             !
    !-------------------------------------------------------------------!
    ! Necessary conditions:                                             !
    !                                                                   !
    ! 1) Petsc initialise cannot previously have been called            !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Toni Collis, V1.1, 30th May 2014                                  !
    !===================================================================!


    implicit none

    integer,          intent(inout) :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    if(ierr /= 0 ) then
       write(*,*) 'PetscInitialize failed. Error code: ', ierr
       stop
    end if

  end subroutine initialise
  

  subroutine setup_MPI(Ndim,dims,comm2d_quasiperiodic,&
       my_id,neighbours, &
       maxl, maxm, maxn, &
       global_dim_x, global_dim_y, global_dim_z, &
       sx,ex,sy,ey, &
       n_local_x,n_local_y,coords, &
       sz_uv, ez_uv, sz_w, ez_w, sz_p, ez_p, &
       n_local_z_uv, n_local_z_w, n_local_z_p, &
       petsc_y_ranges, ksp, da)
    !-------------------------------------------------------------------!
    ! subroutine   S E T U P __ M P I                                   !
    !                                                                   !
    ! 1) Call MPI comm size and rank                                    !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Arguments:                                                        !
    !                                                                   !
    ! 1) ierr - error code produce by petsc                             !
    ! 2) my_id - mpi process id, from comm_rank call                    !
    ! 2) num_procs - number of mpi processes, from comm_size call       !
    !-------------------------------------------------------------------!
    ! Necessary conditions:                                             !
    !                                                                   !
    ! 1) MPI must have been initialised (currently done by calling      !
    !  initialise subroutine that initialises MPI through PetSc         !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Toni Collis, V1.1, 30th May 2014                                  !
    !===================================================================!
    implicit none

    integer,                  intent(in)    :: Ndim
    integer, dimension(Ndim), intent(in)    :: dims
    integer,                  intent(in)    :: maxl, maxm, maxn
    integer,                  intent(in)    :: global_dim_x, global_dim_y, global_dim_z
    integer,                  intent(inout) :: comm2d_quasiperiodic
    integer,                  intent(inout) :: my_id
    integer, dimension(6),    intent(inout) :: neighbours
    integer,                  intent(inout) :: sx,sy,ex,ey
    integer,                  intent(inout) :: n_local_x,n_local_y
    integer, dimension(Ndim), intent(inout) :: coords
    integer,                  intent(inout) :: sz_uv, ez_uv, sz_w, ez_w, sz_p, ez_p
    integer,                  intent(inout) :: n_local_z_uv,n_local_z_w,n_local_z_p
   
    !Petsc stuff
    PetscInt, dimension(0:dims(1)-1), intent(inout)  :: petsc_y_ranges
    DM,                       intent(inout) :: da
    KSP,                      intent(inout) :: ksp
    ! The following are expected to be defined in 'pressure_solver'.
    ! PetscInt :: maxl, maxm, maxn

    !Local variables
    logical, dimension(Ndim) :: periodic(Ndim)
    logical                  :: reorder
    integer                  :: ierr
    PetscInt :: dof,stencil_width
    parameter (dof = 1, stencil_width = 1)
    PetscReal :: rtol, abstol, dtol
    PetscInt :: maxits
    KSPType :: solver_type
    PC :: pc
    PCType :: pc_type

    periodic(1) = .false.
    periodic(2) = .true.
    periodic(3) = .false.
    reorder = .false.

    Call mpi_cart_create(PETSC_COMM_WORLD,Ndim,dims,periodic,reorder,comm2d_quasiperiodic,ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_cart_create',my_id)
    Call mpi_comm_rank(  comm2d_quasiperiodic,my_id,ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_comm_rank',my_id)
    Call mpi_cart_coords(comm2d_quasiperiodic,my_id,Ndim,coords,ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_card_coords',my_id)
    Call get_mpi_neighbours(neighbours,comm2d_quasiperiodic)

    call mpi_decomp_2d(sx,ex,sy,ey,n_local_x,n_local_y,maxl,maxm,coords,dims,Ndim)

    petsc_y_ranges = n_local_x
    petsc_y_ranges(0) = n_local_x + 1
    petsc_y_ranges(dims(1)-1) = n_local_x + 1

    call DMDACreate3d(PETSC_COMM_WORLD,                               &
         DM_BOUNDARY_PERIODIC,                                      &
         DM_BOUNDARY_NONE,                                          &
         DM_BOUNDARY_NONE,                                          &
         DMDA_STENCIL_BOX,                                            &
         global_dim_x, global_dim_y+2, global_dim_z+2,                &
         dims(2),dims(1), dims(3),                                    &
         dof, stencil_width,                                          &
         PETSC_NULL_INTEGER,                                          &
         petsc_y_ranges,                                              &
         PETSC_NULL_INTEGER,                                          &
         da, ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call /DMDACreate3d',my_id)
    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPCreate',my_id)
    call KSPSetFromOptions(ksp, ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPSetFromOptions',my_id)
    call KSPSetComputeInitialGuess(ksp, copy_pressure_in_to_petsc, PETSC_NULL_OBJECT, ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPSetComputeInitialGuess',my_id)
    call KSPSetComputeRHS(ksp, compute_rhs, PETSC_NULL_OBJECT, ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPSetComputeRHS',my_id)
    call KSPSetComputeOperators(ksp, compute_matrix, PETSC_NULL_OBJECT, ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPSetComputeOperators',my_id)
    call KSPSetDM(ksp, da, ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPSetDM',my_id)

    call calculate_indices(da, sx, sy, ex, ey, n_local_x, n_local_y, my_id, ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call calculate_indices ',my_id)

    ! Defaults:
    !   rtol = 1.0d-5
    !   abstol = 1.0d-5
    !   dtol = 10000.0
    !   maxits = 10000
    call KSPSetTolerances(ksp,                             &
         PETSC_DEFAULT_REAL,  &
         PETSC_DEFAULT_REAL,  &
         PETSC_DEFAULT_REAL,  &
         PETSC_DEFAULT_INTEGER,           &
         ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPSetTolerances',my_id)

    if (is_master(my_id)) then
       call KSPGetType(ksp, solver_type, ierr)
       if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPGetType',my_id)
       call KSPGetPC(ksp, pc, ierr)
       if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPGetPC',my_id)
       call PCGetType(pc, pc_type, ierr)
       if(ierr /= 0 ) call abort(ierr, 'Failed to call PCGetType',my_id)
       write(*, *) 'Solver = ', solver_type
       write(*, *) 'Preconditioner = ', pc_type
       call KSPGetTolerances(ksp, rtol, abstol, dtol, maxits, ierr)
       if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPGetTolerances',my_id)
       write(*, *) 'rtol =', rtol, 'abstol =', abstol, 'dtol =', dtol, 'maxits =', maxits
       write(*, *) 'PERIODIC, NONE, NONE'
    end if

    sz_uv=0
    ez_uv=maxn-2
    n_local_z_uv=ez_uv-sz_uv+1        
    
    sz_w =0
    ez_w =maxn-1
    n_local_z_w =ez_w -sz_w +1 
    
    sz_p= 0
    ez_p= maxn
    n_local_z_p =ez_p- sz_p +1 
  
  end subroutine setup_MPI


  subroutine get_size_and_rank(ierr, my_id, num_procs)
    !-------------------------------------------------------------------!
    ! subroutine   G E T __ S I Z E __ A N D __ R A N K                 !
    !                                                                   !
    ! 1) Call MPI comm size and rank                                    !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Arguments:                                                        !
    !                                                                   !
    ! 1) ierr - error code produce by petsc                             !
    ! 2) my_id - mpi process id, from comm_rank call                    !
    ! 2) num_procs - number of mpi processes, from comm_size call       !
    !-------------------------------------------------------------------!
    ! Necessary conditions:                                             !
    !                                                                   !
    ! 1) MPI must have been initialised (currently done by calling      !
    !  initialise subroutine that initialises MPI through PetSc         !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Toni Collis, V1.1, 30th May 2014                                  !
    !===================================================================!


    implicit none

    integer,          intent(inout) :: ierr
    integer,          intent(inout) :: my_id
    integer,          intent(inout) :: num_procs

    call MPI_Comm_size(PETSC_COMM_WORLD,num_procs,ierr)
    if(ierr /= 0) call abort(ierr, 'Failed to call MPI_Comm_size',my_id)
    call MPI_Comm_rank(PETSC_COMM_WORLD,my_id,ierr)
    if(ierr /= 0) call abort(ierr, 'Failed to call MPI_Comm_rank',my_id)

  end subroutine get_size_and_rank


  subroutine compute_initial_strides(ex,ey,sx,sy,ex_max,my_id, &
       stride_uv_xz, stride_uv_yz,stride_w_xz,  stride_w_yz,stride_p_xz,  stride_p_yz, &
       sz_uv, ez_uv, sz_w, ez_w, sz_p, ez_p,stride_p_aug1_yz, stride_p_aug1_xz, &
       stride_p_aug2_yz, stride_p_aug2_xz,stride_p_augaug1_xz,stride_p_augaug1_yz, &
       stride_p_augaug2_xz,stride_p_augaug2_yz,stride_p_augaug3_xz,stride_p_augaug3_yz )
    !-------------------------------------------------------------------!
    ! subroutine C O M P U T E __ I N I T I A L __ S T R I D E S        !
    !                                                                   !
    ! 1)                              !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Arguments:                                                        !
    !                                                                   !
    ! 1)                            !
    !-------------------------------------------------------------------!
    ! Necessary conditions:                                             !
    !                                                                   !
    ! 1)      !
    !                                                                   !
    !-------------------------------------------------------------------!
    ! Toni Collis, V1.1, 30th May 2014                                  !
    !===================================================================!


    implicit none
    
    integer,          intent(in)    :: my_id
    integer,          intent(inout) :: ex,ey,sx,sy,ex_max
    integer,          intent(inout) :: stride_uv_xz, stride_uv_yz
    integer,          intent(inout) :: stride_w_xz,  stride_w_yz
    integer,          intent(inout) :: stride_p_xz,  stride_p_yz
    integer,          intent(inout) :: sz_uv, ez_uv, sz_w, ez_w, sz_p, ez_p
    integer,          intent(inout) :: stride_p_aug1_yz, stride_p_aug1_xz
    integer,          intent(inout) :: stride_p_aug2_yz, stride_p_aug2_xz
    integer,          intent(inout) :: stride_p_augaug1_xz,stride_p_augaug1_yz
    integer,          intent(inout) :: stride_p_augaug2_xz,stride_p_augaug2_yz
    integer,          intent(inout) :: stride_p_augaug3_xz,stride_p_augaug3_yz


    ! Local variables
    integer                         :: stride_uv_aug1_yz,stride_uv_aug1_xz
    integer                         :: stride_w_aug1_yz,stride_w_aug1_xz
    integer                         :: stride_uv_aug2_yz,stride_uv_aug2_xz
    integer                         :: stride_w_aug2_yz,stride_w_aug2_xz
    integer                         :: ierr

    ! Pick out max value ex -- signal for a right-hand boundary.
    Call mpi_allreduce(ex,ex_max,1,mpi_integer,mpi_max,PETSC_COMM_WORLD,ierr)
    if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_allreduce for max on ex',my_id)

    if(is_master(my_id))then
     write(*,*) '*******************************************'
     write(*,*) '** For channel BCs: ex_max=', ex_max
     write(*,*) '*******************************************'
    end if

    ! First set of strides - unchanged
    call get_stride_p(stride_uv_yz,stride_uv_xz,sx,ex,sy,ey,sz_uv,ez_uv)
    call get_stride_p(stride_w_yz, stride_w_xz, sx,ex,sy,ey,sz_w, ez_w )
    call get_stride_p(stride_p_yz,stride_p_xz,sx,ex,sy,ey,sz_p, ez_p)

    ! Second set of strides: first- and second-order halos for augmented arrays
    call get_stride_p_aug1(stride_uv_aug1_yz,stride_uv_aug1_xz,sx,ex,sy,ey,sz_uv,ez_uv)
    call get_stride_p_aug2(stride_uv_aug2_yz,stride_uv_aug2_xz,sx,ex,sy,ey,sz_uv,ez_uv)

    call get_stride_p_aug1(stride_w_aug1_yz,stride_w_aug1_xz,sx,ex,sy,ey,sz_w,ez_w)
    call get_stride_p_aug2(stride_w_aug2_yz,stride_w_aug2_xz,sx,ex,sy,ey,sz_w,ez_w)

    call get_stride_p_aug1(stride_p_aug1_yz,stride_p_aug1_xz,sx,ex,sy,ey,sz_p,ez_p)
    call get_stride_p_aug2(stride_p_aug2_yz,stride_p_aug2_xz,sx,ex,sy,ey,sz_p,ez_p)

    ! Third set of strides: first-, second-, and third-order halos for augmented arrays
    call get_stride_p_augaug1(stride_p_augaug1_yz,stride_p_augaug1_xz,sx,ex,sy,ey,sz_p,ez_p)
    call get_stride_p_augaug2(stride_p_augaug2_yz,stride_p_augaug2_xz,sx,ex,sy,ey,sz_p,ez_p)
    call get_stride_p_augaug2(stride_p_augaug3_yz,stride_p_augaug3_xz,sx,ex,sy,ey,sz_p,ez_p)

  end subroutine compute_initial_strides

end module tpls_configure_mpi_petsc
