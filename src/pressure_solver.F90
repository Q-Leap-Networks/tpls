!> Pressure solver subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module pressure_solver

  use petsc

  use tpls_state

  implicit none

#include "petsc/finclude/petscdef.h"

  ! Included from tpls_state
  ! PetscInt :: maxl, maxm, maxn
  ! Note that the first two dimensions are swapped so that as
  ! far as PETSc is concerned  it is the first dimension that is
  ! periodic.  
  ! PetscInt :: global_dim_x, global_dim_y, global_dim_z
  ! PetscScalar, allocatable, dimension(:, :, :) :: pres
  ! PetscScalar :: dx, dy, dz, dt
  ! PetscScalar, allocatable, dimension(:, :, :) :: RHS_p, u3
  ! PetscScalar, allocatable, dimension(:) :: u_inlet

  PetscInt, private :: x_max,         y_max,         z_max
  PetscInt, private :: x_start,       y_start,       z_start
  PetscInt, private :: x_width,       y_width,       z_width
  PetscInt, private :: x_start_ghost, y_start_ghost, z_start_ghost
  PetscInt, private :: x_width_ghost, y_width_ghost, z_width_ghost
  PetscInt, private :: x_max_ghost,   y_max_ghost,   z_max_ghost

contains

  subroutine calculate_indices(da, sx, sy, ex, ey, n_local_x, n_local_y, rank, ierr)
    implicit none

    ! Arguments
    DM, intent(inout) :: da
    PetscInt, intent(in) :: sx, sy, ex, ey, n_local_x, n_local_y, rank
    PetscInt, intent(inout) :: ierr

    call DMDAGetCorners(da, x_start, y_start, z_start, x_width, y_width, z_width, ierr)
    x_max = x_start + x_width - 1
    y_max = y_start + y_width - 1
    z_max = z_start + z_width - 1

    if (sx==1 .AND. ((sx-1).NE.y_start .OR. ex/=y_max)) then
       write(*, *) 'ID1', rank, 'sx-1, y_start, ex, y_max, n_local_x, y_width', sx-1, y_start, ex, y_max, n_local_x, y_width
    end if
    if (ex==(maxl-1) .AND. (sx.NE.y_start .OR. (ex+1)/=y_max)) then
       write(*, *) 'ID2', rank, 'sx, y_start, ex+1, y_max, n_local_x, y_width', sx, y_start, ex+1, y_max, n_local_x, y_width
    end if
    if ((sx/=1 .AND. ex/=(maxl-1)) .AND. (sx/=y_start .OR. ex/=y_max)) then
       write(*, *) 'ID3', rank, 'sx, y_start, ex, y_max, n_local_x, y_width', sx, y_start, ex, y_max, n_local_x, y_width
    end if
    if (sy/=(x_start+1) .OR. ey/=(x_max+1)) then
       write(*, *) 'ID4', rank, 'sy, x_start, ey, x_max, n_local_y, x_width', sy, x_start, ey, x_max, n_local_y, x_width
    end if

    call DMDAGetGhostCorners(da, x_start_ghost, y_start_ghost, z_start_ghost,  &
         x_width_ghost, y_width_ghost, z_width_ghost, ierr)
    x_max_ghost = x_start_ghost + x_width_ghost - 1
    y_max_ghost = y_start_ghost + y_width_ghost - 1
    z_max_ghost = z_start_ghost + z_width_ghost - 1

    if ((sx-1)/=y_start_ghost .OR. (ex+1)/=y_max_ghost) then
       write(*, *) 'GHOST1', rank, 'sx-1, y_start_ghost, ex+1, y_max_ghost',  &
            sx-1, y_start_ghost, ex+1, y_max_ghost
    end if
    if ((sy-1)/=(x_start_ghost+1) .OR. (ey+1)/=(x_max_ghost+1)) then
       write(*, *) 'GHOST2', rank, 'sy-1, x_start_ghost+1, ey+1, x_max_ghost+1',  &
            sy-1, x_start_ghost+1, ey+1, x_max_ghost+1
    end if
    if (z_start/=z_start_ghost .OR. z_width/=z_width_ghost) then
       write(*, *) 'GHOST3', rank, 'z_start, z_start_ghost, z_width, z_width_ghost',  &
            z_start, z_start_ghost, z_width, z_width_ghost
    end if

  end subroutine calculate_indices


  subroutine output_converged_reason(ksp)

    implicit none

    ! Arguments
    KSP, intent(inout) :: ksp

    ! Local variables
    KSPConvergedReason :: reason
    character(29) :: converged_reason
    PetscInt :: ierr

    call KSPGetConvergedReason(ksp, reason, ierr)
    select case(reason)
    case (-10)
       converged_reason = 'KSP_DIVERGED_INDEFINITE_MAT'
    case (-9)
       converged_reason = 'KSP_DIVERGED_NAN'
    case (-8)
       converged_reason = 'KSP_DIVERGED_INDEFINITE_PC'
    case (-7)
       converged_reason = 'KSP_DIVERGED_NONSYMMETRIC'
    case (-6)
       converged_reason = 'KSP_DIVERGED_BREAKDOWN_BICG'
    case (-5)
       converged_reason = 'KSP_DIVERGED_BREAKDOWN'
    case (-4)
       converged_reason = 'KSP_DIVERGED_DTOL'
    case (-3)
       converged_reason = 'KSP_DIVERGED_ITS'
    case (-2)
       converged_reason = 'KSP_DIVERGED_NULL'
    case (-1)
       converged_reason = 'UNUSED'
    case (0)
       converged_reason = 'KSP_CONVERGED_ITERATING'
    case (1)
       converged_reason = 'KSP_CONVERGED_RTOL_NORMAL'
    case (2)
       converged_reason = 'KSP_CONVERGED_RTOL'
    case (3)
       converged_reason = 'KSP_CONVERGED_ATOL'
    case (4)
       converged_reason = 'KSP_CONVERGED_ITS'
    case (5)
       converged_reason = 'KSP_CONVERGED_CG_NEG_CURVE'
    case (6)
       converged_reason = 'KSP_CONVERGED_CG_CONSTRAINED'
    case (7)
       converged_reason = 'KSP_CONVERGED_STEP_LENGTH'
    case (8)
       converged_reason = 'KSP_CONVERGED_HAPPY_BREAKDOWN'
    case (9)
       converged_reason = 'KSP_CONVERGED_ATOL_NORMAL'
    end select
    write(*, *) 'reason =', converged_reason

  end subroutine output_converged_reason


  !> Copy the data generated by the other parts of TPLS into
  !! locations where they can be accessed by PETSc. 
  !! The global pressure array in the non-PETSc code is
  !! \code
  !! pres_global(0:maxl, 0:maxm, 0:maxn). 
  !! \endcode
  !! This includes a bounday layer on all sides even though the
  !!  second dimension is periodic as the periodic structure is
  !!  implemented by copying data one interior face of the cuboid  to
  !!  the opposite face in the boundary layer. 
  !! So the interior is
  !! \code
  !! (1:maxl-1, 1:maxm-1, 1:maxn-1)
  !! \endcode
  !! The indexing of the data in the PETSc code is different for two reasons.
  !! - PETSc lays out data on processes differently from the way
  !!   that it is done in the non-PETSc code which necessitates
  !!   swapping of the first two dimensions. 
  !! - PETSc is instructed to impose the periodic boundary condition
  !!  behind the scenes (through specifying DMDA_BOUNDARY_PERIODIC
  !!  for the appropriate dimesnion). The ghost points that are
  !!   required to do this are not visible in the global vector, but
  !!    they do appear in the local vectors. 
  !! Consequently the interior in the PETSc code is
  !! \code
  !! (0:maxm-2, 1:maxl -1, 1:maxn-1)
  !! \endcode
  !! which may written as
  !! \code
  !! (0:global_dim_x-1, 1:global_dim_y, 1:global_dim_z) 
  !! \endcode
  !! Including the explicit boundaries we have
  !! \code
  !! (0:global_dim_x-1, 0:global_dim_y+1, 0:global_dim_z+1) 
  !! \endcode
  !! or 
  !! \code
  !! (0:maxm-2, 0:maxl, 0:maxn) 
  !! \endcode
  !! Including the ghost points supplied by PETSc we have
  !! \code
  !! (-1:global_dim_x, 0:global_dim_y+1, 0:global_dim_z+1) 
  !! \endcode
  !! or 
  !! \code
  !! (-1:maxm-1, 0:maxl, 0:maxn)
  !! \endcode
  !! Let a pressure datum have coordinates (i1, j1, k1) in the non
  !! -PETSc code and (i2, j2, k2) in the PETSc code, then
  !! \code
  !! i2 = j1 - 1
  !! j2 = i1
  !! k2 = k1, or
  !! i1 = j2
  !! j1 = i2+1
  !! k1 = k2.
  !! \endcode
  subroutine copy_pressure_in_to_petsc(ksp, x, ctx, ierr)

    implicit none

    ! Arguments
    KSP, intent(inout) :: ksp
    Vec, intent(inout) :: x  ! x is a global vector.
    PetscErrorCode, intent(inout) ::  ierr
    PetscInt, intent(inout) :: ctx

    ! Local variables
    DM :: dm
    PetscScalar, pointer, dimension(:, :, :) :: x_3da
    PetscInt :: i, j, k
    call KSPGetDM(ksp,dm,ierr)
    call DMDAVecGetArrayF90(dm, x, x_3da, ierr)

    do k = z_start, z_max
       do j = y_start, y_max
          do i = x_start, x_max
             x_3da(i, j, k) = pres(j, i+1, k)
          end do
       end do
    end do

    call DMDAVecRestoreArrayF90(dm, x, x_3da, ierr)
  end subroutine copy_pressure_in_to_petsc


  subroutine compute_rhs(ksp, b, ctx, ierr)

    implicit none

    ! Arguments
    KSP, intent(inout) :: ksp
    Vec, intent(inout) :: b  ! b is a global vector.
    PetscInt, intent(inout) :: ctx
    PetscErrorCode, intent(inout) :: ierr

    ! Local variables
    PetscScalar, pointer, dimension(:, :, :) :: b_3da
    PetscInt :: i, j, k
    DM :: dm

    call KSPGetDM(ksp, dm, ierr)

    call DMDAVecGetArrayF90(dm, b, b_3da, ierr)

    ! One could introduce a new variable, b_rhs, shaped like RHS_p
    !  into main_ns_hybrid.F90. 
    !          b_rhs = -RHS_p
    !          if (sx==1) then
    !            do k = 1, maxn-1
    !              do j = sy, ey
    !                b_rhs(0, j, k)=(dx/dt)*(u_inlet(k-1)-u3(0, j-1, k-1))/dy
    !              end do
    !            end do
    !          end if
    ! and then use it here thus
    !        do k = z_start, z_max
    !          do j = y_start, y_max
    !            do i = x_start, x_max
    !              b_3da(i, j, k) = b_rhs(j, i+1, k)
    !            end do
    !          end do
    !        end do
    ! This may be easier to read but the following uses less space.
    do k = z_start, z_max
       do j = y_start, y_max
          do i = x_start, x_max
             b_3da(i, j, k) = -RHS_p(j, i+1, k)
          end do
       end do
    end do

    if (y_start==0) then
       do k = 1, global_dim_z
          do i = x_start, x_max
             b_3da(i, 0, k) = (dx/dt)*(u_inlet(k-1)-u3(0, i ,k-1))/dy
          end do
       end do
    end if

    call DMDAVecRestoreArrayF90(dm, b, b_3da, ierr)

    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

  end subroutine compute_rhs


  subroutine compute_matrix(ksp, A, B, ctx, ierr)

    implicit none

    ! Arguments
    KSP, intent(inout) :: ksp
    Mat, intent(inout) :: A, B
    PetscInt, intent(inout) :: ctx
    PetscErrorCode, intent(inout) :: ierr

    ! Local variables
    PetscInt :: i, j, k
    PetscScalar    :: v(7)
    MatStencil     :: row(4), col(4, 7)
    MatNullSpace   :: nullspace

    do k = z_start, z_max
       do j = y_start, y_max
          do i = x_start, x_max
             row(MatStencil_i) = i
             row(MatStencil_j) = j
             row(MatStencil_k) = k
             ! Deal with the edges of the cuboid.
             ! The edges with i=0 and i=(global_dim_x-1) are taken
             !  care of by the periodic boundary condition. 
             if ((j==0 .AND. k==0) .OR.                  &
                  (j==(global_dim_y+1) .AND. k==0)) then
                v(1) = 1/dz
                col(MatStencil_i, 1) = i
                col(MatStencil_j, 1) = j
                col(MatStencil_k, 1) = k
                v(2) = -1/dz
                col(MatStencil_i, 2) = i
                col(MatStencil_j, 2) = j
                col(MatStencil_k, 2) = k+1
                call MatSetValuesStencil(B, 1, row, 2, col, v, INSERT_VALUES, ierr)
             else if ((j==0 .AND. k==(global_dim_z+1)) .OR.                  &
                  (j==(global_dim_y+1) .AND. k==(global_dim_z+1))) then
                v(1) = -1/dz
                col(MatStencil_i, 1) = i
                col(MatStencil_j, 1) = j
                col(MatStencil_k, 1) = k
                v(2) = 1/dz
                col(MatStencil_i, 2) = i
                col(MatStencil_j, 2) = j
                col(MatStencil_k, 2) = k-1
                call MatSetValuesStencil(B, 1, row, 2, col, v, INSERT_VALUES, ierr)
                ! Deal with the faces of the cuboid, excluding the
                !  edges and the corners. 
                ! The faces i=0 and i=(global_dim_x-1) are taken care
                !  of by the periodic boundary condition. 
             else if (j==0) then
                ! Von Neumann boundary condition on y=0 boundary.
                v(1) = 1/dy
                col(MatStencil_i, 1) = i
                col(MatStencil_j, 1) = j
                col(MatStencil_k, 1) = k
                v(2) = -1/dy
                col(MatStencil_i, 2) = i
                col(MatStencil_j, 2) = j+1
                col(MatStencil_k, 2) = k
                call MatSetValuesStencil(B, 1, row, 2, col, v, INSERT_VALUES, ierr)
             else if (j==(global_dim_y+1)) then
                ! Von Neumann boundary condition on y=(global_dim_y+1) boundary.
                v(1) = -1/dy
                col(MatStencil_i, 1) = i
                col(MatStencil_j, 1) = j
                col(MatStencil_k, 1) = k
                v(2) = 1/dy
                col(MatStencil_i, 2) = i
                col(MatStencil_j, 2) = j-1
                col(MatStencil_k, 2) = k
                call MatSetValuesStencil(B, 1, row, 2, col, v, INSERT_VALUES, ierr)
             else if (k==0) then
                ! Von Neumann boundary condition on z=0 boundary.
                v(1) = 1/dz
                col(MatStencil_i, 1) = i
                col(MatStencil_j, 1) = j
                col(MatStencil_k, 1) = k
                v(2) = -1/dz
                col(MatStencil_i, 2) = i
                col(MatStencil_j, 2) = j
                col(MatStencil_k, 2) = k+1
                call MatSetValuesStencil(B, 1, row, 2, col, v, INSERT_VALUES, ierr)
             else if (k==(global_dim_z+1)) then
                ! Von Neumann boundary conditions on z=(global_dim_z+1) boundary.
                v(1) = -1/dz
                col(MatStencil_i, 1) = i
                col(MatStencil_j, 1) = j
                col(MatStencil_k, 1) = k
                v(2) = 1/dz
                col(MatStencil_i, 2) = i
                col(MatStencil_j, 2) = j
                col(MatStencil_k, 2) = k-1
                call MatSetValuesStencil(B, 1, row, 2, col, v, INSERT_VALUES, ierr)
             else
                ! Deal with the interior.
                ! Laplacian in 3D.
                v(1) = -1/dz**2
                col(MatStencil_i, 1) = i
                col(MatStencil_j, 1) = j
                col(MatStencil_k, 1) = k-1
                v(2) = -1/dy**2
                col(MatStencil_i, 2) = i
                col(MatStencil_j, 2) = j-1
                col(MatStencil_k, 2) = k
                v(3) = -1/dx**2
                col(MatStencil_i, 3) = i-1
                col(MatStencil_j, 3) = j
                col(MatStencil_k, 3) = k
                v(4) = 2*(1/dx**2+1/dy**2+1/dz**2)
                col(MatStencil_i, 4) = i
                col(MatStencil_j, 4) = j
                col(MatStencil_k, 4) = k
                v(5) = -1/dx**2
                col(MatStencil_i, 5) = i+1
                col(MatStencil_j, 5) = j
                col(MatStencil_k, 5) = k
                v(6) = -1/dy**2
                col(MatStencil_i, 6) = i
                col(MatStencil_j, 6) = j+1
                col(MatStencil_k, 6) = k
                v(7) = -1/dz**2
                col(MatStencil_i, 7) = i
                col(MatStencil_j, 7) = j
                col(MatStencil_k, 7) = k+1
                call MatSetValuesStencil(B, 1, row, 7, col, v, INSERT_VALUES, ierr)
             end if
          end do
       end do
    end do

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY, ierr)

    if (A.ne.B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
       call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    endif

    call MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, (/PETSC_NULL_OBJECT/),  &
         nullspace, ierr)
    call MatSetNullSpace(A, nullspace, ierr)
    call MatNullSpaceDestroy(nullspace, ierr)

  end subroutine compute_matrix


  subroutine copy_pressure_out_of_petsc(ksp, ierr)

    implicit none

    ! Arguments
    KSP, intent(in) :: ksp
    PetscErrorCode, intent(out) ::  ierr

    ! Local variables
    DM :: dm

    Vec :: global_vector
    Vec :: local_vector

    PetscScalar, pointer, dimension(:, :, :) :: local_vector_3da
    PetscScalar, pointer, dimension(:, :, :) :: global_vector_3da

    PetscInt :: i, j, k

    call KSPGetDM(ksp, dm, ierr)
    call KSPGetSolution(ksp, global_vector, ierr)
    call DMCreateLocalVector(dm, local_vector, ierr)
    call DMGlobalToLocalBegin(dm, global_vector, INSERT_VALUES, local_vector, ierr)
    call DMGlobalToLocalEnd(dm, global_vector, INSERT_VALUES, local_vector, ierr)

    call DMDAVecGetArrayF90(dm, local_vector, local_vector_3da, ierr)
    call DMDAVecGetArrayF90(dm, global_vector, global_vector_3da, ierr)

    do k = z_start_ghost, z_max_ghost
       do j = y_start_ghost, y_max_ghost
          do i = x_start_ghost, x_max_ghost
             pres(j, i+1, k) = local_vector_3da(i, j, k)
          end do
       end do
    end do

    call DMDAVecRestoreArrayF90(dm, global_vector, global_vector_3da, ierr)
    call DMDAVecRestoreArrayF90(dm, local_vector, local_vector_3da, ierr)

    ! As call VecDestroy(global_vector, ierr) causes a memory error,
    ! the following call is required to prevent a memory leak.
    call VecDestroy(local_vector, ierr)
  end subroutine copy_pressure_out_of_petsc

end module pressure_solver
