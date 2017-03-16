!> TPLS global data and state.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune,
!! Toni Collis, David Scott, Peter Spelt, Mike Jackson.
!! @version $Revision: 252 $
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt, The University of
!! Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module tpls_state

#include "petsc/finclude/petscdef.h"

  use petsc 

  implicit none

  !> Number of domain grid points in X(l) dimension.
  PetscInt :: maxl
  !> Number of domain grid points in Y(m) dimension.
  PetscInt :: maxm
  !> Number of domain grid points in Z(n) dimension.
  PetscInt :: maxn
  !> Number of domain grid points in Y(m) dimension - 1.
  !! Note that the X(l) and Y(m) dimensions are swapped so that as
  !! far as PETSc is concerned  it is the first dimension that is 
  !! periodic.  
  PetscInt :: global_dim_x
  !> Number of domain grid points in X(l) dimension - 1.
  !! Note that the X(l) and Y(m) dimensions are swapped so that as
  !! far as PETSc is concerned  it is the first dimension that is 
  !! periodic.  
  PetscInt :: global_dim_y
  !> Number of domain grid points in Z(n) dimension - 1.
  PetscInt :: global_dim_z
  !> Pressure.
  PetscScalar, allocatable, dimension(:, :, :) :: pres
  !> dx.
  PetscScalar :: dx
  !> dy.
  PetscScalar :: dy
  !> dz.
  PetscScalar :: dz
  !> dt.
  PetscScalar :: dt
  !> Velocity in X(l) dimension.
  PetscScalar, allocatable, dimension(:, :, :) :: u3
  !> RHS pressure.
  PetscScalar, allocatable, dimension(:, :, :) :: RHS_p
  !> Velocity in X(l) dimension at inlet.
  PetscScalar, allocatable, dimension(:) :: u_inlet

end module tpls_state
