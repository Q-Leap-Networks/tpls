!> Cahn-Hilliard (DIM) subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

! David Scott:
! In deriving this code, the following correspondences between the
!  variables of the initial DIM code and those of the TPLS code, have
!  been assumed. 
! fi(i, j, k) <-> phi(i-1, j-1, k-1)
! fiold(i, j, k) <-> phiOld(i-1, j-1, k-1)
! finew(i, j, k) <-> phiNew(i-1, j-1, k-1)
! d2fi(i, j, k) <-> d2phi(i-1, j-1, k-1)
! RHS(i, j, k) <-> RHS0(i-1, j-1, k-1) (renamed conv2_phi)
! source(i, j, k) <-> source0(i-1, j-1, k-1)
! bulkE(i, j, k) <-> bulkE0(i-1, j-1, k-1)
! dfx(i, j, k) <-> dfx0(i-1, j-1, k-1)
! dfy(i, j, k) <-> dfy0(i, j-1, k-1)
! dfz(i, j, k) <-> dfz0(i, j-1, k-1)
! u(i, j, k) <-> u0(i-1, j-1, k-1)
! v(i, j, k) <-> v0(i, j-1, k-1)
! w(i, j, k) <-> w0(i, j-1, k-1)
! ConvPrevC(i, j, k) <-> ConvPrevPhi(i-1, j-1, k-1)

module cahn_hilliard_solver

  use tpls_mpi

  implicit none

contains

  !> Calculate DIM.
  !! IMPORTANT: VERSION 1 according to notes (Dimensional Analysis
  !!  summary and CH notes) 
  !! Note that sz_p = 0 and ez_p = maxn.
  subroutine dim(maxl, maxn, dx, dy, dz, dt, epn, Pe, phi, phiOld, conv2_phi, &
       ConvPrevPhi, phiNew, sx, ex, sy, ey, sz_p, ez_p, height,               &
       stride_p_xz, stride_p_yz, neighbours, comm2d_quasiperiodic,            &
       max_iteration_dim)

    implicit none

    ! Arguments
    integer, intent(in) :: maxl, maxn
    integer, intent(in) :: sx, ex, sy, ey, sz_p, ez_p
    integer, intent(in) :: stride_p_xz, stride_p_yz
    integer, intent(in) :: neighbours(6)
    integer, intent(in) :: comm2d_quasiperiodic
    integer, intent(in) :: max_iteration_dim

    double precision, intent(in) :: dx, dy, dz, dt
    double precision, intent(in) :: Pe, epn
    double precision, intent(in) :: height
    double precision, intent(in) :: phi(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)
    double precision, intent(in) :: phiOld(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)
    double precision, intent(in) :: conv2_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    double precision, intent(inout) :: ConvPrevPhi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    double precision, intent(inout) :: phiNew(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)

    ! Local variables
    integer :: i, j, k
    integer :: iteration

    double precision :: bulkE0(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)
    double precision :: d2phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    double precision :: source0(sx:ex,sy:ey,sz_p:ez_p)

    double precision, parameter :: nmda = 0.25d0
    double precision, parameter :: tao = 0.25d0

    double precision :: c1
    double precision :: a1, a2, a3, a4, a5, a6
    double precision :: temp1, temp2, temp3, temp4, temp5, rhoB

    double precision :: some
    double precision :: diagonal
    double precision :: epn2
    double precision, parameter :: alfa1 = 1.1d0

    integer :: num_threads
    integer :: OMP_GET_NUM_THREADS
    integer :: tid
    integer :: OMP_GET_THREAD_NUM

    double precision :: z_val
    double precision :: res
    double precision, allocatable :: res_vec(:)
    integer :: alloc_stat

    !$OMP PARALLEL SHARED(num_threads)
    !$OMP MASTER
    num_threads = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
    !$OMP END PARALLEL

    allocate(res_vec(0:num_threads-1), STAT=alloc_stat)
    if (alloc_stat.ne.0) STOP '** Not enough memory **'

    ! Calculations.

    epn2 = epn*epn

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j, k)
    do i = sx, ex
       do j = sy, ey
          do k = 1, maxn-1

             bulkE0(i, j, k) = (1.d0 - 3.d0*phi(i, j, k) + 2.d0*phi(i, j, k)*phi(i, j, k))*phi(i, j, k)/2.d0

             d2phi(i, j, k) = (phi(i+1, j, k) - 2.d0*phi(i, j, k) + phi(i-1, j, k))/dx/dx &
                  + (phi(i, j+1, k) - 2.d0*phi(i, j, k) + phi(i, j-1, k))/dy/dy &
                  + (phi(i, j, k+1) - 2.d0*phi(i, j, k) + phi(i, j, k-1))/dz/dz

             !source0(i, j, k) = -ConvPrevPhi(i, j, k)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO


    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j)
    do i = sx, ex
       do j = sy, ey
          ! z direction
          ! non-periodic
          bulkE0(i, j, 0) = bulkE0(i, j, 1)
          bulkE0(i, j, maxn) = bulkE0(i, j, maxn-1)
          d2phi(i, j, 0) = d2phi(i, j, 1)
          d2phi(i, j, maxn) = d2phi(i, j, maxn-1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! Exchange first order halos.
    call exchange2d(bulkE0, stride_p_xz, stride_p_yz, neighbours,     &
         ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)
    call exchange2d(d2phi, stride_p_xz, stride_p_yz, neighbours,      &
         ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)

    ! TODO - WHAT SHOULD THE BOUNDARY CONDITIONS BE AT THE INLET?
    if (sx==1) then
       !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j, k)
       do j = sy, ey
          do k = 1, maxn-1
             ! x direction
             ! non-periodic
             d2phi(0, j, k) = d2phi(1, j, k)
             bulkE0(0, j, k) = bulkE0(1, j, k)
          enddo
       enddo
       !$OMP END PARALLEL DO
    endif
    if (ex==(maxl-1)) then
       !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j, k)
       do j = sy, ey
          do k = 1, maxn-1
             ! x direction
             ! non-periodic
             d2phi(maxl, j, k) = d2phi(maxl-1, j, k)
             bulkE0(maxl, j, k) = bulkE0(maxl-1, j, k)
          enddo
       enddo
       !$OMP END PARALLEL DO
    endif

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j, k, c1, a1, a2, a3, a4, a5, a6, &
    !$OMP& temp1, temp2, temp3, temp4, temp5, rhoB)
    do i = sx, ex
       do j = sy, ey
          do k = 1, maxn-1
             c1 = (phi(i-1, j, k) + phi(i, j, k))/2.d0
             a1 = dsqrt(c1*c1*(1.d0 - dabs(c1))*(1.d0 - dabs(c1)))
             c1 = (phi(i+1, j, k) + phi(i, j, k))/2.d0
             a2 = dsqrt(c1*c1*(1.d0 - dabs(c1))*(1.d0 - dabs(c1)))
             c1 = (phi(i, j-1, k) + phi(i, j, k))/2.d0
             a3 = dsqrt(c1*c1*(1.d0 - dabs(c1))*(1.d0 - dabs(c1)))
             c1 = (phi(i, j+1, k) + phi(i, j, k))/2.d0
             a4 = dsqrt(c1*c1*(1.d0 - dabs(c1))*(1.d0 - dabs(c1)))
             if (k+1.eq.2) then
                c1 = phi(i, j, k)
             else
                c1 = (phi(i, j, k-1) + phi(i, j, k))/2.d0
             endif
             a5 = dsqrt(c1*c1*(1.d0 - dabs(c1))*(1.d0 - dabs(c1)))
             c1 = (phi(i, j, k+1) + phi(i, j, k))/2.d0
             a6 = dsqrt(c1*c1*(1.d0 - dabs(c1))*(1.d0 - dabs(c1)))

             temp1 = (a2*(bulkE0(i+1, j, k) - bulkE0(i, j, k)) &
                  - a1*(bulkE0(i, j, k) - bulkE0(i-1, j, k)))/dx/dx &
                  + (a4*(bulkE0(i, j+1, k) - bulkE0(i, j, k)) &
                  - a3*(bulkE0(i, j, k) - bulkE0(i, j-1, k)))/dy/dy &
                  + (a6*(bulkE0(i, j, k+1) - bulkE0(i, j, k)) &
                  - a5*(bulkE0(i, j, k) - bulkE0(i, j, k-1)))/dz/dz

             temp2 = (d2phi(i+1, j, k) - 2.d0*d2phi(i, j, k) + d2phi(i-1, j, k))/dx/dx &
                  + (d2phi(i, j+1, k) - 2.d0*d2phi(i, j, k) + d2phi(i, j-1, k))/dy/dy &
                  + (d2phi(i, j, k+1) - 2.d0*d2phi(i, j, k) + d2phi(i, j, k-1))/dz/dz

             temp3 = (a2*(d2phi(i+1, j, k) - d2phi(i, j, k)) &
                  - a1*(d2phi(i, j, k) - d2phi(i-1, j, k)))/dx/dx &
                  + (a4*(d2phi(i, j+1, k) - d2phi(i, j, k)) &
                  - a3*(d2phi(i, j, k) - d2phi(i, j-1, k)))/dy/dy &
                  + (a6*(d2phi(i, j, k+1) - d2phi(i, j, k)) &
                  - a5*(d2phi(i, j, k) - d2phi(i, j, k-1)))/dz/dz

             temp4 = (temp1 - epn2*temp3)/epn - nmda*tao/2.d0*d2phi(i, j, k) + epn2*nmda/2.d0*temp2

             temp5 = temp4/Pe - conv2_phi(i, j, k)

             !              ConvPrevPhi(i, j, k) = temp5
             !
             !              source0(i, j, k) = source0(i, j, k) + temp5*2.d0
             source0(i, j, k) = temp5*2.d0 - ConvPrevPhi(i, j, k)
             ConvPrevPhi(i, j, k) = temp5

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO


    a1 = 1.d0/(dx**2.d0) + 1.d0/(dy**2.d0) + 1.d0/(dz**2.d0)
    a2 = 2.d0/(dx**4.d0) + 2.d0/(dy**4.d0) + 2.d0/(dz**4.d0)
    diagonal = 1.5d0 + nmda/2.d0/Pe*dt*(2.d0*tao*a1 + epn2*(4.d0*a1*a1+a2))

    ! SOR.
    do iteration = 1, max_iteration_dim

       res = 0.d0

       !$OMP PARALLEL WORKSHARE
       forall(i = sx:ex, j = sy:ey, k = 1:maxn-1)

          ! d2phi: grad2(C)
          d2phi(i, j, k) = (phiNew(i+1, j, k) - 2.d0*phiNew(i, j, k) + phiNew(i-1, j, k))/dx/dx &
               + (phiNew(i, j+1, k) - 2.d0*phiNew(i, j, k) + phiNew(i, j-1, k))/dy/dy &
               + (phiNew(i, j, k+1) - 2.d0*phiNew(i, j, k) + phiNew(i, j, k-1))/dz/dz

       endforall
       !$OMP END PARALLEL WORKSHARE


       !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j)
       do i = sx, ex
          do j = sy, ey
             ! Non-Periodic in z
             d2phi(i, j, 0) = d2phi(i, j, 1)
             d2phi(i, j, maxn) = d2phi(i, j, maxn-1)
          enddo
       enddo
       !$OMP END PARALLEL DO

       ! Exchange first order halos.
       call exchange2d(d2phi, stride_p_xz, stride_p_yz, neighbours,     &
            ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)

       ! TODO - WHAT SHOULD THE BOUNDARY CONDITIONS BE AT THE INLET?
       if (sx==1) then
          !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j, k)
          do j = sy, ey
             do k = 1, maxn-1
                ! Non-Periodic in x
                d2phi(0, j, k) = d2phi(1, j, k)
             enddo
          enddo
          !$OMP END PARALLEL DO
       endif
       if (ex==(maxl-1)) then
          !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j, k)
          do j = sy, ey
             do k = 1, maxn-1
                ! Non-Periodic in x
                d2phi(maxl, j, k) = d2phi(maxl-1, j, k)
             enddo
          enddo
          !$OMP END PARALLEL DO
       endif

       !$OMP PARALLEL WORKSHARE
       forall(i = sx:ex, j = sy:ey, k = 1:maxn-1)

          ! bulkE0: grad4(C)
          bulkE0(i, j, k) = (d2phi(i+1, j, k) - 2.d0*d2phi(i, j, k) + d2phi(i-1, j, k))/dx/dx &
               + (d2phi(i, j+1, k) - 2.d0*d2phi(i, j, k) + d2phi(i, j-1, k))/dy/dy &
               + (d2phi(i, j, k+1) - 2.d0*d2phi(i, j, k) + d2phi(i, j, k-1))/dz/dz

       endforall
       !$OMP END PARALLEL WORKSHARE

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, some, res, tid, temp5)
       tid = OMP_GET_THREAD_NUM()
       res = 0.d0

       !$OMP DO
       do i = sx, ex
          do j = sy, ey
             do k = 1, maxn-1
                some = 2.d0*phi(i, j, k) - 0.5d0*phiOld(i, j, k) &
                     - 1.5d0*phiNew(i, j, k) + dt*source0(i, j, k) &
                     + dt*nmda/2.d0/Pe*(tao*d2phi(i, j, k) - epn2*bulkE0(i, j, k))

                temp5 = phiNew(i, j, k)

                phiNew(i, j, k) = phiNew(i, j, k) + alfa1*some/diagonal

                if (phiNew(i, j, k).lt.0.d0) phiNew(i, j, k) = 0.d0
                if (phiNew(i, j, k).gt.1.d0) phiNew(i, j, k) = 1.d0

                res = dmax1(res, dabs(phiNew(i, j, k)-temp5))
             enddo
          enddo
       enddo
       !$OMP END DO

       res_vec(tid) = res
       !$OMP END PARALLEL

       !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j)
       do i = sx, ex
          do j = sy, ey
             ! z direction
             ! Non-periodic
             phiNew(i, j, 0) = phiNew(i, j, 1)
             phiNew(i, j, maxn) = phiNew(i, j, maxn-1)
          enddo
       enddo
       !$OMP END PARALLEL DO

       ! Exchange first order halos.
       call exchange2d(phiNew, stride_p_xz, stride_p_yz, neighbours,     &
            ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)

       ! TODO - WHAT SHOULD THE BOUNDARY CONDITIONS BE AT THE INLET?
       if (sx==1) then
          !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j, k)
          do j = sy, ey
             do k = 1, maxn-1
                ! x direction
                ! Non-periodic
                z_val = k*dz - 0.5d0*dz
                phiNew(0,j,k)=0.5d0+0.5d0*tanh((z_val-height)/(2.d0*(2d0**0.5d0)*epn))
             enddo
          enddo
          !$OMP END PARALLEL DO
       endif
       if (ex==(maxl-1)) then
          !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j, k)
          do j = sy, ey
             do k = 1, maxn-1
                ! x direction
                ! Non-periodic
                phiNew(maxl, j, k) = phiNew(maxl-1, j, k)
             enddo
          enddo
          !$OMP END PARALLEL DO
       endif

       res = 0.d0

       do tid = 0, num_threads-1

          if (res_vec(tid).gt.res) then
             res = res_vec(tid)
          endif

       enddo

    enddo ! End of SOR loop.

    ! write(6, *) 'Iterations required to converge (CH) = ', iteration

    deallocate(res_vec)
    if (alloc_stat.ne.0) STOP '** Unsuccessful deallocation **'

    return
  end subroutine dim

end module cahn_hilliard_solver
