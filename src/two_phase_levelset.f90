!> Two-phase level-set subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module two_phase_levelset

  use levelset
  use tpls_mpi

  implicit none

contains

  ! Note that sz_p = 0 and ez_p = maxn.
  subroutine do_tpls(phi2, phi3, conv0_phi, conv1_phi, conv2_phi, &
       sx, ex, sy, ey,                                            &
       sz_p, ez_p, maxn, ex_max, stride_p_xz, stride_p_yz,        &
       max_iteration_levelset, neighbours, comm2d_quasiperiodic,  &
       height, err2,                                              &
       dx, dz, dt)

    implicit none

    ! Arguments.
    integer, intent(in) :: sx, ex, sy, ey, sz_p, ez_p
    integer, intent(in) :: maxn, ex_max
    integer, intent(in) :: stride_p_xz, stride_p_yz
    integer, intent(in) :: max_iteration_levelset
    integer, intent(in) :: neighbours(6)
    integer, intent(in) :: comm2d_quasiperiodic
    double precision, intent(in) :: dx, dz, dt
    double precision, intent(in) :: height
    double precision, intent(inout) :: err2
    double precision, intent(inout) :: phi2(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    double precision, intent(inout) :: conv0_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    double precision, intent(inout) :: conv1_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    double precision, intent(inout) :: conv2_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    double precision, intent(inout) :: phi3(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)

    ! Local variables.
    integer :: i, j, k
    integer :: startx, endx
    integer :: iteration_sor
    integer :: istatus_mx(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)
    double precision :: z_val
    double precision :: D_mx(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)
    double precision :: sign_mx(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)
    double precision :: phi_reinit(sx-1:ex+1, sy-1:ey+1, sz_p:ez_p)
    double precision :: RHS_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)

    !$omp parallel workshare
    !! AB3:
    RHS_phi = (23.d0/12.d0)*conv2_phi - (4.d0/3.d0)*conv1_phi + (5.d0/12.d0)*conv0_phi
    !$omp end parallel workshare

    ! Non-augmented RHS of phi: exchange first-order halos only (default now is with corners).
    call exchange2d(RHS_phi, stride_p_xz, stride_p_yz, neighbours, &
         ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)

    !$omp parallel workshare
    ! phi3 and ph2 correspond to finew and fi in the DIM code.
    forall (i = sx:ex, j = sy:ey, k = 1:maxn-1)
       phi3(i, j, k) = phi2(i, j, k) - dt*RHS_phi(i, j, k)
    end forall
    !$omp end parallel workshare

    if(sx==1)then
       do k = 1, maxn-1
          z_val = k*dz - 0.5d0*dz
          phi3(0, sy:ey, k) = (z_val-height)
       end do
    end if
    if(ex==ex_max)then
       do k = 1, maxn-1
          phi3(ex_max+1, sy:ey, k) = phi3(ex_max, sy:ey, k)
       end do
    end if

    if (sx==1) then
       startx = 0
    else
       startx = sx
    endif
    if (ex==ex_max) then
       endx = ex_max+1
    else
       endx = ex
    endif
    phi3(startx:endx, sy:ey, 0) = 2.d0*phi3(startx:endx, sy:ey, 1) - phi3(startx:endx, sy:ey, 2)
    phi3(startx:endx, sy:ey, maxn) = 2.d0*phi3(startx:endx, sy:ey, maxn-1) - phi3(startx:endx, sy:ey, maxn-2)

    ! Non-augmented phi: exchange first-order halos only.
    ! Here, phi3 means the same thing as phi0 in the paper of Smereka and Russo.

    call exchange2d(phi3, stride_p_xz, stride_p_yz, neighbours,  &
         ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)

    call get_levelset_status(phi3, D_mx, istatus_mx, sign_mx, sx, ex, sy, ey, maxn, dx)

    phi_reinit=phi3
    do iteration_sor = 1, max_iteration_levelset
       call exchange2d(phi_reinit, stride_p_xz, stride_p_yz, neighbours, &
            ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)
       call do_levelset_iteration(phi3, phi_reinit, D_mx, istatus_mx, sign_mx, &
            sx, ex, sy, ey, maxn, ex_max, dx, height, err2)
    end do

    ! Now update phi3 to its redistanced value.
    !$omp parallel workshare
    phi3 = phi_reinit
    !$omp end parallel workshare

  end subroutine do_tpls

end module two_phase_levelset
