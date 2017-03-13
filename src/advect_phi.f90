!> Advect-phi subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune,
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module advect_phi

  use tpls_mpi
  use weno

  implicit none

contains
  
  subroutine do_advect_phi(phi2, u2, v2, w2, conv2_phi,      &
       sx, ex, sy, ey, sz_p, ez_p, sz_uv, sz_w, ez_uv, ez_w, &
       maxl, maxn, ex_max, neighbours, comm2d_quasiperiodic, &
       dx, dy, dz,                                           &
       stride_p_augaug1_xz, stride_p_augaug1_yz,             &
       stride_p_augaug2_xz, stride_p_augaug2_yz,             &
       stride_p_augaug3_xz, stride_p_augaug3_yz,             &
       stride_uv_xz, stride_uv_yz,                           &
       stride_w_xz, stride_w_yz)

    implicit none

    ! Arguments.
    integer, intent(in) :: sx, ex, sy, ey, sz_p, ez_p, sz_uv, sz_w, ez_uv, ez_w
    integer, intent(in) :: maxl, maxn, ex_max
    integer, intent(in) :: neighbours(6)
    integer, intent(in) :: comm2d_quasiperiodic
    integer, intent(in) :: stride_p_augaug1_xz, stride_p_augaug1_yz
    integer, intent(in) :: stride_p_augaug2_xz, stride_p_augaug2_yz
    integer, intent(in) :: stride_p_augaug3_xz, stride_p_augaug3_yz
    integer, intent(in) :: stride_uv_xz, stride_uv_yz, stride_w_xz, stride_w_yz

    double precision, intent(in) :: dx, dy, dz
    double precision, intent(inout) :: phi2(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
    double precision, intent(inout) :: u2(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    double precision, intent(inout) :: v2(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    double precision, intent(inout) :: w2(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)
    double precision, intent(out) :: conv2_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)

    ! Local variables.
    integer :: i, j, k

    double precision :: phi_augaug(sx-3:ex+3,sy-3:ey+3,sz_p:ez_p) 
    double precision :: dfx(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
    double precision :: dfy(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv) 
    double precision :: dfz(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w) 

    ! Advect phi and reinitialise.

    ! Initialise augmented phi.
    phi_augaug(sx:ex, sy:ey, sz_p:ez_p) = phi2(sx:ex, sy:ey, sz_p:ez_p)

    ! Initialise augmented phi at x-boundaries.
    if(sx==1)then
       phi_augaug(0, sy:ey, sz_p:ez_p) = phi2(0, sy:ey, sz_p:ez_p)
    end if
    if(ex==ex_max)then
       phi_augaug(ex_max+1, sy:ey, sz_p:ez_p) = phi2(ex_max+1, sy:ey, sz_p:ez_p)
    end if

    ! Augmented phi: Exchange first-order halos.
    call exchange2d_augaug1(phi_augaug, stride_p_augaug1_xz, stride_p_augaug1_yz, &
         neighbours, ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)  

    ! Augmented phi: Exchange second-order halos.
    call exchange2d_augaug2(phi_augaug, stride_p_augaug2_xz, stride_p_augaug2_yz, &
         neighbours, ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)

    ! Augmented phi: Exchange third-order halos.
    call exchange2d_augaug3(phi_augaug, stride_p_augaug3_xz, stride_p_augaug3_yz, &
         neighbours, ex, ey, ez_p, sx, sy, sz_p, comm2d_quasiperiodic)

    ! Obtain upwinded derivatives of phi.
    ! For the interpolation of phi on to the u grid, quantities from
    !  i=0 to i=ex_max are needed. 
    ! For the other grids, quantities from i=1 to i=ex_max are
    !  needed.
    ! Extra if/else statements in the following subroutine take care
    !  of this.

    call WENOAppoxConvection(phi_augaug, u2, v2, w2, sx, ex, sy, ey, maxl, maxn, dfx, dfy, dfz)

    ! Non-augmented phi fluxes: exchange first-order halos only.  
    call exchange2d(dfx, stride_uv_xz, stride_uv_yz, neighbours,  &
         ex, ey, ez_uv, sx, sy, sz_uv, comm2d_quasiperiodic)
    call exchange2d(dfy, stride_uv_xz, stride_uv_yz, neighbours,  &
         ex, ey, ez_uv, sx, sy, sz_uv, comm2d_quasiperiodic)
    call exchange2d(dfz, stride_w_xz, stride_w_yz, neighbours,  &
         ex, ey, ez_w, sx, sy, sz_w, comm2d_quasiperiodic)

    ! Non-augmented velocities: exchange first-order halos only.     
    call exchange2d(u2, stride_uv_xz, stride_uv_yz, neighbours,  &
         ex, ey, ez_uv, sx, sy, sz_uv, comm2d_quasiperiodic)
    call exchange2d(v2,stride_uv_xz, stride_uv_yz, neighbours,  &
         ex, ey, ez_uv, sx, sy, sz_uv, comm2d_quasiperiodic)
    call exchange2d(w2, stride_w_xz, stride_w_yz, neighbours,  &
         ex, ey, ez_w, sx, sy, sz_w, comm2d_quasiperiodic)

    !$omp parallel workshare
    conv2_phi=0.d0

    ! Work directly on the phi-grid.

    ! conv2_phi correponds to RHS in the DIM code.
    forall (i = sx:ex, j = sy:ey, k = 1:maxn-1)
       conv2_phi(i, j, k) = (u2(i, j-1, k-1)*dfx(i, j-1, k-1) - u2(i-1, j-1, k-1)*dfx(i-1, j-1, k-1))/dx &
            + (v2(i, j, k-1)*dfy(i, j, k-1) - v2(i, j-1, k-1)*dfy(i, j-1, k-1))/dy         &
            + (w2(i, j-1, k)*dfz(i, j-1, k) - w2(i, j-1, k-1)*dfz(i, j-1, k-1))/dz
    end forall
    !$omp end parallel workshare
  end subroutine do_advect_phi

end module advect_phi
