!> Level-set subroutines.
!! The reinitialization step follows Russo and Smereka JCP 163, 51.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module levelset

  implicit none

contains

  subroutine get_levelset_status(phi, D_mx, istatus, sign_mx, sx, ex, sy, ey, &
      n, h)

    implicit none

    ! Arguments
    integer, intent(in) :: sx, ex, sy, ey, n

    double precision, intent(in) :: phi(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(out) :: sign_mx(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(out) :: D_mx(sx-1:ex+1, sy-1:ey+1, 0:n)

    integer, intent(out) :: istatus(sx-1:ex+1, sy-1:ey+1, 0:n)

    ! Local variables
    integer :: i, j, k

    double precision :: h, a1, a2, a3, a4, a5, a6, a7

    !$omp parallel default(shared), private(i, j, k, a1, a2, a3, a4, a5, a6, a7)

    !$omp workshare
    D_mx = 0.d0
    istatus = 0
    sign_mx = 0.d0
    !$omp end workshare

    !$omp do
    do k = 1, n-1
       do j = sy, ey
          do i = sx, ex

             ! Compute sign function of phi
             if (phi(i, j, k).lt.0.d0) then
                sign_mx(i, j, k) = -1d0
             else if (phi(i, j, k).gt.0.d0) then
                sign_mx(i, j, k) = 1d0
             else
                sign_mx(i, j, k) = 0.d0
             endif

             ! Compute indicator function \Sigma_{\Delta x} (=istatus)

             a1 = phi(i, j, k)*phi(i+1, j, k)
             a2 = phi(i, j, k)*phi(i-1, j, k)
             a3 = phi(i, j, k)*phi(i, j+1, k)
             a4 = phi(i, j, k)*phi(i, j-1, k)
             a5 = phi(i, j, k)*phi(i, j, k+1)
             a6 = phi(i, j, k)*phi(i, j, k-1)

             if (dmin1(a1, a2, a3, a4, a5, a6).le.0.d0) then
                istatus(i, j, k) = 1
             else
                istatus(i, j, k) = 0
             endif

             ! Compute distance function (D_mx)

             if (istatus(i, j, k).eq.1) then
                a1 = dsqrt((phi(i+1, j, k) - phi(i-1, j, k))**2 + (phi(i, j+1, k) - phi(i, j-1, k))**2 &
                     + (phi(i, j, k+1) - phi(i, j, k-1))**2)/2d0
                a2 = dabs(phi(i+1, j, k) - phi(i, j, k))
                a3 = dabs(phi(i-1, j, k) - phi(i, j, k))
                a4 = dabs(phi(i, j+1, k) - phi(i, j, k))
                a5 = dabs(phi(i, j-1, k) - phi(i, j, k))
                a6 = dabs(phi(i, j, k+1) - phi(i, j, k))
                a7 = dabs(phi(i, j, k-1) - phi(i, j, k))
                D_mx(i, j, k) = h*phi(i, j, k)/dmax1(a1, a2, a3, a4, a5, a6, a7, h)
             endif

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine get_levelset_status


  subroutine do_levelset_iteration(phi, phi_reinit, D_mx, istatus, sign_mx, &
      sx, ex, sy, ey, n, ex_max, h, height, error_val)

    implicit none

    ! Arguments
    integer, intent(in) :: sx, ex, sy, ey, n
    integer, intent(in) :: ex_max
    integer, intent(in) :: istatus(sx-1:ex+1, sy-1:ey+1, 0:n)

    double precision, intent(in) :: h ! Equivalent to "dx" in the main code.
    double precision, intent(in) :: height
    double precision, intent(in) :: phi(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(in) :: sign_mx(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(in) :: D_mx(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(inout) :: phi_reinit(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(out) :: error_val

    ! Local variables
    integer :: i, j, k

    double precision :: phi_reinit_updated(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision :: dtau, dc
    double precision :: aP, aM, bP, bM, cP, cM, dP, dM, eP, eM, fP, fM, z, dz, dy

    ! Fictitious time
    dtau = 0.3d0*h

    !$omp parallel default(shared), private(i, j, k, aP, bM, cP, dM, eP, fM, &
    !$omp& bP, aM, dP, cM, fP, eM)
    !$omp do
    do k = 1, n-1
       do j = sy, ey
          do i = sx, ex

             if (istatus(i, j, k).eq.1) then
                dc = sign_mx(i, j, k)*dabs(phi_reinit(i, j, k)) - D_mx(i, j, k)
             else

                dc = 0.d0

                if (phi(i, j, k).gt.0.d0) then
                   aP = dmax1(phi_reinit(i, j, k) - phi_reinit(i-1, j, k), 0.d0)
                   bM = dmin1(phi_reinit(i+1, j, k) - phi_reinit(i, j, k), 0.d0)
                   cP = dmax1(phi_reinit(i, j, k) - phi_reinit(i, j-1, k), 0.d0)
                   dM = dmin1(phi_reinit(i, j+1, k) - phi_reinit(i, j, k), 0.d0)
                   eP = dmax1(phi_reinit(i, j, k) - phi_reinit(i, j, k-1), 0.d0)
                   fM = dmin1(phi_reinit(i, j, k+1) - phi_reinit(i, j, k), 0.d0)
                   dc = sign_mx(i, j, k)*(dsqrt(dmax1(aP*aP, bM*bM) + dmax1(cP*cP, dM*dM) &
                        + dmax1(eP*eP, fM*fM)) - h)
                else if (phi(i, j, k).lt.0.d0) then
                   bP = dmax1(phi_reinit(i+1, j, k) - phi_reinit(i, j, k), 0.d0)
                   aM = dmin1(phi_reinit(i, j, k) - phi_reinit(i-1, j, k), 0.d0)
                   dP = dmax1(phi_reinit(i, j+1, k) - phi_reinit(i, j, k), 0.d0)
                   cM = dmin1(phi_reinit(i, j, k) - phi_reinit(i, j-1, k), 0.d0)
                   fP = dmax1(phi_reinit(i, j, k+1) - phi_reinit(i, j, k), 0.d0)
                   eM = dmin1(phi_reinit(i, j, k) - phi_reinit(i, j, k-1), 0.d0)
                   dc = sign_mx(i, j, k)*(dsqrt((dmax1(bP*bP, aM*aM)) + (dmax1(fP*fP, eM*eM)) &
                        + (dmax1(dP*dP, cM*cM))) - h)
                endif

             endif

             phi_reinit_updated(i, j, k) = phi_reinit(i, j, k) - dtau*dc/h

          end do
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    dz = h
    dy = h
    if (sx==1) then
       !$omp parallel default(shared), private(j, k)
       !$omp do
       do k = 1, n-1
          do j = sy, ey
             z = k*dz - 0.5d0*dz
             phi_reinit_updated(0, j, k) = z - height
          end do
       end do
       !$omp end do
       !$omp end parallel
    end if

    if (ex==ex_max) then
       !$omp parallel default(shared), private(k)
       !$omp do
       do k = 1, n-1
          phi_reinit_updated(ex_max+1, sy:ey, k) = phi_reinit_updated(ex_max, sy:ey, k)
       end do
       !$omp end do
       !$omp end parallel
    end if

    !$omp parallel workshare
    phi_reinit_updated(sx:ex, sy:ey, 0) = 2.d0*phi_reinit_updated(sx:ex, sy:ey, 1) &
         - phi_reinit_updated(sx:ex, sy:ey, 2)
    phi_reinit_updated(sx:ex, sy:ey, n) = 2.d0*phi_reinit_updated(sx:ex, sy:ey, n-1) &
         - phi_reinit_updated(sx:ex, sy:ey, n-2)
    !$omp end parallel workshare

    if (sx==1) then
       !$omp parallel workshare
       phi_reinit_updated(0:ex, sy:ey, 0) = 2.d0*phi_reinit_updated(0:ex, sy:ey, 1) &
            - phi_reinit_updated(0:ex, sy:ey, 2)
       phi_reinit_updated(0:ex, sy:ey, n) = 2.d0*phi_reinit_updated(0:ex, sy:ey, n-1) &
            - phi_reinit_updated(0:ex, sy:ey, n-2)
       !$omp end parallel workshare
    end if

    if (ex==ex_max) then
       !$omp parallel workshare
       phi_reinit_updated(sx:ex_max+1, sy:ey, 0) = 2.d0*phi_reinit_updated(sx:ex_max+1, sy:ey, 1) &
            - phi_reinit_updated(sx:ex_max+1, sy:ey, 2)
       phi_reinit_updated(sx:ex_max+1, sy:ey, n) = 2.d0*phi_reinit_updated(sx:ex_max+1, sy:ey, n-1) &
            - phi_reinit_updated(sx:ex_max+1, sy:ey, n-2)
       !$omp end parallel workshare
    end if

    error_val = 0d0
    !$omp parallel default(shared), private(i, j, k), reduction(MAX:error_val)

    !$omp do
    do k = 0, n
       do j = sy, ey
          do i = sx, ex
             error_val = error_val + (phi_reinit(i, j, k) - phi_reinit_updated(i, j, k))**2
          end do
       end do
    end do
    !$omp end do

    !$omp workshare
    phi_reinit = phi_reinit_updated
    !$omp end workshare nowait

    !$omp end parallel
    return
  end subroutine do_levelset_iteration


  !> Get CSF on phi-grid.
  subroutine get_csf(fx_csf, fy_csf, fz_csf, fi, sx, ex, sy, ey, &
      n, ex_max, dx, dy, dz, scap, smooth_width)

    implicit none

    ! Subroutine Arguments
    integer, intent(in) :: sx, ex, sy, ey, n, ex_max

    double precision, intent(in) :: dx, dy, dz, scap, smooth_width
    double precision, intent(in) :: fi(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(out) :: fx_csf(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(out) :: fy_csf(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision, intent(out) :: fz_csf(sx-1:ex+1, sy-1:ey+1, 0:n)

    ! Local Variables
    integer :: i, j, k, im1, ip1, jm1, jp1

    double precision :: curvature(sx-1:ex+1, sy-1:ey+1, 0:n)
    double precision :: term1, term2, term3
    double precision :: nrm_val, max_kappa
    double precision :: dfi_x, dfi_y, dfi_z
    double precision :: dirac_phi, phi_val
    double precision :: d2fi_x, d2fi_y, d2fi_z
    double precision :: dfi_x_up, dfi_x_down, d2fi_xz
    double precision :: dfi_y_up, dfi_y_down, d2fi_yz
    double precision :: dfi_x_in, dfi_x_out, d2fi_xy
    !      double precision, parameter :: pi = 3.14159265359d0
    double precision, parameter :: pi = 3.1415926535

    max_kappa = 1.d0/(dmin1(dx, dy, dz))

    !$omp parallel default(shared), private(i, j, k, ip1, im1, jp1, jm1,                &
    !$omp& d2fi_x, d2fi_y, d2fi_z, dfi_x, dfi_y, dfi_z, dfi_x_up, dfi_x_down, d2fi_xz,  &
    !$omp& dfi_y_up, dfi_y_down, d2fi_yz,                                               &
    !$omp& dfi_x_out, dfi_x_in, d2fi_xy, nrm_val, term1, term2, term3, phi_val, dirac_phi)

    !$omp do
    ! First, compute curvature
    do k = 1, n-1
       do j = sy, ey
          do i = sx, ex

             ip1 = i+1
             im1 = i-1
             jp1 = j+1
             jm1 = j-1

             d2fi_x = (fi(ip1, j, k) - 2.d0*fi(i, j, k) + fi(im1, j, k))/(dx**2.d0)
             d2fi_y = (fi(i, jp1, k) - 2.d0*fi(i, j, k) + fi(i, jm1, k))/(dy**2.d0)
             d2fi_z = (fi(i, j, k+1) - 2.d0*fi(i, j, k) + fi(i, j, k-1))/(dz**2.d0)

             dfi_x = (fi(ip1, j, k) - fi(im1, j, k))/(2.d0*dx)
             dfi_y = (fi(i, jp1, k) - fi(i, jm1, k))/(2.d0*dy)
             dfi_z = (fi(i, j, k+1) - fi(i, j, k-1))/(2.d0*dz)

             dfi_x_up = (fi(ip1, j, k+1) - fi(im1, j, k+1))/(2.d0*dx)
             dfi_x_down = (fi(ip1, j, k-1) - fi(im1, j, k-1))/(2.d0*dx)
             d2fi_xz = (dfi_x_up - dfi_x_down)/(2.d0*dz)

             dfi_y_up = (fi(i, jp1, k+1) - fi(i, jm1, k+1))/(2.d0*dy)
             dfi_y_down = (fi(i, jp1, k-1) - fi(i, jm1, k-1))/(2.d0*dy)
             d2fi_yz = (dfi_y_up - dfi_y_down)/(2.d0*dz)

             dfi_x_out = (fi(ip1, jp1, k) - fi(im1, jp1, k))/(2.d0*dx)
             dfi_x_in = (fi(ip1, jm1, k) - fi(im1, jm1, k))/(2.d0*dx)
             d2fi_xy = (dfi_x_out - dfi_x_in)/(2.d0*dy)

             nrm_val = dfi_x**2 + dfi_y**2 + dfi_z**2

             if (nrm_val.lt.1d-10) then
                curvature(i, j, k) = 0.d0
             else
                term1 = (d2fi_x + d2fi_y + d2fi_z)/(nrm_val**0.5d0)
                term2 = dfi_x*dfi_x*d2fi_x + dfi_y*dfi_y*d2fi_y + dfi_z*dfi_z*d2fi_z
                term3 = dfi_x*dfi_y*d2fi_xy + dfi_x*dfi_z*d2fi_xz + dfi_y*dfi_z*d2fi_yz
                term2 = term2/(nrm_val**1.5d0)
                term3 = 2.d0*term3/(nrm_val**1.5d0)
                curvature(i, j, k) = -(term1 - term2 - term3)
             end if

             if (curvature(i, j, k).gt.max_kappa) then
                curvature(i, j, k) = max_kappa
             end if

          end do
       end do
    end do
    !$omp end do

    ! curvature(:, :, n) = 2.d0*curvature(:, :, n-1) - curvature(:, :, n-2)
    ! curvature(:, :, 0) = 2.d0*curvature(:, :, 1) - curvature(:, :, 2)

    !$omp do
    do k = 1, n-1
       do j = sy, ey
          do i = sx, ex

             ip1 = i+1
             im1 = i-1
             jp1 = j+1
             jm1 = j-1

             dfi_x = (fi(ip1, j, k) - fi(im1, j, k))/(2.d0*dx)
             dfi_y = (fi(i, jp1, k) - fi(i, jm1, k))/(2.d0*dy)
             dfi_z = (fi(i, j, k+1) - fi(i, j, k-1))/(2.d0*dz)

             phi_val = fi(i, j, k)
             if (phi_val.lt.(-smooth_width)) then
                dirac_phi = 0.d0
             elseif (phi_val.gt.smooth_width) then
                dirac_phi = 0.d0
             else
                dirac_phi = (1.d0/(2.d0*smooth_width)) + (1.d0/(2.d0*smooth_width))*cos(pi*phi_val/smooth_width)
             end if

             fx_csf(i, j, k) = dirac_phi*scap*curvature(i, j, k)*dfi_x
             fy_csf(i, j, k) = dirac_phi*scap*curvature(i, j, k)*dfi_y
             fz_csf(i, j, k) = dirac_phi*scap*curvature(i, j, k)*dfi_z

          end do
       end do
    end do
    !$omp end do

    !$omp workshare
    fx_csf(sx:ex, sy:ey, 0) = 2.d0*fx_csf(sx:ex, sy:ey, 1) - fx_csf(sx:ex, sy:ey, 2)
    fx_csf(sx:ex, sy:ey, n) = 2.d0*fx_csf(sx:ex, sy:ey, n-1) - fx_csf(sx:ex, sy:ey, n-2)

    fy_csf(sx:ex, sy:ey, 0) = 2.d0*fy_csf(sx:ex, sy:ey, 1) - fy_csf(sx:ex, sy:ey, 2)
    fy_csf(sx:ex, sy:ey, n) = 2.d0*fy_csf(sx:ex, sy:ey, n-1) - fy_csf(sx:ex, sy:ey, n-2)

    fz_csf(sx:ex, sy:ey, 0) = 2.d0*fz_csf(sx:ex, sy:ey, 1) - fz_csf(sx:ex, sy:ey, 2)
    fz_csf(sx:ex, sy:ey, n) = 2.d0*fz_csf(sx:ex, sy:ey, n-1) - fz_csf(sx:ex, sy:ey, n-2)
    !$omp end workshare

    !$omp end parallel

    if (sx==1) then
       !$omp parallel workshare
       fx_csf(0, sy:ey, :) = 2.d0*fx_csf(1, sy:ey, :) - fx_csf(2, sy:ey, :)
       !$omp end parallel workshare
    end if

    if (ex==ex_max) then
       !$omp parallel workshare
       fx_csf(ex_max+1, sy:ey, :) = 2.d0*fx_csf(ex_max, sy:ey, :) - fx_csf(ex_max-1, sy:ey, :)
       !$omp end parallel workshare
    end if

    if (sx==1) then
       !$omp parallel workshare
       fy_csf(0, sy:ey, :) = 2.d0*fy_csf(1, sy:ey, :) - fy_csf(2, sy:ey, :)
       !$omp end parallel workshare
    end if

    if (ex==ex_max) then
       !$omp parallel workshare
       fy_csf(ex_max+1, sy:ey, :) = 2.d0*fy_csf(ex_max, sy:ey, :) - fy_csf(ex_max-1, sy:ey, :)
       !$omp end parallel workshare
    end if

    if (sx==1) then
       !$omp parallel workshare
       fz_csf(0, sy:ey, :) = 2.d0*fz_csf(1, sy:ey, :) - fz_csf(2, sy:ey, :)
       !$omp end parallel workshare
    end if

    if (ex==ex_max) then
       !$omp parallel workshare
       fz_csf(ex_max+1, sy:ey, :) = 2.d0*fz_csf(ex_max, sy:ey, :) - fz_csf(ex_max-1, sy:ey, :)
       !$omp end parallel workshare
    end if

  end subroutine get_csf

end module levelset
