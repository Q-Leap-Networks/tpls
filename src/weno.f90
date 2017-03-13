!> Weighted Essentially Non-Oscillatory, WENO, functions and subroutines.
!! Fifth order accurate WENO scheme designed by Jiang and C.-W. Shu 1996.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune,
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module weno

  implicit none

contains

  ! Compute the upwinding convective term.
  subroutine WENOAppoxConvection(phi, u, v, w, sx, ex, sy, ey, l, n, dfx,dfy, dfz)
    !call      WENOAppoxConvection(phi_augaug, u2, v2, w2, sx, ex, sy, ey, maxl, maxn, dfx, dfy, dfz)
    implicit none

    ! Arguments
    integer, intent(in) :: sx, ex, sy, ey, l, n

    double precision, intent(in) :: phi(sx-3:ex+3, sy-3:ey+3, 0:n)
    double precision, intent(in) :: u(sx-1:ex+1, sy-1:ey+1, 0:n-2)
    double precision, intent(in) :: v(sx-1:ex+1, sy-1:ey+1, 0:n-2)
    double precision, intent(in) :: w(sx-1:ex+1, sy-1:ey+1, 0:n-1)

    double precision, intent(out) :: dfx(sx-1:ex+1, sy-1:ey+1, 0:n-2)
    double precision, intent(out) :: dfy(sx-1:ex+1, sy-1:ey+1, 0:n-2)
    double precision, intent(out) :: dfz(sx-1:ex+1, sy-1:ey+1, 0:n-1)

    ! Local variables
    integer :: i, j, k, im2, im1, ip2, ip3, jm2, jm1, jp2, jp3

    double precision :: v1, v2, v3, v4, v5

    !$omp parallel default(shared), private(i, j, k, im2, im1, ip2, ip3, v1, v2, v3, v4, v5)
    !$omp do
    do k = 0, n-2
       do j = sy, ey
          do i = sx, ex

             if((i.ge.3).and.(i.le.(l-3))) then

                im2 = i-2
                im1 = i-1
                ip2 = i+2
                ip3 = i+3

                if(u(i, j, k).ge.0)then
                   v1 = phi(im2, j+1, k+1)
                   v2 = phi(im1, j+1, k+1)
                   v3 = phi(i, j+1, k+1)
                   v4 = phi(i+1, j+1, k+1)
                   v5 = phi(ip2, j+1, k+1)
                else
                   v1 = phi(ip3, j+1, k+1)
                   v2 = phi(ip2, j+1, k+1)
                   v3 = phi(i+1, j+1, k+1)
                   v4 = phi(i, j+1, k+1)
                   v5 = phi(im1, j+1, k+1)
                endif

                dfx(i, j, k) = WENO5(v1, v2, v3, v4, v5)

             else

                if(u(i, j, k).ge.0.d0)then
                   dfx(i, j, k) = phi(i, j+1, k+1)
                else
                   dfx(i, j, k) = phi(i+1, j+1, k+1)
                endif

             end if

          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! Take care of boundary point at i=0 as this is the only one needed

    if(sx==1)then
       !$omp parallel default(shared), private(j, k)
       !$omp do
       do k = 0, n-2
          do j = sy, ey
             if(u(0, j, k).ge.0.d0)then
                dfx(0, j, k) = phi(0, j+1, k+1)
             else
                dfx(0, j, k) = phi(1, j+1, k+1)
             endif
          end do
       end do
       !$omp end do
       !$omp end parallel
    end if

    !$omp parallel default(shared), private(i, j, k, jm2, jm1, jp2, jp3, v1, v2, v3, v4, v5)

    !$omp do
    do k = 0, n-2
       do j = sy, ey
          do i = sx, ex

             jm2 = j-2
             jm1 = j-1
             jp2 = j+2
             jp3 = j+3

             if(v(i, j, k).ge.0.d0)then
                v1 = phi(i, jm2, k+1)
                v2 = phi(i, jm1, k+1)
                v3 = phi(i, j, k+1)
                v4 = phi(i, j+1, k+1)
                v5 = phi(i, jp2, k+1)
             else
                v1 = phi(i, jp3, k+1)
                v2 = phi(i, jp2, k+1)
                v3 = phi(i, j+1, k+1)
                v4 = phi(i, j, k+1)
                v5 = phi(i, jm1, k+1)
             endif

             dfy(i, j, k) = WENO5(v1, v2, v3, v4, v5)

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do  
    do k = 0, n-1
       do j = sy, ey
          do i = sx, ex

             if((k.ge.3).and.(k.le.(n-3)))then

                if(w(i, j, k).ge.0)then
                   v1 = phi(i, j+1, k-2)
                   v2 = phi(i, j+1, k-1)
                   v3 = phi(i, j+1, k)
                   v4 = phi(i, j+1, k+1)
                   v5 = phi(i, j+1, k+2)
                else
                   v1 = phi(i, j+1, k+3)
                   v2 = phi(i, j+1, k+2)
                   v3 = phi(i, j+1, k+1)
                   v4 = phi(i, j+1, k)
                   v5 = phi(i, j+1, k-1)
                endif

                dfz(i, j, k) = WENO5(v1, v2, v3, v4, v5)

             else

                dfz(i, j, k) = (phi(i, j+1, k)+phi(i, j+1, k+1))/2.d0

             endif

          enddo
       enddo
    enddo
    !$omp end do

    !$omp end parallel

    return
  end subroutine WENOAppoxConvection


  function WENO5(v1, v2, v3, v4, v5)

    implicit none

    ! Arguments
    double precision :: WENO5
    double precision :: v1, v2, v3, v4, v5

    ! Local variables
    double precision, parameter :: epslon = 1.e-7
    double precision :: t1, t2, t3
    double precision :: s1, s2, s3
    double precision :: w1, w2, w3
    double precision :: a1, a2, a3

    t1 = v1 - 2.d0*v2 + v3
    t2 = v1 - 4.d0*v2 + 3.d0*v3
    s1 = 13.d0/12.d0*t1*t1 + 0.25d0*t2*t2

    t1 = v2 - 2.d0*v3 + v4
    t2 = v2 - v4
    s2 = 13.d0/12.d0*t1*t1 + 0.25d0*t2*t2

    t1 = v3 - 2.d0*v4 + v5
    t2 = 3.d0*v3 - 4.d0*v4 + v5
    s3 = 13.d0/12.d0*t1*t1 + 0.25d0*t2*t2

    a1 = 0.1d0/(epslon+s1)**2
    a2 = 0.6d0/(epslon+s2)**2
    a3 = 0.3d0/(epslon+s3)**2
    w1 = a1/(a1+a2+a3)
    w2 = a2/(a1+a2+a3)
    w3 = a3/(a1+a2+a3)

    t1 = v1/3.d0 - 7.d0*v2/6.d0 + 11.d0*v3/6.d0
    t2 = -v2/6.d0 + 5.d0*v3/6.d0 + v4/3.d0
    t3 = v3/3.d0 + 5.d0*v4/6.d0 - v5/6.d0

    WENO5 = w1*t1 + w2*t2 + w3*t3

  end function WENO5

end module weno
