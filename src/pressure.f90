!> Pressure subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module pressure

  implicit none

contains

  subroutine get_source_pres(rhs,u3,v3,w3,ex,ey,sx,sy,maxn,dx,dy,dz,dt)

    implicit none

    ! Arguments
    integer,intent(in) :: ex,ey,sx,sy,maxn

    double precision,intent(in) :: dx,dy,dz,dt
    double precision,intent(inout) :: rhs(sx-1:ex+1,sy-1:ey+1,0:maxn)
    double precision,intent(in) :: u3(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: v3(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: w3(sx-1:ex+1,sy-1:ey+1,0:maxn-1)

    ! Local variables
    integer :: i,j,k

    !$omp parallel default(shared), private(i,j,k)

    !$omp workshare
    rhs=0.d0
    !$omp end workshare

    !$omp do  
    do k=1,maxn-1
       do j=sy,ey
          do i=sx,ex
             rhs(i,j,k)=((u3(i,j-1,k-1)-u3(i-1,j-1,k-1))/dx)+ ((v3(i,j,k-1)-v3(i,j-1,k-1))/dy)+ &
                  ((w3(i,j-1,k)-w3(i,j-1,k-1))/dz)
          end do
       end do
    end do
    !$omp end do

    !$omp workshare
    rhs=rhs/dt
    !$omp end workshare

    !$omp end parallel

    return
  end subroutine get_source_pres


  subroutine  get_uvw_pres(u3,v3,w3,p,ex,ey,sx,sy,maxn,ex_max,dx,dy,dz,dt,u_inlet)

    implicit none

    ! Arguments
    integer,intent(in) :: ex,ey,sx,sy,maxn,ex_max

    double precision,intent(in) :: dx,dy,dz,dt
    double precision,intent(in) :: p(sx-1:ex+1,sy-1:ey+1,0:maxn)
    double precision,intent(inout) :: u3(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(inout) :: v3(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(inout) :: w3(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision,intent(in) :: u_inlet(0:maxn-2)

    ! Local variables
    integer :: i,j,k
    integer :: sx_loc,ex_loc

    if(sx==1)then
       sx_loc=sx-1
    else
       sx_loc=sx
    end if

    if(ex==ex_max)then
       ex_loc=ex+1
    else
       ex_loc=ex
    end if

    !$omp parallel default(shared), private(i,j,k)
    !$omp do  
    do k=0,maxn-2
       do j=sy,ey
          do i=sx_loc,ex_loc

             if(i==0)then
                u3(i,j,k)=u_inlet(k)
             elseif(i==ex_max+1)then
                ! last u3 point is completely fictitious and exists only for accounting reasons.
                u3(i,j,k)=0.d0
             else
                u3(i,j,k)=u3(i,j,k)-(dt/dx)*(p(i+1,j+1,k+1)-p(i,j+1,k+1))
             end if

          end do
       end do
    end do
    !$omp end do

    !$omp do  
    do k=0,maxn-2
       do j=sy,ey
          do i=sx_loc,ex_loc

             v3(i,j,k)=v3(i,j,k)-(dt/dy)*(p(i,j+1,k+1)-p(i,j,k+1))

          end do
       end do
    end do
    !$omp end do

    !$omp do  
    do k=0,maxn-1
       do j=sy,ey
          do i=sx_loc,ex_loc

             w3(i,j,k)=w3(i,j,k)-(dt/dz)*(p(i,j+1,k+1)-p(i,j+1,k))

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine get_uvw_pres

end module pressure
