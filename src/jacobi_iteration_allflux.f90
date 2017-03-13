!> JACOBI iteration subroutines.
!!
!! Subroutines for JACOBI iteration for the Crank--Nicholson and
!!  pressure steps. These subroutines use flux-conservative
!!   differencing, and the cell-averaged viscosities are copied from
!!    the master subroutine. 
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Toni Collis, Iain Bethune, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module jacobi_iteration_allflux

  implicit none

contains

  subroutine do_jacobi_u(u3,u3_old,RHS,viscosity,dx,dy,dz,dt, &
       ex,ey,sx,sy,maxn,ex_max)

    implicit none

    ! Arguments
    integer,intent(in) :: sx,sy,ex,ey
    integer,intent(in) :: maxn,ex_max

    double precision,intent(inout) :: u3(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: u3_old(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: RHS(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn)
    double precision,intent(in) :: dx,dy,dz,dt

    ! Local variables
    integer :: i,j,k,ex_loc

    double precision :: diag_val,ax,ay,az
    double precision :: mu_plushalf_x_val,mu_minushalf_x_val
    double precision :: mu_plushalf_y_val,mu_minushalf_y_val
    double precision :: mu_plushalf_z_val,mu_minushalf_z_val
    double precision :: u_minusz,u_plusz

    ax=dt/(dx*dx)
    ay=dt/(dy*dy)
    az=dt/(dz*dz)

    if(ex==ex_max)then
       ex_loc=ex-1
    else
       ex_loc=ex
    end if

    !$omp parallel default(shared), private(i,j,k,mu_plushalf_x_val,mu_minushalf_x_val, &
    !$omp& mu_plushalf_y_val,mu_minushalf_y_val,mu_plushalf_z_val,mu_minushalf_z_val,   &
    !$omp& diag_val,u_minusz,u_plusz)                                               
    !$omp do 
    do k=0,maxn-2
       do j=sy,ey
          do i=sx,ex_loc

             if(k.eq.0)then
                u_minusz=-2.d0*u3_old(i,j,0)+u3_old(i,j,1)/3.d0
             else
                u_minusz=u3_old(i,j,k-1)
             end if

             if(k.eq.(maxn-2))then
                u_plusz=-2.d0*u3_old(i,j,maxn-2)+u3_old(i,j,maxn-3)/3.d0
             else
                u_plusz=u3_old(i,j,k+1)
             end if

             mu_plushalf_x_val =viscosity(i+1,j+1,k+1)
             mu_minushalf_x_val=viscosity(i,  j+1,k+1)

             mu_plushalf_y_val =(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+2,k+1)+viscosity(i+1,j+2,k+1))/4.d0
             mu_minushalf_y_val=(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j  ,k+1)+viscosity(i+1,j  ,k+1))/4.d0

             mu_plushalf_z_val =(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k+2)+viscosity(i+1,j+1,k+2))/4.d0
             mu_minushalf_z_val=(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k  )+viscosity(i+1,j+1,k  ))/4.d0

             diag_val=1.d0+(ax/2.d0)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                  (ay/2.d0)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                  (az/2.d0)*(mu_plushalf_z_val+mu_minushalf_z_val)

             u3(i,j,k) = (1.d0/diag_val) * ( (ax/2.d0)*(mu_plushalf_x_val*u3_old(i+1,j,k)+mu_minushalf_x_val*u3_old(i-1,j,k)) + &
                  (ay/2.d0)*(mu_plushalf_y_val*u3_old(i,j+1,k)+mu_minushalf_y_val*u3_old(i,j-1,k)) + &
                  (az/2.d0)*(mu_plushalf_z_val*u_plusz    +mu_minushalf_z_val*u_minusz   ) + RHS(i,j,k))

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine do_jacobi_u


  subroutine do_jacobi_v(v3,v3_old,RHS,viscosity,dx,dy,dz,dt, &
       ex,ey,sx,sy,maxn,ex_max)

    implicit none

    ! Arguments
    integer,intent(in) :: sx,sy,ex,ey
    integer,intent(in) :: maxn,ex_max

    double precision,intent(inout) :: v3(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: v3_old(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: RHS(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn)
    double precision,intent(in) :: dx,dy,dz,dt

    ! Local variables
    integer :: i,j,k

    double precision :: diag_val,ax,ay,az
    double precision :: mu_plushalf_x_val,mu_minushalf_x_val
    double precision :: mu_plushalf_y_val,mu_minushalf_y_val
    double precision :: mu_plushalf_z_val,mu_minushalf_z_val
    double precision :: v_minusz,v_plusz,v_plusx,v_minusx

    ax=dt/(dx*dx)
    ay=dt/(dy*dy)
    az=dt/(dz*dz)

    !$omp parallel default(shared), private(i,j,k,mu_plushalf_x_val,mu_minushalf_x_val, &
    !$omp& mu_plushalf_y_val,mu_minushalf_y_val,mu_plushalf_z_val,mu_minushalf_z_val,   &
    !$omp& diag_val,v_minusz,v_plusz,v_minusx,v_plusx)                              
    !$omp do 
    do k=0,maxn-2
       do j=sy,ey
          do i=sx,ex

             if(k.eq.0)then
                v_minusz=-2.d0*v3_old(i,j,0)+v3_old(i,j,1)/3.d0
             else
                v_minusz=v3_old(i,j,k-1)
             end if

             if(k.eq.(maxn-2))then
                v_plusz=-2.d0*v3_old(i,j,maxn-2)+v3_old(i,j,maxn-3)/3.d0
             else
                v_plusz=v3_old(i,j,k+1)
             end if

             v_plusx=v3_old(i+1,j,k)
             v_minusx=v3_old(i-1,j,k)

             mu_plushalf_x_val= (viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i+1,j,k+1))/4.d0
             mu_minushalf_x_val=(viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i-1,j+1,k+1)+viscosity(i-1,j,k+1))/4.d0

             mu_plushalf_y_val= viscosity(i,j+1,k+1)
             mu_minushalf_y_val=viscosity(i,j,  k+1)

             mu_plushalf_z_val= (viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i,j+1,k+2)+viscosity(i,j,k+2))/4.d0
             mu_minushalf_z_val=(viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i,j+1,k  )+viscosity(i,j,k  ))/4.d0

             diag_val=1.d0+(ax/2.d0)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                  (ay/2.d0)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                  (az/2.d0)*(mu_plushalf_z_val+mu_minushalf_z_val)


             v3(i,j,k) = (1.d0/diag_val) * ( (ax/2.d0)*(mu_plushalf_x_val*v_plusx        +mu_minushalf_x_val*v_minusx       ) + &
                  (ay/2.d0)*(mu_plushalf_y_val*v3_old(i,j+1,k)+mu_minushalf_y_val*v3_old(i,j-1,k)) + &
                  (az/2.d0)*(mu_plushalf_z_val*v_plusz    +mu_minushalf_z_val*v_minusz           ) + &
                  RHS(i,j,k))

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    if(sx==1)then
       do k=0,maxn-2
          do j=sy,ey
             ! v3(0,j,k)=-(dx/dy)*(u3(0,j,k)-u3(0,j-1,k)) +(dt/dy)*(pres(1,j+1,k+1)-pres(1,j,k+1))
             v3(0,j,k)=v3(1,j,k)
          end do
       end do
    end if

    if(ex==ex_max)then
       v3(ex_max+1,sy:ey,:)=v3(ex_max,sy:ey,:)
    end if

    return
  end subroutine do_jacobi_v


  subroutine do_jacobi_w(w3,w3_old,RHS,viscosity,dx,dy,dz,dt, &
      ex,ey,sx,sy,maxn,ex_max)

    implicit none

    ! Arguments
    integer,intent(in) :: sx,sy,ex,ey
    integer,intent(in) :: maxn,ex_max

    double precision,intent(inout) :: w3(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision,intent(in) :: w3_old(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision,intent(in) :: RHS(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision,intent(in) :: viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn)
    double precision,intent(in) :: dx,dy,dz,dt

    ! Local variables
    integer :: i,j,k

    double precision :: diag_val,ax,ay,az
    double precision :: mu_plushalf_x_val,mu_minushalf_x_val
    double precision :: mu_plushalf_y_val,mu_minushalf_y_val
    double precision :: mu_plushalf_z_val,mu_minushalf_z_val
    double precision :: w_plusx,w_minusx

    ax=dt/(dx*dx)
    ay=dt/(dy*dy)
    az=dt/(dz*dz)

    !$omp parallel default(shared), private(i,j,k,mu_plushalf_x_val,mu_minushalf_x_val, &
    !$omp& mu_plushalf_y_val,mu_minushalf_y_val,mu_plushalf_z_val,mu_minushalf_z_val,   &
    !$omp& diag_val,w_minusx,w_plusx)                                               
    !$omp do  
    do k=1,maxn-2
       do j=sy,ey
          do i=sx,ex

             mu_plushalf_x_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k)+viscosity(i+1,j+1,k+1))/4.d0
             mu_minushalf_x_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i-1,j+1,k)+viscosity(i-1,j+1,k+1))/4.d0

             mu_plushalf_y_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j+2,k)+viscosity(i,j+2,k+1))/4.d0
             mu_minushalf_y_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j  ,k)+viscosity(i,j  ,k+1))/4.d0

             mu_plushalf_z_val=  viscosity(i,j+1,k+1)
             mu_minushalf_z_val= viscosity(i,j+1,k  )

             w_plusx=w3_old(i+1,j,k)
             w_minusx=w3_old(i-1,j,k)

             diag_val=1.d0+(ax/2.d0)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                  (ay/2.d0)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                  (az/2.d0)*(mu_plushalf_z_val+mu_minushalf_z_val)

             w3(i,j,k) = (1.d0/diag_val) * ( (ax/2.d0)*(mu_plushalf_x_val*w_plusx        +mu_minushalf_x_val*w_minusx       ) + &
                  (ay/2.d0)*(mu_plushalf_y_val*w3_old(i,j+1,k)+mu_minushalf_y_val*w3_old(i,j-1,k)) + &
                  (az/2.d0)*(mu_plushalf_z_val*w3_old(i,j,k+1)+mu_minushalf_z_val*w3_old(i,j,k-1)) + &
                  RHS(i,j,k))

          end do
       end do
    end do
    !$omp end do

    !$omp workshare
    w3(sx:ex,sy:ey,0)=0.d0
    w3(sx:ex,sy:ey,maxn-1)=0.d0
    !$omp end workshare nowait
    !$omp end parallel

    if(sx==1)then
       do k=1,maxn-2
          do j=sy,ey
             ! w3(0,j,k)=(dx/dz)*( (u_inlet(k)-u_inlet(k-1) )  -(u3(0,j,k)-u3(0,j,k-1))) &
             !  +(dt/dz)*(pres(1,j+1,k+1)-pres(1,j+1,k))
             w3(0,j,k)=w3(1,j,k)
          end do
       end do

       w3(0,sy:ey,0)=0.d0
       w3(0,sy:ey,maxn-1)=0.d0

    end if

    if(ex==ex_max)then
       w3(ex_max+1,sy:ey,:)=w3(ex_max,sy:ey,:)
    end if

    return
  end subroutine do_jacobi_w


  subroutine get_difference(ua,ub,ex,ey,ez,sx,sy,sz,diff)
    implicit none

    ! Arguments
    integer, intent(in) :: ex,ey,ez,sx,sy,sz

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez), intent(in) :: ua, ub
    double precision, intent(out) :: diff

    ! Local variables
    integer :: i,j,k

    diff = 0.0d0
    !$omp parallel default(shared), private(i,j,k), reduction(+:diff)
    !$omp do   
    do k=sz,ez
       do j = sy,ey
          do i = sx,ex
             diff = diff + (ua(i,j,k)-ub(i,j,k))**2
          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine get_difference

end module jacobi_iteration_allflux
