!> Momentum subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module momentum_allflux

  implicit none

contains

  !> Compute the viscosity term at level maxn.  
  !! This subroutine uses flux-conservative differencing.
  !! Also, it is the master subroutine where the cell-averaged
  !!  viscosities are computed. 
  !! To compute these cell-averages, the augmented viscosity array is
  !!  needed (MPI stuff). 
  subroutine get_diffusion_all(diff_u,diff_v,diff_w,viscosity,u2,v2,w2, &
       sx,sy,ex,ey,maxn,ex_max,dx,dy,dz)

    implicit none

    ! Arguments
    integer,intent(in) :: sx,sy,ex,ey,maxn,ex_max

    double precision,intent(in) :: u2(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(inout) :: diff_u(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: v2(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(inout) :: diff_v(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: w2(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision,intent(inout) :: diff_w(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision,intent(in) :: viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn)
    double precision,intent(in) :: dx,dy,dz

    ! Local variables
    integer :: i,j,k,ip1,im1,jp1,jm1,ex_loc

    double precision :: mu_duxE,mu_duxW,mu_duyN,mu_duyS,mu_duzI,mu_duzO
    double precision :: mu_dvxE,mu_dvxW,mu_dvyN,mu_dvyS,mu_dvzI,mu_dvzO
    double precision :: mu_dwxE,mu_dwxW,mu_dwyN,mu_dwyS,mu_dwzI,mu_dwzO
    double precision :: u_minusz,u_plusz
    double precision :: v_minusz,v_plusz
    double precision :: v_plusx,v_minusx
    double precision :: w_plusx,w_minusx
    double precision :: mu_plushalf_x_val,mu_minushalf_x_val
    double precision :: mu_plushalf_y_val,mu_minushalf_y_val
    double precision :: mu_plushalf_z_val,mu_minushalf_z_val

    !$omp parallel workshare
    diff_u=0.d0
    !$omp end parallel workshare

    if(ex==ex_max)then
       ex_loc=ex-1
    else
       ex_loc=ex
    end if

    !$omp parallel default(shared), private(i,j,k,ip1,im1,jp1,jm1,                    &
    !$omp& u_minusz,u_plusz,v_minusz,v_plusz,mu_plushalf_x_val,mu_minushalf_x_val,    &
    !$omp& mu_plushalf_y_val,mu_minushalf_y_val,mu_plushalf_z_val,mu_minushalf_z_val, &
    !$omp& mu_duxE,mu_duxW,mu_duyN,mu_duyS,mu_duzO,mu_duzI,                           &
    !$omp& mu_dvxE,mu_dvxW,mu_dvyN,mu_dvyS,mu_dvzO,mu_dvzI,                           &
    !$omp& mu_dwxE,mu_dwxW,mu_dwyN,mu_dwyS,mu_dwzO,mu_dwzI)
    !$omp do  
    do k=0,maxn-2
       do j=sy,ey
          do i=sx,ex_loc

             ip1=i+1
             im1=i-1
             jp1=j+1
             jm1=j-1

             if(k.eq.0)then
                u_minusz=-2.d0*u2(i,j,0)+u2(i,j,1)/3.d0
             else
                u_minusz=u2(i,j,k-1)
             end if

             if(k.eq.(maxn-2))then
                u_plusz=-2.d0*u2(i,j,maxn-2)+u2(i,j,maxn-3)/3.d0
             else
                u_plusz=u2(i,j,k+1)
             end if

             mu_plushalf_x_val =viscosity(i+1,j+1,k+1)
             mu_minushalf_x_val=viscosity(i,  j+1,k+1)

             mu_plushalf_y_val =(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+2,k+1)+viscosity(i+1,j+2,k+1))/4.d0
             mu_minushalf_y_val=(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j  ,k+1)+viscosity(i+1,j  ,k+1))/4.d0

             mu_plushalf_z_val =(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k+2)+viscosity(i+1,j+1,k+2))/4.d0
             mu_minushalf_z_val=(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k  )+viscosity(i+1,j+1,k  ))/4.d0

             mu_duxE= mu_plushalf_x_val*((u2(ip1,j,k)-u2(i,j,k))/dx)
             mu_duxW=mu_minushalf_x_val*((u2(i,j,k)-u2(im1,j,k))/dx)

             mu_duyN= mu_plushalf_y_val*((u2(i,jp1,k)-u2(i,j,k))/dy)
             mu_duyS=mu_minushalf_y_val*((u2(i,j,k)-u2(i,jm1,k))/dy)

             mu_duzO= mu_plushalf_z_val*((u_plusz-u2(i,j,k))/dz)
             mu_duzI=mu_minushalf_z_val*((u2(i,j,k)-u_minusz)/dz)

             diff_u(i,j,k)=((mu_duxE-mu_duxW)/dx)+((mu_duyN-mu_duyS)/dy)+((mu_duzO-mu_duzI)/dz)

          end do
       end do
    end do
    !$omp end do

    !$omp workshare
    diff_v=0.d0
    !$omp end workshare

    !$omp do  
    do k=0,maxn-2
       do j=sy,ey
          do i=sx,ex

             ip1=i+1
             im1=i-1
             jp1=j+1
             jm1=j-1

             if(k.eq.0)then
                v_minusz=-2.d0*v2(i,j,0)+v2(i,j,1)/3.d0
             else
                v_minusz=v2(i,j,k-1)
             end if

             if(k.eq.(maxn-2))then
                v_plusz=-2.d0*v2(i,j,maxn-2)+v2(i,j,maxn-3)/3.d0
             else
                v_plusz=v2(i,j,k+1)
             end if

             v_minusx=v2(i-1,j,k)
             v_plusx=v2(i+1,j,k)

             mu_plushalf_x_val= (viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i+1,j,k+1))/4.d0
             mu_minushalf_x_val=(viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i-1,j+1,k+1)+viscosity(i-1,j,k+1))/4.d0

             mu_plushalf_y_val= viscosity(i,j+1,k+1)
             mu_minushalf_y_val=viscosity(i,j,  k+1)

             mu_plushalf_z_val= (viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i,j+1,k+2)+viscosity(i,j,k+2))/4.d0
             mu_minushalf_z_val=(viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i,j+1,k  )+viscosity(i,j,k  ))/4.d0

             mu_dvxE= mu_plushalf_x_val*((v_plusx    -v2(i,  j,k))/dx)
             mu_dvxW=mu_minushalf_x_val*((v2(i,  j,k)-v_minusx   )/dx)

             mu_dvyN= mu_plushalf_y_val*((v2(i,jp1,k)-v2(i,j,k))/dy)
             mu_dvyS=mu_minushalf_y_val*((v2(i,j,k)-v2(i,jm1,k))/dy)

             mu_dvzO= mu_plushalf_z_val*((v_plusz-v2(i,j,k))/dz)
             mu_dvzI=mu_minushalf_z_val*((v2(i,j,k)-v_minusz)/dz)

             diff_v(i,j,k)=((mu_dvxE-mu_dvxW)/dx)+((mu_dvyN-mu_dvyS)/dy)+((mu_dvzO-mu_dvzI)/dz)

          end do
       end do
    end do
    !$omp end do

    !$omp workshare
    diff_w=0.d0
    !$omp end workshare

    !$omp do  
    do k=1,maxn-2
       do j=sy,ey
          do i=sx,ex

             ip1=i+1
             im1=i-1
             jp1=j+1
             jm1=j-1

             w_minusx=w2(i-1,j,k)
             w_plusx=w2(i+1,j,k)

             mu_plushalf_x_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k)+viscosity(i+1,j+1,k+1))/4.d0
             mu_minushalf_x_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i-1,j+1,k)+viscosity(i-1,j+1,k+1))/4.d0

             mu_plushalf_y_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j+2,k)+viscosity(i,j+2,k+1))/4.d0
             mu_minushalf_y_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j  ,k)+viscosity(i,j  ,k+1))/4.d0

             mu_plushalf_z_val=  viscosity(i,j+1,k+1)
             mu_minushalf_z_val= viscosity(i,j+1,k  )

             mu_dwxE= mu_plushalf_x_val*((w_plusx    -w2(i,  j,k))/dx)
             mu_dwxW=mu_minushalf_x_val*((w2(i,  j,k)-w_minusx   )/dx)

             mu_dwyN= mu_plushalf_y_val*((w2(i,jp1,k)-w2(i,j,  k))/dy)
             mu_dwyS=mu_minushalf_y_val*((w2(i,j,  k)-w2(i,jm1,k))/dy)

             mu_dwzO= mu_plushalf_z_val*((w2(i,j,k+1)-w2(i,j,k  ))/dz)
             mu_dwzI=mu_minushalf_z_val*((w2(i,j,k  )-w2(i,j,k-1))/dz)

             diff_w(i,j,k)=((mu_dwxE-mu_dwxW)/dx)+((mu_dwyN-mu_dwyS)/dy)+((mu_dwzO-mu_dwzI)/dz)

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine get_diffusion_all


  !> Compute the convective term at level maxn.  
  !! There are two parts to the convective term: the ordinary
  !!  convective term, and a part involving gradients of the
  !!   viscosity (the transposed part of the stress tensor). 
  !! The differencing is fully flux-conservative.
  !! To compute the cell-averaged viscosities, the augmented
  !!  viscosity array is needed (MPI stuff). 
  subroutine get_conv_all(conv_u2,conv_v2,conv_w2,viscosity_aug, &
      u2,v2,w2,sx,sy,ex,ey,maxn,ex_max,dx,dy,dz)

    implicit none

    ! Arguments
    integer,intent(in):: sx,ex,sy,ey,ex_max,maxn

    double precision,intent(in) :: dx,dy,dz
    double precision,intent(inout) :: conv_u2(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: u2(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(inout) :: conv_v2(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(in) :: v2(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision,intent(inout) :: conv_w2(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision,intent(in) :: w2(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision,intent(in) :: viscosity_aug(sx-2:ex+2,sy-2:ey+2,0:maxn)

    ! Local variables

    integer :: i,j,k,im1,ip1,jm1
    integer :: ex_loc

    double precision :: ub,vb,tempu,tempv,tempw,temp_conv1,temp_conv2
    double precision :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    double precision :: Xconv,Yconv,Zconv
    double precision :: mu_dudxE,mu_dudxW,mu_dudyE,mu_dudyW,mu_dudzE,mu_dudzW
    double precision :: mu_dvdxN,mu_dvdxS,mu_dvdyN,mu_dvdyS,mu_dvdzN,mu_dvdzS
    double precision :: mu_dwdxO,mu_dwdxI,mu_dwdyO,mu_dwdyI,mu_dwdzO,mu_dwdzI
    double precision :: mu_plushalf_x,mu_minushalf_x
    double precision :: mu_plushalf_y,mu_minushalf_y
    double precision :: mu_plushalf_z,mu_minushalf_z
    double precision :: v_plusx,v_minusx
    double precision :: w_plusx,w_minusx

    ! ******************************  X direction

    !$omp parallel workshare
    conv_u2=0.d0
    !$omp end parallel workshare

    if(ex==ex_max)then
       ex_loc=ex-1
    else
       ex_loc=ex
    end if

    !$omp parallel default(shared), private(i,j,k,im1,ip1,jm1,        &
    !$omp& tempu,tempv,tempw,Xconv,Yconv,Zconv,temp_conv1,temp_conv2, &
    !$omp& mu_plushalf_x,mu_minushalf_x,mu_plushalf_y,mu_minushalf_y, &
    !$omp& mu_plushalf_z,mu_minushalf_z,                              &
    !$omp& dudx,dudy,dudz,ub,vb,                                      &
    !$omp& dvdx,dvdy,dvdz,v_minusx,v_plusx,                           &
    !$omp& dwdx,dwdy,dwdz,w_minusx,w_plusx,                           &
    !$omp& mu_dudxE,mu_dudxW,mu_dvdxN,mu_dvdxS,mu_dwdxO,mu_dwdxI,     &
    !$omp& mu_dudyE,mu_dudyW,mu_dvdyN,mu_dvdyS,mu_dwdyO,mu_dwdyI,     &
    !$omp& mu_dudzE,mu_dudzW,mu_dvdzN,mu_dvdzS,mu_dwdzO,mu_dwdzI)
    !$omp do  
    do k=0, maxn-2
       do j=sy, ey
          do i=sx, ex_loc

             im1=i-1
             ip1=i+1

             tempv=(v2(i,j,k)+v2(ip1,j,k)+v2(i,j+1,k)+v2(ip1,j+1,k))/4.d0
             tempw=(w2(i,j,k)+w2(ip1,j,k)+w2(i,j,k+1)+w2(ip1,j,k+1))/4.d0

             dudx=0.5d0*(u2(i+1,j,k)-u2(im1,j,k))/dx           
             Xconv=dudx*u2(i,j,k)

             ub=u2(i,j+1,k)
             vb=u2(i,j-1,k)
             dudy=0.5d0*(ub-vb)/dy

             Yconv=dudy*tempv

             if(k==0)then
                dudz= (u2(i,j,k)+u2(i,j,k+1))/(2.d0*dz)
             elseif(k==maxn-2)then
                dudz=-(u2(i,j,k)+u2(i,j,k-1))/(2.d0*dz)
             else
                ub=u2(i,j,k+1)
                vb=u2(i,j,k-1)
                dudz=0.5d0*(ub-vb)/dz
             endif

             Zconv=dudz*tempw
             temp_conv1=Xconv+Yconv+Zconv

             mu_plushalf_x=viscosity_aug(i+1,j+1,k+1)
             mu_minushalf_x=viscosity_aug(i,j+1,k+1)

             mu_plushalf_y = (viscosity_aug(i+1,j+2,k+1)+ viscosity_aug(i,j+2,k+1)+ viscosity_aug(i+1,j+1,k+1)+ &
                  viscosity_aug(i,j+1,k+1))/4.d0
             mu_minushalf_y= (viscosity_aug(i+1,j,k+1)  + viscosity_aug(i,j,k+1)  + viscosity_aug(i+1,j+1,k+1)+ &
                  viscosity_aug(i,j+1,k+1))/4.d0

             mu_plushalf_z = (viscosity_aug(i+1,j+1,k+2)+ viscosity_aug(i,j+1,k+2)+ viscosity_aug(i+1,j+1,k+1)+ &
                  viscosity_aug(i,j+1,k+1))/4.d0
             mu_minushalf_z= (viscosity_aug(i+1,j+1,k)  + viscosity_aug(i,j+1,k)  + viscosity_aug(i+1,j+1,k+1)+ &
                  viscosity_aug(i,j+1,k+1))/4.d0

             mu_dudxE=mu_plushalf_x* (u2(i+1,j,k)-u2(i,  j,k))/dx
             mu_dudxW=mu_minushalf_x*(u2(i,  j,k)-u2(i-1,j,k))/dx

             mu_dvdxN=mu_plushalf_y* (v2(i+1,j+1,k)-v2(i,j+1,k))/dx
             mu_dvdxS=mu_minushalf_y*(v2(i+1,j,  k)-v2(i,j,  k))/dx

             mu_dwdxO=mu_plushalf_z* (w2(i+1,j,k+1)-w2(i,j,k+1))/dx
             mu_dwdxI=mu_minushalf_z*(w2(i+1,j,k  )-w2(i,j,k  ))/dx

             temp_conv2=((mu_dudxE-mu_dudxW)/dx)+((mu_dvdxN-mu_dvdxS)/dy)+((mu_dwdxO-mu_dwdxI)/dz)

             conv_u2(i,j,k)=temp_conv1-temp_conv2
          enddo
       enddo
    enddo
    !$omp end do

    ! ******************************  Y direction 

    !$omp workshare
    conv_v2=0.d0
    !$omp end workshare

    !$omp do  
    do k=0, maxn-2
       do j=sy, ey
          do i=sx, ex

             jm1=j-1
             tempu=(u2(i-1,jm1,k)+u2(i-1,j,k)+u2(i,jm1,k)+u2(i,j,k))/4.d0
             tempw=(w2(i,jm1,k)+w2(i,j,k)+w2(i,jm1,k+1)+w2(i,j,k+1))/4.d0

             dvdy=0.5d0*(v2(i,j+1,k)-v2(i,jm1,k))/dy
             Yconv=dvdy*v2(i,j,k)

             v_minusx=v2(i-1,j,k)
             v_plusx=v2(i+1,j,k)

             dvdx=0.5d0*(v_plusx-v_minusx)/dx
             Xconv=dvdx*tempu

             if(k==0)then
                dvdz=(v2(i,j,k)+v2(i,j,k+1))/(2.d0*dz)
             elseif(k==maxn-2)then
                dvdz=-(v2(i,j,k)+v2(i,j,k-1))/(2.d0*dz)
             else
                dvdz=0.5d0*(v2(i,j,k+1)-v2(i,j,k-1))/dz
             endif

             Zconv=dvdz*tempw
             temp_conv1=Xconv+Yconv+Zconv

             mu_plushalf_x = (  viscosity_aug(i+1,j+1,k+1)+   viscosity_aug(i+1,j,k+1) +  viscosity_aug(i,j+1,k+1)+  &
                  viscosity_aug(i,j,k+1))/4.d0
             mu_minushalf_x = (  viscosity_aug(i-1,j+1,k+1)  +   viscosity_aug(i-1,j,k+1)   +  viscosity_aug(i,j+1,k+1)+  &
                  viscosity_aug(i,j,k+1))/4.d0

             mu_plushalf_y= viscosity_aug(i,j+1,k+1)
             mu_minushalf_y=viscosity_aug(i,j,  k+1)

             mu_plushalf_z = (viscosity_aug(i,j+1,k+2) +viscosity_aug(i,j,k+2)  +  viscosity_aug(i,j+1,k+1)+  &
                  viscosity_aug(i,j,k+1))/4.d0
             mu_minushalf_z = (viscosity_aug(i,j+1,k)   +viscosity_aug(i,j,k)    +  viscosity_aug(i,j+1,k+1)+  &
                  viscosity_aug(i,j,k+1))/4.d0

             mu_dudyE=mu_plushalf_x* (u2(i,j,k)-u2(i,j-1,k))/dy
             mu_dudyW=mu_minushalf_x*(u2(i-1,  j,k)-u2(i-1,  j-1,k))/dy

             mu_dvdyN=mu_plushalf_y* (v2(i,j+1,k)-v2(i,j,  k))/dy
             mu_dvdyS=mu_minushalf_y*(v2(i,j,  k)-v2(i,j-1,k))/dy

             mu_dwdyO=mu_plushalf_z* (w2(i,j,k+1)-w2(i,j-1,k+1))/dy
             mu_dwdyI=mu_minushalf_z*(w2(i,j,k  )-w2(i,j-1,k  ))/dy

             temp_conv2=((mu_dudyE-mu_dudyW)/dx)+((mu_dvdyN-mu_dvdyS)/dy)+((mu_dwdyO-mu_dwdyI)/dz)

             conv_v2(i,j,k)=temp_conv1-temp_conv2
          enddo
       enddo
    enddo
    !$omp end do

    ! ******************************  Z direction  

    !$omp workshare
    conv_w2=0.d0
    !$omp end workshare

    !$omp do  
    do k=1, maxn-2
       do j=sy, ey
          do i=sx, ex
             tempu=(u2(i-1,j,k-1)+u2(i-1,j,k)+u2(i,j,k-1)+u2(i,j,k))/4.d0
             tempv=(v2(i,j,k-1)+v2(i,j,k)+v2(i,j+1,k-1)+v2(i,j+1,k))/4.d0

             ub = w2(i,j,k+1)
             vb = w2(i,j,k-1)
             dwdz=0.5*(ub-vb)/dz
             Zconv=dwdz*w2(i,j,k)

             w_minusx=w2(i-1,j,k)
             w_plusx=w2(i+1,j,k)

             dwdx=0.5d0*(w_plusx-w_minusx)/dx
             Xconv=dwdx*tempu

             ub=w2(i,j+1,k)
             vb=w2(i,j-1,k)
             dwdy=0.5d0*(ub-vb)/dy
             Yconv=dwdy*tempv

             temp_conv1=Xconv+Yconv+Zconv

             mu_plushalf_x=(viscosity_aug(i+1,j+1,k+1)+  viscosity_aug(i+1,j+1,k) +  viscosity_aug(i,j+1,k+1)+  &
                  viscosity_aug(i,j+1,k))/4.d0
             mu_minushalf_x=(viscosity_aug(i-1,j+1,k+1)  +  viscosity_aug(i-1,j+1,k)   +  viscosity_aug(i,j+1,k+1)+  &
                  viscosity_aug(i,j+1,k))/4.d0

             mu_plushalf_y= (viscosity_aug(i,j+2,k+1)+  viscosity_aug(i,j+2,k)+  viscosity_aug(i,j+1,k+1)+  &
                  viscosity_aug(i,j+1,k))/4.d0
             mu_minushalf_y= (viscosity_aug(i,j,k+1)  +  viscosity_aug(i,j,k)  +  viscosity_aug(i,j+1,k+1)+  &
                  viscosity_aug(i,j+1,k))/4.d0

             mu_plushalf_z=viscosity_aug(i,j+1,k+1)
             mu_minushalf_z=viscosity_aug(i,j+1,k)

             mu_dudzE=mu_plushalf_x* (u2(i,j,k)-u2(i,j,k-1))/dz
             mu_dudzW=mu_minushalf_x*(u2(i-1,  j,k)-u2(i-1,  j,k-1))/dz

             mu_dvdzN=mu_plushalf_y* (v2(i,j+1,k)-v2(i,j+1,k-1))/dz
             mu_dvdzS=mu_minushalf_y*(v2(i,j,  k)-v2(i,j  ,k-1))/dz

             mu_dwdzO=mu_plushalf_z* (w2(i,j,k+1)-w2(i,j,k  ))/dz
             mu_dwdzI=mu_minushalf_z*(w2(i,j,k  )-w2(i,j,k-1))/dz

             temp_conv2=((mu_dudzE-mu_dudzW)/dx)+((mu_dvdzN-mu_dvdzS)/dy)+((mu_dwdzO-mu_dwdzI)/dz)

             conv_w2(i,j,k)=temp_conv1-temp_conv2
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine get_conv_all


  subroutine get_viscosity(viscosity,phi,ex,ey,sx,sy,maxn,smooth_width,Re,mu_minus,mu_plus)

    implicit none

    ! Arguments
    integer,intent(in) :: sx,ex,sy,ey,maxn

    double precision,intent(inout) :: viscosity(sx-1:ex+1,sy-1:ey+1,0:maxn)
    double precision,intent(in) :: phi(sx-1:ex+1,sy-1:ey+1,0:maxn)
    double precision,intent(in):: smooth_width,Re,mu_minus,mu_plus

    ! Local variables
    integer :: i,j,k

    double precision :: heaviside_phi
    double precision :: pi=3.1415926535
    double precision :: phi_val,temp

    !$omp parallel default(shared), private(i,j,k, &
    !$omp& phi_val,heaviside_phi,temp)
    !$omp do  
    do k=0,maxn
       do j=sy,ey
          do i=sx,ex

             phi_val=phi(i,j,k)
             if(phi_val.lt.-smooth_width)then
                heaviside_phi=0.d0
             elseif(phi_val.gt.smooth_width)then
                heaviside_phi=1.d0
             else
                heaviside_phi=0.5d0+(phi_val/(2.d0*smooth_width))+(1.d0/(2.d0*pi))*sin(pi*phi_val/smooth_width)
             end if

             !temp=(1.d0/mu_minus)+((1.d0/mu_plus)-(1.d0/mu_minus))*heaviside_phi
             !viscosity(i,j,k)=(1.d0/temp)/Re

             temp=mu_plus*heaviside_phi+mu_minus*(1.d0-heaviside_phi)
             viscosity(i,j,k)=temp/Re

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine get_viscosity

end module momentum_allflux
