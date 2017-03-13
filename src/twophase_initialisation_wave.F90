!> Two-phase level-set initial conditions subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune,
!! Toni Collis, David Scott, Peter Spelt, Mike Jackson.
!! @version $Revision: 328 $
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt, The University of
!! Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module twophase_initialisation_wave

  use option_names

  implicit none

contains

  !> Get the initial conditions of the velocity, pressure, viscosity
  !! and level-set function at time t=0. 
  subroutine get_initial_twophase(u, v, w, p, phi, visc, &
    l, m, n, dx, dy, dz, Reynolds, dpdl, mu_minus, mu_plus, &
    height, smooth_width_scale, smooth_width, Pe, epn,&
    idm)
    implicit none
    integer, intent(in) :: l !< Number of domain grid points in X(l) dimension.
    integer, intent(in) :: m !< Number of domain grid points in Y(m) dimension.
    integer, intent(in) :: n !< Number of domain grid points in Z(n) dimension.
    double precision, intent(inout) :: u(0:l,0:m,0:n-2)  !< Velocity in X(l) dimension.
    double precision, intent(inout) :: v(0:l,0:m,0:n-2)  !< Velocity in Y(m) dimension.
    double precision, intent(inout) :: w(0:l,0:m,0:n-1)  !< Velocity in Z(n) dimension.
    double precision, intent(inout) :: p(0:l,0:m,0:n)    !< Pressure.
    double precision, intent(inout) :: phi(0:l,0:m,0:n)  !< Level-set function.
    double precision, intent(inout) :: visc(0:l,0:m,0:n) !< Viscosity.
    double precision, intent(out) :: dx !< dx
    double precision, intent(out) :: dy !< dy
    double precision, intent(out) :: dz !< dz
    double precision, intent(in)  :: Reynolds     !< Reynolds number.
    double precision, intent(in)  :: dpdl         !< Pressure gradient.
    double precision, intent(in)  :: mu_minus     !< Viscosity of the lower fluid.
    double precision, intent(in)  :: mu_plus      !< Viscosity of the upper fluid.
    double precision, intent(in)  :: height       !< Interface height.
    double precision, intent(in)  :: smooth_width_scale !< Smooth width scale.
    double precision, intent(out) :: smooth_width !< Smooth width.
    double precision, intent(out) :: Pe  !< PE (DIM)
    double precision, intent(out) :: epn !< epn (DIM)
    character(len=5), intent(in)  :: idm !< Interface detection method.

    integer :: i, j, k
    integer :: ck1, ck2
    integer :: iteration_init
    double precision :: pi
    parameter(pi=3.14159265359)
    double precision :: rand_phase(1:8,1:8), rand_phase_val
    double precision :: sum_phase
    double precision :: Lx, Ly, Lz, x_val, y_val, z_val
    double precision :: phi_val, heaviside_phi
    double precision :: u2_oned(0:n-2), visc_oned(0:n)
    double precision :: u_minusz, u_plusz, mu_minusk, mu_plusk
    double precision :: RHS, relax, diag_val, residual,temp
    integer :: num_iterations
    parameter(num_iterations=10000)

    ! Initialise parameters.
    dz=1.d0/dble(n-1)
    dx=dz
    dy=dz
    smooth_width=smooth_width_scale*dx
    ! DIM value.
    epn=0.5*dz
    Pe=1.d0/(epn*epn)

    Lz=1.d0
    Lx=dx*(l-1)
    Ly=dy*(m-1)

    ! Create initial velocities and pressure.
    do k=0,n
      z_val=dz*k-(dz/2.d0)
      phi_val=z_val-height
      if(phi_val.lt.-smooth_width)then
        heaviside_phi=0.d0
      elseif(phi_val.gt.smooth_width)then
        heaviside_phi=1.d0
      else
        heaviside_phi=0.5d0+(phi_val/(2.d0*smooth_width))+ &
          (1.d0/(2.d0*pi))*sin(pi*phi_val/smooth_width)
      end if
      ! temp=(1.d0/mu_minus)+((1.d0/mu_plus)-(1.d0/mu_minus))*heaviside_phi
      ! visc_oned(k)=(1.d0/temp)/Reynolds
      temp=mu_plus*heaviside_phi+mu_minus*(1.d0-heaviside_phi)
      visc_oned(k)=temp/Reynolds
    end do ! k

    ! Create numerical base state.
    ! Solves (d/dz)[mu(d/dz)u]-dpdx=0, where we invert the diffusion operator.
    u2_oned=0.d0
    RHS=dpdl
    relax=1.2d0

    do iteration_init=1,num_iterations
      do k=0,n-2
        if(k.eq.0) then
          u_minusz=-2.d0*u2_oned(0)+u2_oned(1)/3.d0
        else
          u_minusz=u2_oned(k-1)
        end if
        if(k.eq.(n-2)) then
          u_plusz=-2.d0*u2_oned(n-2)+u2_oned(n-3)/3.d0
        else
          u_plusz=u2_oned(k+1)
        end if
        mu_minusk=(visc_oned(k)+visc_oned(k+1))/2.d0
        mu_plusk=(visc_oned(k+2)+visc_oned(k+1))/2.d0
        diag_val=(mu_plusk+mu_minusk)/(dz*dz)
        ! residual = (1.d0/diag_val) * ( ((mu_plusk*u_plusz+mu_minusk*u_minusz)/(dz*dz)) -RHS)-u2_oned(k)
        ! u2_oned(k) = u2_oned(k)+relax*residual       
        residual = (1.d0/diag_val) * (((mu_plusk*u_plusz+mu_minusk*u_minusz)/(dz*dz)) - RHS)
        u2_oned(k) = residual    
      end do ! k
    end do ! iteration_init

    do k=0,n-2
      do j=0,m
        do i=0,l
          x_val=i*dx
          y_val=(j+0.5d0)*dy
          u(i,j,k)=u2_oned(k)
          v(i,j,k)=0.d0
        end do ! i
      end do ! j
    end do ! k
  
    do k=0,n-1
      do j=0,m
        do i=0,l
          x_val=i*dx+(dx/2.d0)
          y_val=(j+0.5d0)*dy
          w(i,j,k)=0.d0
        end do ! i
      end do ! j
    end do ! k
  
    do k=0,n
      do j=0,m
        do i=0,l
          z_val=k*dz-(dz/2.d0)
          x_val=i*dx-(dx/2.d0)
          y_val=j*dx-(dy/2.d0)
          p(i,j,k)=-x_val
        end do ! i
      end do ! j
    end do ! k
    
    do k=0,n
      do j=0,m
        do i=0,l
        visc(i,j,k)=visc_oned(k)
        end do ! i
      end do ! j
    end do ! k

    ! Generate arrays of random numbers.
    ! Initialise random number generator to default state for deterministic 
    ! results.
    call random_seed()
    do ck2=1,8
      do ck1=1,8
        call random_number(rand_phase_val) 
        rand_phase(ck1,ck2)=2.d0*pi*rand_phase_val
      end do ! ck1
    end do ! ck2
  
    do k=0,n
      do j=0,m
        do i=0,l
          z_val=k*dz-(dz/2.d0)
          x_val=i*dx-(dx/2.d0)
          y_val=j*dx-(dy/2.d0)
          sum_phase=0.d0
          do ck2=1,8
            do ck1=1,8
              sum_phase=sum_phase+ &
                cos(ck1*(2.d0*pi/Lx)*x_val+(ck2-1)*(2.d0*pi/Ly)*y_val+rand_phase(ck1,ck2))
            end do ! ck1
          end do ! ck2   
          sum_phase=(0.02/9.d0)*sum_phase
          if (idm==level_set_method) then 
            ! LSM
            phi(i,j,k)=z_val-(height+sum_phase)
          else 
            ! DIM
             temp=z_val-(height+sum_phase)
             phi(i,j,k)=0.5d0+0.5d0*tanh(temp/(2.d0*(2d0**0.5d0)*epn))
          end if
        end do ! i
      end do ! j
    end do ! k
    return
  end subroutine get_initial_twophase

end module twophase_initialisation_wave
