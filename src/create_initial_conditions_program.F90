!> Create TPLS initial conditions and save as TPLS input files.
!!
!! Usage: running as follows displays available options.
!! \code
!! $ ./create_initial_conditions [-f INPUT_OPTIONS_FILE] [-d OUTPUT_DIRECTORY]
!! \endcode
!!
!! TPLS initial conditions are created then saved as the following
!! TPLS input files:
!! - initial_u.dat - Velocity in X(l) dimension.
!! - initial_v.dat - Velocity in Y(m) dimension.
!! - initial_w.dat - Velocity in Z(n) dimension.
!! - initial_pressure.dat - Pressure.
!! - initial_phi.dat - Level-set function.
!! - initial_viscosity.dat - Viscosity.
!! - initial_config.dat - TPLS configuration used to create the above.
!! 
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune,
!! Toni Collis, David Scott, Peter Spelt, Mike Jackson.
!! @version $Revision: 328 $
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt, The University of
!! Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

program create_initial_conditions_program

  use option_names
  use options_utils
  use tpls_configure
  use tpls_error_check
  use tpls_io
  use twophase_initialisation_wave

  implicit none

  ! Configuration options file.
  character(len=60) :: default_filename
  parameter(default_filename='initial_config.opt')
  character(len=256) :: filename
  integer :: num_options
  character(len=80), dimension(:,:), allocatable :: options

  ! Initial conditions configuration directory.
  character(len=60) :: default_config_dir
  parameter(default_config_dir='.')
  character(len=256) :: config_dir

  ! Initial conditions.
  double precision,allocatable,dimension(:,:,:) :: u,v,w
  double precision,allocatable,dimension(:,:,:) :: pressure,phi,viscosity
  integer          :: maxl,maxm,maxn
  double precision :: Re,dpdl,mu_minus,mu_plus,smooth_width_scale,height
  double precision :: smooth_width,scap,omega_max
  double precision :: dx,dy,dz,dt
  character(len=5) :: idm

  ! DIM-specific initial conditions.
  double precision :: epn, Pe 

  ! Additional options to save.
  integer :: num_header_lines
  character(len=256),dimension(:),allocatable :: header
  integer :: num_initial_conditions
  character(len=80), dimension(:,:), allocatable :: initial_conditions

  ! General purpose.
  integer :: stdout
  parameter(stdout=6)
  character(len=256) :: message
  integer :: ierr
  logical :: file_exists

  write(*,'(A,I6,A8)') 'TPLS initial conditions (version $Revision: 328 $)'

  ! Read command-line options and options file.

  call get_configuration(default_filename,default_config_dir,config_dir,&
       options,num_options,ierr,message)
  if (ierr/=0) call abort(ierr, message)

  ! Read options.

  call get_grid_configuration(options,num_options,maxl,maxm,maxn,&
       ierr,message,.true.)
  if (ierr/=0) call abort(ierr, message)
  call get_fluid_configuration(options,num_options,Re,dpdl,mu_minus,mu_plus,&
       height,scap,omega_max,smooth_width_scale,.true.,ierr,message,.true.)
  if (ierr/=0) call abort(ierr, message)
  call get_timestep(options,num_options,dt,ierr,message,.true.)
  if (ierr/=0) call abort(ierr, message)
  call get_interface_detection_method(options,num_options,idm,&
       ierr,message,.true.)
  if (ierr/=0) call abort(ierr, message)

  ! Create initial conditions.

  allocate(u(0:maxl,0:maxm,0:maxn-2))
  allocate(v(0:maxl,0:maxm,0:maxn-2))
  allocate(w(0:maxl,0:maxm,0:maxn-1))
  allocate(pressure(0:maxl,0:maxm,0:maxn))
  allocate(phi(0:maxl,0:maxm,0:maxn)) 
  allocate(viscosity(0:maxl,0:maxm,0:maxn))

  call get_initial_twophase(u,v,w,pressure,phi,viscosity,&
       maxl,maxm,maxn,dx,dy,dz,&
       Re,dpdl,mu_minus,mu_plus,height,smooth_width_scale,smooth_width,&
       Pe,epn,idm)

  ! Save initial conditions.

  write(*,'(A,A)') 'Saving initial configuration files to: ',trim(config_dir)
  call get_filepath(config_dir,initial_u_file,filename)
  write(*,'(A)') trim(filename)
  call save_data(filename,'Initial velocity (u)',u,0,maxl,0,maxm,0,maxn-2)
  call get_filepath(config_dir,initial_v_file,filename)
  write(*,'(A)') trim(filename)
  call save_data(filename,'Initial velocity (v)',v,0,maxl,0,maxm,0,maxn-2)
  call get_filepath(config_dir,initial_w_file,filename)
  write(*,'(A)') trim(filename)
  call save_data(filename,'Initial velocity (w)',w,0,maxl,0,maxm,0,maxn-1)
  call get_filepath(config_dir,initial_pressure_file,filename)
  write(*,'(A)') trim(filename)
  call save_data(filename,'Initial pressure',pressure,0,maxl,0,maxm,0,maxn)
  call get_filepath(config_dir,initial_phi_file,filename)
  write(*,'(A)') trim(filename)
  call save_data(filename,'Initial level-set (phi)',phi,0,maxl,0,maxm,0,maxn)
  call get_filepath(config_dir,initial_viscosity_file,filename)
  write(*,'(A)') trim(filename)
  call save_data(filename,'Initial viscosity',viscosity,0,maxl,0,maxm,0,maxn)

  ! Create initial_conditions with values created by get_initial_twophase.

  allocate(initial_conditions(20,2))
  num_initial_conditions=0
  call put_typed_option(initial_conditions,num_initial_conditions,opt_maxl,maxl)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_maxm,maxm)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_maxn,maxn)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_re,Re)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_dpdl,dpdl)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_mu_plus,mu_plus)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_mu_minus,mu_minus)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_height,height)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_scap,scap)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_omega_max,omega_max)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_smooth_width,smooth_width)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_dx,dx)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_dy,dy)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_dz,dz)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_dt,dt)
  call put_typed_option(initial_conditions,num_initial_conditions,opt_idm,idm)
  if (idm==diffuse_interface_method) then
     call put_typed_option(initial_conditions,num_initial_conditions,opt_pe,Pe)
     call put_typed_option(initial_conditions,num_initial_conditions,opt_epn,epn)
  end if

  ! Create header.

  num_header_lines=2
  allocate(header(2))
  header(1)='TPLS initial conditions'
  header(2)='Auto-generated by tpls_io $Revision: 328 $'

  ! Save initial conditions,

  call get_filepath(config_dir,initial_config_file,filename)
  write(*,'(A)') trim(filename)
  call save_options(filename,initial_conditions,num_initial_conditions,header,2)

  ! Clean-up.

  deallocate(u)
  deallocate(v)
  deallocate(w)
  deallocate(pressure)
  deallocate(phi)
  deallocate(viscosity)
  deallocate(options)
  deallocate(header)
  deallocate(initial_conditions)

end program create_initial_conditions_program
