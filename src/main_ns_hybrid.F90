!------------------------------------------------------------!
! TPLS version 2.0                                           !
! main_ns_hybrid.F90                                         !
! Copyright (c) 2013-2015, Prashant Valluri,                 !
!                          Lennon O Naraigh, Iain Bethune,   !
!                          Toni Collis, David Scott,         !
!                          Peter Spelt,                      !
!                          The University of Edinburgh.      !
! All rights reserved.                                       !
!                                                            !
! This program is distributed under the BSD License          !
! See LICENSE.txt for details                                !
!------------------------------------------------------------!

program tpls_program

  use petsc

  use advect_phi
  use cahn_hilliard_solver 
  use grid_utils
  use jacobi_iteration_allflux
  use levelset
  use momentum_allflux
  use option_names
  use options_utils
  use pressure
  use pressure_solver
  use sor_iteration_allflux
  use tpls_configure
  use tpls_configure_mpi_petsc
  use tpls_io
  use tpls_mpi
  use tpls_mpi_error_check
  use tpls_state
  use two_phase_levelset

  implicit none

#include "finclude/petscdef.h"

  ! Included from tpls_mpi
  ! integer :: master_id
  ! parameter(master_id=0)

  ! Included from tpls_state
  ! PetscInt :: maxl, maxm, maxn
  ! Note that the first two dimensions are swapped so that as
  ! far as PETSc is concerned  it is the first dimension that is
  ! periodic.  
  ! PetscInt :: global_dim_x, global_dim_y, global_dim_z
  ! PetscScalar, allocatable, dimension(:, :, :) :: pres
  ! PetscScalar :: dx, dy, dz, dt
  ! PetscScalar, allocatable, dimension(:, :, :) :: RHS_p, u3
  ! PetscScalar, allocatable, dimension(:) :: u_inlet

  ! Configuration options file.
  character(len=60) :: default_filename
  parameter(default_filename='tpls_config.opt')
  character(len=256) :: filename
  integer :: num_options
  character(len=80), dimension(:,:), allocatable :: options

  ! Initial conditions configuration directory.
  character(len=60) :: default_config_dir
  parameter(default_config_dir='.')
  character(len=256) :: config_dir

  ! Initial conditions configuration options.
  integer :: num_initial_conditions
  character(len=80), dimension(:,:), allocatable :: initial_conditions

  ! Process grid bounds.
  integer :: my_id 
  integer :: num_procs
  integer :: num_procs_x, num_procs_y ! User configuration / best guess.
  integer :: num_procs_z
  parameter (num_procs_z=1)
  integer :: Ndim
  parameter (Ndim=3)
  integer :: etag,stag
  parameter (etag=1,stag=2)
  logical :: auto_decomp

  ! User configuration.

  ! Momentum equation solver configuration.
  integer :: max_iteration_mom_u, max_iteration_mom_v, max_iteration_mom_w
  ! Level-set equation solver configuration.
  integer :: max_iteration_levelset
  ! DIM equation solver configuration.
  integer :: max_iteration_dim
  ! phi dat file output frequency.
  integer :: phi_dat_frequency
  ! uvw dat file output frequency.
  integer :: uvw_dat_frequency
  ! uvw dat file output frequency.
  integer :: backup_frequency
  ! Backup file types - true if to back up to this file type
  logical :: backuptext
  logical :: backuphdf5
  ! Number of timesteps.
  integer :: n_timesteps
  ! Fluid flow.
  double precision :: maxu
  ! Interface detection method.
  character(len=5) :: idm

  ! Initial conditions configuration.

  ! PetscInt :: maxl, maxm, maxn
  ! PetscScalar :: dx, dy, dz, dt
  double precision :: Re,dpdl
  double precision :: mu_minus, mu_plus
  double precision :: smooth_width, height, scap, omega_max
  double precision :: epn, Pe ! DIM
  double precision, allocatable, dimension(:,:,:) :: u_global
  double precision, allocatable, dimension(:,:,:) :: v_global
  double precision, allocatable, dimension(:,:,:) :: w_global
  double precision, allocatable, dimension(:,:,:) :: pres_global
  double precision, allocatable, dimension(:,:,:) :: phi_global
  double precision, allocatable, dimension(:,:,:) :: visc_global

  ! MPI
  integer :: status(MPI_STATUS_SIZE)
  integer :: dims(Ndim), coords(Ndim)
  integer :: comm2d_quasiperiodic
  logical :: periodic(Ndim)
  logical :: reorder
  integer :: neighbours(6)
  integer :: sender_id
  integer :: n_local_x, n_local_y, n_local_z_uv, n_local_z_w, n_local_z_p
  integer :: sx, sy, ex, ey
  integer :: sz_uv, ez_uv, sz_w, ez_w, sz_p, ez_p
  integer :: sx_o, ex_o, sy_o, ey_o
  integer :: sx_f, ex_f, sy_f, ey_f
  integer :: stride_uv_xz, stride_uv_yz
  integer :: stride_w_xz,  stride_w_yz
  integer :: stride_p_xz,  stride_p_yz
  integer :: stride_p_aug1_xz, stride_p_aug1_yz
  integer :: stride_p_aug2_xz, stride_p_aug2_yz
  integer :: stride_p_augaug1_xz, stride_p_augaug1_yz
  integer :: stride_p_augaug2_xz, stride_p_augaug2_yz
  integer :: stride_p_augaug3_xz, stride_p_augaug3_yz
  integer :: px, py, pz

  ! TPLS local
  double precision, allocatable, dimension(:,:,:) :: u2,v2,w2
  double precision, allocatable, dimension(:,:,:) :: v3,w3
  double precision, allocatable, dimension(:,:,:) :: u3_old,v3_old,w3_old
  double precision, allocatable, dimension(:,:,:) :: conv0_u,conv0_v,conv0_w
  double precision, allocatable, dimension(:,:,:) :: conv1_u,conv1_v,conv1_w
  double precision, allocatable, dimension(:,:,:) :: conv2_u,conv2_v,conv2_w
  double precision, allocatable, dimension(:,:,:) :: RHS_u,RHS_v,RHS_w
  double precision, allocatable, dimension(:,:,:) :: tempr_u,tempr_v,tempr_w
  double precision, allocatable, dimension(:,:,:) :: diffusion_u,diffusion_v,diffusion_w
  double precision, allocatable, dimension(:,:,:) :: pres_old,tempr_p
  double precision, allocatable, dimension(:,:,:) :: output_u
  double precision, allocatable, dimension(:,:,:) :: output_w
  double precision, allocatable, dimension(:,:,:) :: output_p
  double precision, allocatable, dimension(:,:,:) :: csf_u2,csf_v2,csf_w2
  double precision, allocatable, dimension(:,:,:) :: csf_u1,csf_v1,csf_w1
  double precision, allocatable, dimension(:,:,:) :: csf_u0,csf_v0,csf_w0
  double precision, allocatable, dimension(:,:,:) :: u2_aug,v2_aug,w2_aug
  double precision, allocatable, dimension(:,:,:) :: viscosity1,viscosity2,viscosity2_aug,tempr_visc
  double precision, allocatable, dimension(:,:,:) :: phi2,phi3,tempr_phi
  double precision, allocatable, dimension(:,:,:) :: conv2_phi,conv1_phi,conv0_phi
  double precision, allocatable, dimension(:,:,:) :: ConvPrevPhi,sign_mx
  double precision, allocatable, dimension(:,:,:) :: fx_csf,fy_csf,fz_csf
  double precision, allocatable, dimension(:,:,:) :: conv0_u_global,conv1_u_global
  double precision, allocatable, dimension(:,:,:) :: conv0_v_global,conv1_v_global
  double precision, allocatable, dimension(:,:,:) :: conv0_w_global,conv1_w_global
  double precision, allocatable, dimension(:,:,:) :: conv0_phi_global,conv1_phi_global
  double precision, allocatable, dimension(:,:,:) :: pres_old_global
  double precision, allocatable, dimension(:,:,:) :: phi_old ! DIM
  double precision :: z_val
  double precision :: err1,err2
  double precision :: cfl
  double precision :: heaviside_phi,temp_visc,phi_val
  double precision :: pi
  parameter (pi=3.14159265359)

  integer :: i,j,k,ex_max

  ! PETSc
  DM :: da
  KSP :: ksp
  PetscInt, allocatable, dimension(:)  :: petsc_y_ranges
  PetscInt :: its
  PetscReal :: rnorm

  ! TPLS operation 
  integer :: iteration
  integer :: iteration_sor
  integer :: backup_counter
  double precision :: t0,t1,t2,t3
  double precision :: io_t0,io_t1, io_t2,io_t3,io_t4,io_t5   
  double precision :: t_init,t_iter,t_iter_rate,t_fina,t_total
  double precision :: t_text_bak,t_text_reg,t_text_comms,t_text_tot
  double precision :: t_netcdf_reg,t_netcdf_bak,t_netcdf_tot
  double precision :: text_back_time,text_reg_time,text_comms_time
  double precision :: hdf5_back_time, hdf5_reg_time
  double precision :: total_text_time, total_hdf5_time

  ! General purpose
  integer :: stdout
  parameter(stdout=6)
  character(len=256) :: message
  integer :: ierr
  logical :: file_exists
  character(len=256) ::strg
  integer :: count
  integer :: backup

  t0 = mpi_wtime() 

  ! ******************************************************************************************

  ! Default parameter values.

  backup_counter=0

  ! ****************************************************************************************

  ! Initialisation of petsc
  call initialise(ierr)

  ! Set up MPI.
  ! After set up master_id will be responsible for printing
  ! configuration options or configuration errors. Other processes
  ! will abort silently.  
  call get_size_and_rank(ierr, my_id, num_procs)

  write(*,'(A,I6,A8)') &
       'TPLS (version $Revision: 328 $) process ', my_id, ' started'
  if (is_master(my_id)) then
     write(*,'(A,I6)') 'Number of processes available: ', num_procs
  end if
  flush stdout

  ! Read command-line options and options file.

  call get_configuration(default_filename,default_config_dir,config_dir,&
       options,num_options,ierr,message)
  if (ierr/=0) call abort(ierr, message, my_id)

  ! Load initial conditions configuration file.

  call get_filepath(config_dir,initial_config_file,filename)
  write(*,'(A,A)') 'Loading initial conditions configuration file: ',trim(filename)
  call load_options(filename,initial_conditions,num_initial_conditions,ierr)
  if (ierr==option_file_not_found_err) then
     write(message,'(A,A)') 'Error: Initial conditions configuration file not found: ',trim(filename)
     call abort(config_file_not_found_err,message,my_id)
  else if (ierr/=0) then
     write(message,'(A,A)') 'Error: Loading initial conditions configurations file: ',trim(filename)
     call abort(ierr,message,my_id)
  end if

  ! Read initial conditions configuration.

  call get_grid_configuration(initial_conditions,num_initial_conditions,&
       maxl,maxm,maxn,ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  call get_fluid_configuration(initial_conditions,num_initial_conditions,&
       Re,dpdl,mu_minus,mu_plus,height,scap,omega_max,smooth_width,.false.,&
       ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  call get_dgrid_configuration(initial_conditions,num_initial_conditions,&
       dx,dy,dz,ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  call get_timestep(initial_conditions,num_initial_conditions,dt,&
       ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  call get_interface_detection_method(initial_conditions,&
       num_initial_conditions,idm,ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  if (idm==diffuse_interface_method) then
     call get_dim_configuration(initial_conditions,&
          num_initial_conditions,Pe,epn,ierr,message,is_master(my_id))
     if (ierr/=0) call abort(ierr,message,my_id)
  end if
  deallocate(initial_conditions)

  ! Read configuration.

  call get_process_configuration(options,num_options,&
       auto_decomp,num_procs,num_procs_x,num_procs_y,num_procs_z,&
       maxl,maxm,maxn,ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  call get_solver_configuration(options,num_options,&
       max_iteration_mom_u,max_iteration_mom_v,max_iteration_mom_w, &
       max_iteration_levelset,n_timesteps,ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  call get_output_configuration(options,num_options,backuptext, backuphdf5,&
       phi_dat_frequency,uvw_dat_frequency,backup_frequency, &
       ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  call get_misc_configuration(options,num_options,maxu,&
       ierr,message,is_master(my_id))
  if (ierr/=0) call abort(ierr,message,my_id)
  if (idm==diffuse_interface_method) then
     call get_dim_solver_configuration(options,num_options,&
          max_iteration_dim,ierr,message,is_master(my_id))
     if (ierr/=0) call abort(ierr,message,my_id)
  end if

  deallocate(options)

  flush stdout

  ! Finish initialisation based on configuration options.

  ! Note that the first two dimensions are swapped so that as
  ! far as PETSc is concerned  it is the first dimension that is
  ! periodic.  
  global_dim_x = maxm-1
  global_dim_y = maxl-1
  global_dim_z = maxn-1
  allocate(u_inlet(0:global_dim_z-1))

  allocate(petsc_y_ranges(0:num_procs_x-1))
  cfl=maxu*dt/dx

  allocate(u_global(0:maxl,0:maxm,0:maxn-2))
  allocate(v_global(0:maxl,0:maxm,0:maxn-2))
  allocate(w_global(0:maxl,0:maxm,0:maxn-1))
  allocate(pres_global(0:maxl,0:maxm,0:maxn))
  allocate(phi_global(0:maxl,0:maxm,0:maxn)) 
  allocate(visc_global(0:maxl,0:maxm,0:maxn))

  allocate(conv0_u_global(0:maxl,0:maxm,0:maxn-2))
  allocate(conv1_u_global(0:maxl,0:maxm,0:maxn-2))
  allocate(conv0_v_global(0:maxl,0:maxm,0:maxn-2))
  allocate(conv1_v_global(0:maxl,0:maxm,0:maxn-2))
  allocate(conv0_w_global(0:maxl,0:maxm,0:maxn-1))
  allocate(conv1_w_global(0:maxl,0:maxm,0:maxn-1))
  allocate(conv0_phi_global(0:maxl,0:maxm,0:maxn))
  allocate(conv1_phi_global(0:maxl,0:maxm,0:maxn))
  allocate(pres_old_global(0:maxl,0:maxm,0:maxn))

  dims(1) = num_procs_x
  dims(2) = num_procs_y
  dims(3) = num_procs_z

  periodic(1) = .false.
  periodic(2) = .true.
  periodic(3) = .false.
  reorder = .false.

  call setup_MPI(Ndim,dims,comm2d_quasiperiodic,&
       my_id,neighbours, &
       maxl, maxm, maxn, &
       global_dim_x, global_dim_y, global_dim_z, &
       sx,ex,sy,ey, &
       n_local_x,n_local_y,coords, &
       sz_uv, ez_uv, sz_w, ez_w, sz_p, ez_p, &
       n_local_z_uv, n_local_z_w, n_local_z_p, &
       petsc_y_ranges, ksp,da)

  call compute_initial_strides(ex,ey,sx,sy,ex_max,my_id, &
       stride_uv_xz, stride_uv_yz,stride_w_xz,  stride_w_yz,stride_p_xz,  stride_p_yz, &
       sz_uv, ez_uv, sz_w, ez_w, sz_p, ez_p,stride_p_aug1_yz, stride_p_aug1_xz, &
       stride_p_aug2_yz, stride_p_aug2_xz,stride_p_augaug1_xz,stride_p_augaug1_yz, &
       stride_p_augaug2_xz,stride_p_augaug2_yz,stride_p_augaug3_xz,stride_p_augaug3_yz )

  ! ******************************************************************************************
  ! Compute initial value of all variables

  if (is_master(my_id)) then
     write(*,*) 'dt=',dt
     write(*,*) 'cfl=',cfl
     write(*,*) 'n_timesteps=', n_timesteps
     write(*,*) '(dx,dy,dz)=',dx,dy,dz
  end if

  if (is_master(my_id)) then
     write (*,'(A)') 'Loading initial values'
     call load_configuration_data(config_dir,initial_u_file,maxl,maxm,maxn-2,u_global,ierr,message)
     if (ierr/=0) call abort(ierr,message,my_id)
     call load_configuration_data(config_dir,initial_v_file,maxl,maxm,maxn-2,v_global,ierr,message)
     if (ierr/=0) call abort(ierr,message,my_id)
     call load_configuration_data(config_dir,initial_w_file,maxl,maxm,maxn-1,w_global,ierr,message)
     if (ierr/=0) call abort(ierr,message,my_id)
     call load_configuration_data(config_dir,initial_pressure_file,maxl,maxm,maxn,pres_global,ierr,message)
     if (ierr/=0) call abort(ierr,message,my_id)
     call load_configuration_data(config_dir,initial_phi_file,maxl,maxm,maxn,phi_global,ierr,message)
     if (ierr/=0) call abort(ierr,message,my_id)
     call load_configuration_data(config_dir,initial_viscosity_file,maxl,maxm,maxn,visc_global,ierr,message)
     if (ierr/=0) call abort(ierr,message,my_id)

     ! The subroutine "channel_output" is redundant.
     write(filename,'(A,A)')'AA_initial_phi','.dat'
     write(*,*) 'Outputing initial phi'
     call channel_output_phi(phi_global,maxl,maxm,maxn,dx,dy,dz,filename)
  end if

  call mpi_bcast(u_global, (maxl+1)*(maxm+1)*(maxn-1),mpi_double_precision,0,PETSC_COMM_WORLD,ierr)
  if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_bcast in intiial compute on u_global',my_id)
  call mpi_bcast(v_global, (maxl+1)*(maxm+1)*(maxn-1),mpi_double_precision,0,PETSC_COMM_WORLD,ierr)
  if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_bcast in intiial compute on v_global',my_id)
  call mpi_bcast(w_global, (maxl+1)*(maxm+1)*(maxn  ),mpi_double_precision,0,PETSC_COMM_WORLD,ierr)
  if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_bcast in intiial compute on w_global',my_id)
  call mpi_bcast(pres_global,(maxl+1)*(maxm+1)*(maxn+1),mpi_double_precision,0,PETSC_COMM_WORLD,ierr)
  if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_bcast in intiial compute on pres_global',my_id)
  call mpi_bcast(phi_global,(maxl+1)*(maxm+1)*(maxn+1),mpi_double_precision,0,PETSC_COMM_WORLD,ierr)
  if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_bcast in intiial compute on phi_global',my_id)

  u_inlet(:)=u_global(0,0,:)

  ! ********************* Initial conditions for pressure 

  allocate(      pres(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),       &
       pres_old(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),       &
       RHS_p(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),       &
       tempr_p(1:n_local_x,1:n_local_y,1:n_local_z_p), &
       output_p(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)         )

  pres(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)         =pres_global(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
  pres_old=0.d0*pres

  ! ********************* Initial conditions for phi and various other things

  allocate(       phi2(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),       phi3(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),   &
       phi_old(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),                                           &
       conv2_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),  conv1_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),   &
       conv0_phi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p), ConvPrevPhi(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),  &
       fx_csf(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),     fy_csf(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),   &
       fz_csf(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),                                              &
       tempr_phi(1:n_local_x,1:n_local_y,1:n_local_z_p))

  allocate(       sign_mx(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)    )

  allocate(     viscosity1(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),   &
       viscosity2(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p),   &
       viscosity2_aug(sx-2:ex+2,sy-2:ey+2,sz_p:ez_p),   &
       tempr_visc(1:n_local_x,1:n_local_y,1:n_local_z_p))

  phi2(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)         =phi_global(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
  phi3=0.d0
  phi_old = phi2

  conv2_phi=0.d0
  conv1_phi=0.d0
  conv0_phi=0.d0
  ConvPrevPhi = 0.d0
  sign_mx=0.d0

  viscosity1(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)=visc_global(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
  viscosity2(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)=visc_global(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)
  viscosity2_aug(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)=visc_global(sx-1:ex+1,sy-1:ey+1,sz_p:ez_p)

  ! ********************* Initial conditions in x(u) direction 

  allocate(          u2(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),       u3(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),    &
       u3_old(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),                                               &
       conv0_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),  conv1_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),    &
       conv2_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),                                               &
       csf_u0(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),                                                &
       csf_u1(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),    csf_u2(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),    &
       diffusion_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),    RHS_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),    &
       tempr_u(1:n_local_x,1:n_local_y,1:n_local_z_uv),                                            &
       output_u(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)                                                   )

  ! Augmented array          
  allocate(u2_aug(sx-2:ex+2,sy-2:ey+2,sz_uv:ez_uv) )

  u2(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)    =u_global(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
  u3(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)    =0.d0*u2
  u3_old(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)=0.d0*u2

  conv0_u=0.d0
  conv1_u=0.d0
  conv2_u=0.d0

  u2_aug(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)=u_global(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)

  diffusion_u=0.d0
  RHS_u=   0.d0

  csf_u0=0.d0
  csf_u1=0.d0
  csf_u2=0.d0

  ! ********************* Initial conditions in y(v) direction 

  allocate(          v2(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),     v3(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),      &
       v3_old(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),                                               &
       conv0_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv), conv1_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),     &
       conv2_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),                                               &
       csf_v0(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),                                                &
       csf_v1(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),    csf_v2(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),    &
       diffusion_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),   RHS_v(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv),     &
       tempr_v(1:n_local_x,1:n_local_y,1:n_local_z_uv))

  ! Augmented array            
  allocate(v2_aug(sx-2:ex+2,sy-2:ey+2,sz_uv:ez_uv) )

  v2(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)    =v_global(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)
  v3(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)    =0.d0*v2
  v3_old(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)=0.d0*v2

  conv0_v=0.d0
  conv1_v=0.d0
  conv2_v=0.d0

  v2_aug(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)=v_global(sx-1:ex+1,sy-1:ey+1,sz_uv:ez_uv)

  diffusion_v=0.d0
  RHS_v=   0.d0

  csf_v0=0.d0
  csf_v1=0.d0
  csf_v2=0.d0

  ! ********************* Initial conditions in z(w) direction 

  allocate(          w2(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),      w3(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),       &
       w3_old(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),                                               &
       conv0_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w), conv1_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),       &
       conv2_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),                                               &
       csf_w0(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),                                                &
       csf_w1(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),    csf_w2(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),      &
       diffusion_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),   RHS_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w),       &
       tempr_w(1:n_local_x,1:n_local_y,1:n_local_z_w),                                           &
       output_w(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)                                                   )

  ! Augmented array             
  allocate(w2_aug(sx-2:ex+2,sy-2:ey+2,sz_w:ez_w) )

  w2(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)    =w_global(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)
  w3(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)    =0.d0*w2
  w3_old(sx-1:ex+1,sy-1:ey+1,sz_w:ez_w)=0.d0*w2

  conv0_w=0.d0
  conv1_w=0.d0
  conv2_w=0.d0

  w2_aug=0.d0

  diffusion_w=0.d0
  RHS_w=   0.d0

  csf_w0=0.d0
  csf_w1=0.d0
  csf_w2=0.d0

  ! ****************************************************************************************  
  ! Getting coords of base node in Cartesian topology.

  If (is_master(my_id)) then
     pz=0
     do py=0,dims(2)-1
        do px=0,dims(1)-1

           if((pz==0).and.(py==0).and.(px==0))then
              sx_o=sx
              ex_o=ex
              sy_o=sy
              ey_o=ey 
           end if

        end do
     end do
  end if

  ! *************************************************************************************************
  !       Time loop
  ! *************************************************************************************************

  !Set timers to zero
  text_back_time=0
  text_reg_time=0
  text_comms_time=0
  hdf5_back_time=0
  hdf5_reg_time=0
  total_text_time=0

  t1 = mpi_wtime()
  do iteration=1,n_timesteps

     call do_advect_phi(phi2, u2, v2, w2, conv2_phi,         &
          sx, ex, sy, ey, sz_p, ez_p, sz_uv, sz_w, ez_uv, ez_w, &
          maxl, maxn, ex_max, neighbours, comm2d_quasiperiodic, &
          dx, dy, dz,                                           &
          stride_p_augaug1_xz, stride_p_augaug1_yz,             &
          stride_p_augaug2_xz, stride_p_augaug2_yz,             &
          stride_p_augaug3_xz, stride_p_augaug3_yz,             &
          stride_uv_xz, stride_uv_yz,                           &
          stride_w_xz, stride_w_yz)

     if (idm==level_set_method) then
        call do_tpls(phi2, phi3, conv0_phi, conv1_phi, conv2_phi, sx, ex, sy, ey,             &
             sz_p, ez_p, maxn, ex_max, stride_p_xz, stride_p_yz,                      &
             max_iteration_levelset, neighbours, comm2d_quasiperiodic, height, err2,  &
             dx, dz, dt)

        CALL MPI_ALLREDUCE(err2, err1,  1, MPI_DOUBLE_PRECISION,   MPI_SUM,PETSC_COMM_WORLD, ierr)
        if(ierr /= 0 ) call abort(ierr, 'Failed to call MPI_ALLREDUCE in iteration loop',my_id)
        if (is_master(my_id)) then
           Write(*,*)'Iteration=',iteration, 'phi-residual is ', err1
        end if
     else 
        ! DIM - diffuse_interface_method
        call dim(maxl, maxn, dx, dy, dz, dt, epn, Pe, phi2+0.5, phi_old+0.5, conv2_phi,       &
             ConvPrevPhi, phi3, sx, ex, sy, ey, sz_p, ez_p, height,               &
             stride_p_xz, stride_p_yz, neighbours, comm2d_quasiperiodic, max_iteration_dim)
        phi3=phi3-0.5
     endif

     ! Non-augmented phi: We need to use phi3 in the computation of the CSF,
     ! according to Kang et al.

     ! Should this exchange be inside the 'tpls' subroutine?
     call exchange2d(phi3,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
     call get_csf(fx_csf,fy_csf,fz_csf,phi3,sx,ex,sy,ey,maxn,ex_max,dx,dy,dz,scap,smooth_width)

     ! Non-augmented CSFs: exchange first-order halos only (default now is with corners). 
     call exchange2d(fx_csf,stride_p_xz,stride_p_yz, &
          neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
     call exchange2d(fy_csf,stride_p_xz,stride_p_yz, &
          neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
     call exchange2d(fz_csf,stride_p_xz,stride_p_yz, &
          neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)

     !$omp parallel workshare
     forall (i = sx:ex, j = sy:ey, k = 0:maxn-2)
        csf_u2(i,j,k)=(fx_csf(i+1,j+1,k+1)+fx_csf(i,j+1,k+1))/2.d0
     end forall

     forall (i = sx:ex, j = sy:ey, k = 0:maxn-2)
        csf_v2(i,j,k)=(fy_csf(i,j+1,k+1)+fy_csf(i,j,k+1))/2.d0
     end forall

     forall (i = sx:ex, j = sy:ey, k = 0:maxn-1)
        csf_w2(i,j,k)=(fz_csf(i,j+1,k+1)+fz_csf(i,j+1,k))/2.d0
     end forall
     !$omp end parallel workshare

     ! *****        viscosity business      **************************************************

     ! compute viscosity for lsm
     if(idm==level_set_method) then
        call get_viscosity(viscosity2,phi2,ex,ey,sx,sy,maxn,smooth_width,Re,mu_minus,mu_plus)
     else
        ! compute viscosity for dim
        call get_viscosity(viscosity2,phi2-0.5,ex,ey,sx,sy,maxn,smooth_width,Re,mu_minus,mu_plus)
     end if

     ! Initialise augmented viscosity
     !$omp parallel workshare
     viscosity2_aug(sx:ex,sy:ey,sz_p:ez_p)=viscosity2(sx:ex,sy:ey,sz_p:ez_p)
     !$omp end parallel workshare

     ! Initialise augmented viscosity at boundary points.

     if(sx==1)then
        do k=0,maxn
           do j=sy,ey
              if (idm==level_set_method) then
                 ! for lsm
                 phi_val=phi2(0,j,k)
              else
                 ! for dim
                 phi_val=phi2(0,j,k)-0.5
              endif
              if(phi_val.lt.-smooth_width)then
                 heaviside_phi=0.d0
              elseif(phi_val.gt.smooth_width)then
                 heaviside_phi=1.d0
              else
                 heaviside_phi=0.5d0+(phi_val/(2.d0*smooth_width))+(1.d0/(2.d0*pi))*sin(pi*phi_val/smooth_width)
              end if

              temp_visc=mu_plus*heaviside_phi+mu_minus*(1.d0-heaviside_phi)
              viscosity2_aug(0,j,k)=temp_visc/Re
              viscosity2(0,j,k)=temp_visc/Re
           end do
        end do
     end if

     if(ex==ex_max)then
        do k=0,maxn
           do j=sy,ey
              if (idm==level_set_method) then
                 ! for lsm
                 phi_val=phi2(ex_max+1,j,k)
              else
                 ! for dim
                 phi_val=phi2(ex_max+1,j,k)-0.5
              endif

       
              if(phi_val.lt.-smooth_width)then
                 heaviside_phi=0.d0
              elseif(phi_val.gt.smooth_width)then
                 heaviside_phi=1.d0
              else
                 heaviside_phi=0.5d0+(phi_val/(2.d0*smooth_width))+(1.d0/(2.d0*pi))*sin(pi*phi_val/smooth_width)
              end if

              temp_visc=mu_plus*heaviside_phi+mu_minus*(1.d0-heaviside_phi)
              viscosity2_aug(ex_max+1,j,k)=temp_visc/Re
              viscosity2(ex_max+1,j,k)=temp_visc/Re
           end do
        end do
     end if

     ! Exchange first-order halos
     call exchange2d_aug1(viscosity2_aug,stride_p_aug1_xz,stride_p_aug1_yz, &
          neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)

     ! Exchange second-order halos 
     call exchange2d_aug2(viscosity2_aug,stride_p_aug2_xz,stride_p_aug2_yz, &
          neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)

     ! Exchange halos for non-augmented viscosity arrays
     call exchange2d(viscosity2,stride_p_xz,stride_p_yz, &
          neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)
     call exchange2d(viscosity1,stride_p_xz,stride_p_yz, &
          neighbours,ex,ey,ez_p,sx,sy,sz_p, comm2d_quasiperiodic)

     ! ***** Convective and diffusive terms      *********************************************

     call exchange2d(u2,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv, comm2d_quasiperiodic)
     call exchange2d(v2,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv, comm2d_quasiperiodic)
     call exchange2d(w2, stride_w_xz, stride_w_yz,neighbours,ex,ey,ez_w, sx,sy,sz_w,  comm2d_quasiperiodic)

     call get_conv_all(conv2_u,conv2_v,conv2_w,viscosity2_aug,u2,v2,w2,sx,sy,ex,ey,maxn,ex_max,dx,dy,dz)
     ! diffusive term - uses flux-conservative differencing.
     call get_diffusion_all(diffusion_u,diffusion_v,diffusion_w,viscosity2_aug,u2,v2,w2, &
          sx,sy,ex,ey,maxn,ex_max,dx,dy,dz)

     !! AB2:
     !! RHS_u=u2-dt*((3.d0/2.d0)*conv2_u-(1.d0/2.d0)*conv1_u)+(dt/2.d0)*diffusion_u+dt*((3.d0/2.d0)*csf_u2-(1.d0/2.d0)*csf_u1)
     !! RHS_v=v2-dt*((3.d0/2.d0)*conv2_v-(1.d0/2.d0)*conv1_v)+(dt/2.d0)*diffusion_v+dt*((3.d0/2.d0)*csf_v2-(1.d0/2.d0)*csf_v1)
     !! RHS_w=w2-dt*((3.d0/2.d0)*conv2_w-(1.d0/2.d0)*conv1_w)+(dt/2.d0)*diffusion_w+dt*((3.d0/2.d0)*csf_w2-(1.d0/2.d0)*csf_w1)

     !! AB3:
     !$omp parallel workshare
     RHS_u=u2-dt*((23.d0/12.d0)*conv2_u-(4.d0/3.d0)*conv1_u+(5.d0/12.d0)*conv0_u)+(dt/2.d0)*diffusion_u + &
          dt*((23.d0/12.d0)*csf_u2-(4.d0/3.d0)*csf_u1+(5.d0/12.d0)*csf_u0)

     RHS_v=v2-dt*((23.d0/12.d0)*conv2_v-(4.d0/3.d0)*conv1_v+(5.d0/12.d0)*conv0_v)+(dt/2.d0)*diffusion_v + &
          dt*((23.d0/12.d0)*csf_v2-(4.d0/3.d0)*csf_v1+(5.d0/12.d0)*csf_v0)

     RHS_w=w2-dt*((23.d0/12.d0)*conv2_w-(4.d0/3.d0)*conv1_w+(5.d0/12.d0)*conv0_w)+(dt/2.d0)*diffusion_w + &
          dt*((23.d0/12.d0)*csf_w2-(4.d0/3.d0)*csf_w1+(5.d0/12.d0)*csf_w0)

     ! Solving helmholtz operator at level n+1 for u
     ! Initial guess for u3:

     u3=0.d0
     !$omp end parallel workshare
     do iteration_sor=1,(max_iteration_mom_u/2)
        u3_old=u3

        call exchange2d(u3_old,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
        call do_jacobi_u(u3,u3_old,RHS_u,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,ex_max)

        if(sx==1)then
           u3(0,sy:ey,:)=u3(1,sy:ey,:)
        end if

        if(ex==ex_max)then
           u3(ex_max,sy:ey,:)=u3(ex_max-1,sy:ey,:)+                  &
                (dt/dx)*(pres(ex_max+1,sy:ey,1:ez_p-1)+pres(ex_max-1,sy:ey,1:ez_p-1)  &
                -2.d0*pres(ex_max,sy:ey,1:ez_p-1))

           ! last u3 point is completely fictitious and exists only for
           ! accounting reasons.
           u3(ex_max+1,sy:ey,:)=0.d0
        end if

        call exchange2d(u3,stride_uv_xz,stride_uv_yz,neighbours,  &
             ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
        call do_sor_u(u3,RHS_u,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,ex_max,0,iteration_sor)
        call exchange2d(u3,stride_uv_xz,stride_uv_yz,neighbours,  &
             ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
        call do_sor_u(u3,RHS_u,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,ex_max,1,iteration_sor)

        ! Special conditions on intermediate velocity u*: du/dx=0 at the inlet
        ! and the outlet.
        if(sx==1)then
           u3(0,sy:ey,:)=u3(1,sy:ey,:)
        end if

        if(ex==ex_max)then
           u3(ex_max,sy:ey,:)=u3(ex_max-1,sy:ey,:)+                  &
                (dt/dx)*(pres(ex_max+1,sy:ey,1:ez_p-1)+pres(ex_max-1,sy:ey,1:ez_p-1)  &
                -2.d0*pres(ex_max,sy:ey,1:ez_p-1))

           ! last u3 point is completely fictitious and exists only for
           ! accounting reasons.
           u3(ex_max+1,sy:ey,:)=0.d0
        end if

     end do

     call get_difference(u3,u3_old,ex,ey,ez_uv,sx,sy,sz_uv,err2)
     CALL MPI_ALLREDUCE(err2, err1,  1, MPI_DOUBLE_PRECISION,   MPI_SUM,PETSC_COMM_WORLD, ierr)
     if(ierr /= 0 ) call abort(ierr, 'Failed to call MPI_ALLREDUCE on err2 and err1 after get_difference',my_id)
     if (is_master(my_id)) then
        Write(*,*)'Iteration=',iteration, 'u---residual is ', err1
     end if

     ! Need the next step if unlet conditions on v and w are to involve pressure differences at level n:
     call exchange2d(pres,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

     !$omp parallel workshare
     v3=0.d0
     !$omp end parallel workshare
     do iteration_sor=1,(max_iteration_mom_v/2)
        v3_old=v3
        call exchange2d(v3_old,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
        call do_jacobi_v(v3,v3_old,RHS_v,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,ex_max)

        call exchange2d(v3,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
        call do_sor_v(v3,RHS_v,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,0,iteration_sor)
        call exchange2d(v3,stride_uv_xz,stride_uv_yz,neighbours,ex,ey,ez_uv,sx,sy,sz_uv,comm2d_quasiperiodic)
        call do_sor_v(v3,RHS_v,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,1,iteration_sor)

        if(sx==1)then
           do k=0,maxn-2
              do j=sy,ey
                 v3(0,j,k)=v3(1,j,k)
              end do
           end do
        end if

        if(ex==ex_max)then
           v3(ex_max+1,sy:ey,:)=v3(ex_max,sy:ey,:)
        end if

     end do

     call get_difference(v3,v3_old,ex,ey,ez_uv,sx,sy,sz_uv,err2)
     CALL MPI_ALLREDUCE(err2, err1,  1, MPI_DOUBLE_PRECISION,   MPI_SUM,PETSC_COMM_WORLD, ierr)
     if(ierr /= 0 ) call abort(ierr, 'Failed to call MPI_ALLREDUCE on err2 and err1 after get_difference',my_id)
     if (is_master(my_id)) then
        Write(*,*)'Iteration=',iteration, 'v---residual is ', err1
     end if

     !$omp parallel workshare
     w3=0.d0
     !$omp end parallel workshare
     do iteration_sor=1,(max_iteration_mom_w/2)
        w3_old=w3
        call exchange2d(w3_old,stride_w_xz,stride_w_yz,neighbours,ex,ey,ez_w,sx,sy,sz_w,comm2d_quasiperiodic)
        call do_jacobi_w(w3,w3_old,RHS_w,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,ex_max)

        call exchange2d(w3,stride_w_xz,stride_w_yz,neighbours,ex,ey,ez_w,sx,sy,sz_w,comm2d_quasiperiodic)
        call do_sor_w(w3,RHS_w,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,0,iteration_sor)
        call exchange2d(w3,stride_w_xz,stride_w_yz,neighbours,ex,ey,ez_w,sx,sy,sz_w,comm2d_quasiperiodic)
        call do_sor_w(w3,RHS_w,viscosity2_aug,dx,dy,dz,dt,ex,ey,sx,sy,maxn,1,iteration_sor)

        !$omp parallel workshare
        w3(sx:ex,sy:ey,0)=0.d0
        w3(sx:ex,sy:ey,maxn-1)=0.d0
        !$omp end parallel workshare

        if(sx==1)then
           do k=1,maxn-2
              do j=sy,ey
                 w3(0,j,k)=w3(1,j,k)
              end do
           end do

           w3(0,sy:ey,0)=0.d0
           w3(0,sy:ey,maxn-1)=0.d0

        end if

        if(ex==ex_max)then
           w3(ex_max+1,sy:ey,:)=w3(ex_max,sy:ey,:)
        end if
     end do

     call get_difference(w3,w3_old,ex,ey,ez_w,sx,sy,sz_w,err2)
     CALL MPI_ALLREDUCE(err2, err1,  1, MPI_DOUBLE_PRECISION,   MPI_SUM,PETSC_COMM_WORLD, ierr)
     if(ierr /= 0 ) call abort(ierr, 'Failed to call MPI_ALLREDUCE on err2, err1 after get_difference',my_id)
     if (is_master(my_id)) then
        Write(*,*)'Iteration=',iteration, 'w---residual is ', err1
     end if

     call exchange2d(u3,   stride_uv_xz, stride_uv_yz, neighbours,ex,ey,ez_uv,sx,sy,sz_uv, comm2d_quasiperiodic)
     call exchange2d(v3,   stride_uv_xz, stride_uv_yz, neighbours,ex,ey,ez_uv,sx,sy,sz_uv, comm2d_quasiperiodic)
     call exchange2d(w3,    stride_w_xz, stride_w_yz,  neighbours,ex,ey, ez_w,sx,sy, sz_w, comm2d_quasiperiodic)

     call get_source_pres(RHS_p,u3,v3,w3,ex,ey,sx,sy,maxn,dx,dy,dz,dt)
     call exchange2d(RHS_p,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)

     !$omp parallel workshare
     pres=pres_old
     !$omp end parallel workshare

     call KSPSolve(ksp, PETSC_NULL_OBJECT, PETSC_NULL_OBJECT, ierr)
     if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPSolve in iteration loop ',my_id)

     if (is_master(my_id)) then
        call output_converged_reason(ksp)
        call KSPGetResidualNorm(ksp, rnorm, ierr)
        if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPGetResidualNorm when checking convergence',my_id)
        write(*, *) 'approximate, preconditioned, residual norm', rnorm
        call KSPGetIterationNumber(ksp, its, ierr)
        if(ierr /= 0 ) call abort(ierr, 'Failed to call KSPGetIterationNumber when checking convergence',my_id)
        write(*, *) 'iterations =', its
     end if

     call copy_pressure_out_of_petsc(ksp, ierr)
     if(ierr /= 0 ) call abort(ierr, 'Failed to call copy_pressure_out_of_petsc',my_id)

     ! This exchange is not needed because 'copy_pressure_out_of_petsc' takes care of the ghost points.
     !call exchange2d(pres,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez_p,sx,sy,sz_p,comm2d_quasiperiodic)
     call get_uvw_pres(u3,v3,w3,pres,ex,ey,sx,sy,maxn,ex_max,dx,dy,dz,dt,u_inlet)


     ! ***************************************************************************
     io_t0 = mpi_wtime()
     if(backuphdf5) then
        ! Periodic output of phi file          
        if(mod(iteration,phi_dat_frequency).eq.0)then
           filename = 'phi'
           !           if(is_master(my_id)) write(*,*) 'Outputting to ',filename

           call output_phi_hdf5(filename,sx,ex,sy,ey,sz_p,ez_p,&
                phi2(sx:ex,sy:ey,sz_p:ez_p),num_procs_x,num_procs_y,&
                num_procs_z,my_id, iteration)
           ! call output_3D_hdf5(filename,sx,ex,sy,ey,sz_uv,ez_uv,&
           !      phi2(sx:ex,sy:ey,sz_p:ez_p),num_procs_x,num_procs_y,&
           !      num_procs_z,my_id,iteration)
        end if

        ! **************************************************************************                  
        ! Periodic output of uvw file          
        if(mod(iteration,uvw_dat_frequency).eq.0)then
           call output_uvw_hdf5(sx,ex,sy,ey,sz_uv,ez_uv,sz_w,ez_w,sz_p,ez_p,&
                output_u(sx:ex,sy:ey,sz_uv:ez_uv),v2(sx:ex,sy:ey,sz_uv:ez_uv),&
                output_w(sx:ex,sy:ey,sz_w:ez_w),output_p(sx:ex,sy:ey,sz_p:ez_p),&
                num_procs_x,num_procs_y, num_procs_z,my_id,iteration)
        end if

        ! **************************************************************************
        io_t1 = mpi_wtime()
        t_netcdf_reg = io_t1-io_t0

        if(mod(iteration,backup_frequency).eq.0)then
           backup_counter=backup_counter+1

           if(mod(backup_counter,2).eq.0)then
              write(strg,*)'0'
              backup=0
           else
              write(strg,*)'1'
              backup=1
           end if
           !if(is_master(my_id)) write(*,*) 'Outputting to ',filename
           call output_backup_hdf5(backup,dx,dy,dz,sx,ex,sy,ey,sz_p,ez_p,sz_uv,&
                ez_uv,sz_w,ez_w,u2(sx:ex,sy:ey,sz_uv:ez_uv),v2(sx:ex,sy:ey,sz_uv:ez_uv),&
                w2(sx:ex,sy:ey,sz_w:ez_w),conv0_u(sx:ex,sy:ey,sz_uv:ez_uv),&
                conv1_u(sx:ex,sy:ey,sz_uv:ez_uv),conv0_v(sx:ex,sy:ey,sz_uv:ez_uv),&
                conv1_v(sx:ex,sy:ey,sz_uv:ez_uv),conv0_w(sx:ex,sy:ey,sz_w:ez_w),&
                conv1_w(sx:ex,sy:ey,sz_w:ez_w),pres_old(sx:ex,sy:ey,sz_p:ez_p),&
                pres(sx:ex,sy:ey,sz_p:ez_p),phi2(sx:ex,sy:ey,sz_p:ez_p),&
                conv0_phi(sx:ex,sy:ey,sz_p:ez_p),conv1_phi(sx:ex,sy:ey,sz_p:ez_p),&
                num_procs_x,num_procs_y, num_procs_z,my_id)
        end if
     end if

     ! ****************************************************************************
     io_t2 = mpi_wtime()
     t_netcdf_bak = io_t2-io_t1
     t_netcdf_tot = io_t2-io_t0

     if(backuptext) then

        ! Periodic output to files          
        if(mod(iteration,phi_dat_frequency).eq.0)then

           ! ***********************************************
           ! gather local pressures to master node

           output_p=pres

           If (is_master(my_id)) Then
              pres_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       pres_global(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)=output_p(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank when gathering local pressures', my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv when gathering local pressures on sy_f',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to callmpi_recv when gathering local pressures on ey_f',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,&
                            'Failed to call mpi_recv when gathering local pressures on  sx_f',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv when gathering local pressures on ex_f',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv when gathering local pressures on n_local_x',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv when gathering local pressures on n_local_y',my_id)

                       Call mpi_recv(tempr_p,n_local_x*n_local_y*n_local_z_p,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv when gathering local pressures on tempr_p',my_id)

                       do k=1,n_local_z_p
                          do j=1,n_local_y
                             do i=1,n_local_x
                                pres_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_p(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend when sending pressures for sy_1',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend when sending pressures for ey_1',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend when sending pressures for sx_1',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend when sending pressures for ex_1',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend when sending pressures for n_local_x',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend when sending pressures for n_local_y',my_id)

              Call mpi_ssend(output_p(sx:ex,sy:ey,sz_p:ez_p),n_local_x*n_local_y*n_local_z_p,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call when sending pressures for output_p',my_id)

           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then            
              pres_global(0,:,:)=pres_global(1,:,:)
              pres_global(maxl,:,:)=pres_global(maxl-1,:,:)

              pres_global(:,0,:)=pres_global(:,maxm-1,:)
              pres_global(:,maxm,:)=pres_global(:,1,:)
              pres_global(:,:,maxn)=pres_global(:,:,maxn-1)
              pres_global(:,:,0)=pres_global(:,:,1)
           end if


           ! ***********************************************
           ! gather local pressures (old pressure) to master node

           If (is_master(my_id)) Then
              pres_old_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       pres_old_global(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)=pres_old(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank in gather local pressures from comm2d',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv when gather local pressures sy_f',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv when gather local pressures on ey_f',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv when gather local pressures on sx_f',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv when gather local pressures on ex_f',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv when gather local pressures on n_local_x',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv when gather local pressures on n_local_y',my_id)

                       Call mpi_recv(tempr_p,n_local_x*n_local_y*n_local_z_p,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv when gather local pressures on tempr_p',my_id)

                       do k=1,n_local_z_p
                          do j=1,n_local_y
                             do i=1,n_local_x
                                pres_old_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_p(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during local pressure gather',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ey during local pressure gather',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sx during local pressure gather',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ex during local pressure gather',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_x during local pressure gather',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during local pressure gather',my_id)

              Call mpi_ssend(pres_old(sx:ex,sy:ey,sz_p:ez_p),n_local_x*n_local_y*n_local_z_p,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend pres_old during local pressure gather',my_id)

           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then            
              pres_old_global(0,:,:)=pres_old_global(1,:,:)
              pres_old_global(maxl,:,:)=pres_old_global(maxl-1,:,:)

              pres_old_global(:,0,:)=pres_old_global(:,maxm-1,:)
              pres_old_global(:,maxm,:)=pres_old_global(:,1,:)
              pres_old_global(:,:,maxn)=pres_old_global(:,:,maxn-1)
              pres_old_global(:,:,0)=pres_old_global(:,:,1)
           end if

           ! ***********************************************
           ! gather local phase field to master node

           If (is_master(my_id)) Then
              phi_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       phi_global(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)=phi2(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather of local phase field',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sy_f during gather of local phase field',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ey_f during gather of local phase field',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather of local phase f7ield ',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ex_f during gather of local phase field',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_x during gather of local phase field',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_y during gather of local phase field ',my_id)

                       Call mpi_recv(tempr_phi,n_local_x*n_local_y*n_local_z_p,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on tempr_phi during gather of local phase field',my_id)

                       do k=1,n_local_z_p
                          do j=1,n_local_y
                             do i=1,n_local_x
                                phi_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_phi(i,j,k)
                             end do
                          end do
                       end do
                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during gather of local phase field ',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call  mpi_ssend on ey during gather of local phase field',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call  mpi_ssend on sx during gather of local phase field',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call  mpi_ssend on ex during gather of local phase field',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call  mpi_ssend on n_local_x during gather of local phase field',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather of local phase field ',my_id)

              Call mpi_ssend(phi2(sx:ex,sy:ey,sz_p:ez_p),n_local_x*n_local_y*n_local_z_p,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on phi2 during gather of local phase field ',my_id)

           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then

              do k=1,maxn-1
                 do j=1,maxm-1
                    z_val=k*dz-0.5d0*dz
                    phi_global(0,j,k)=z_val-height
                 end do
              end do

              phi_global(maxl,:,:)=phi_global(maxl-1,:,:)

              phi_global(:,0,:)=phi_global(:,maxm-1,:)
              phi_global(:,maxm,:)=phi_global(:,1,:)
              phi_global(:,:,maxn)=phi_global(:,:,maxn-1)
              phi_global(:,:,0)=phi_global(:,:,1)
           end if

           ! ***********************************************
           ! gather local phase field convective derivatives to master node

           If (is_master(my_id)) Then
              conv0_phi_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       conv0_phi_global(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)=conv0_phi(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr,& 
                            'Failed to call mpi_cart_rank during gather local phase field convective derivatives',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sy_f during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ey_f during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather local phase field convective derivatives ',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ex_f during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer, sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_x during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer, sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_y during gather local phase field convective derivatives',my_id)

                       Call mpi_recv(tempr_phi,n_local_x*n_local_y*n_local_z_p,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on tempr_phi during gather local phase field convective derivatives',my_id)

                       do k=1,n_local_z_p
                          do j=1,n_local_y
                             do i=1,n_local_x
                                conv0_phi_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_phi(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during gather local phase field convective derivatives',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ey during gather local phase field convective derivatives ',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sx during gather local phase field convective derivatives ',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ex during gather local phase field convective derivatives ',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_x during gather local phase field convective derivatives ',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local phase field convective derivatives ',my_id)

              Call mpi_ssend(conv0_phi(sx:ex,sy:ey,sz_p:ez_p),n_local_x*n_local_y*n_local_z_p,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on conv0_phi during gather local phase field convective derivatives',my_id)

           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              conv0_phi_global(0,:,:)=conv0_phi_global(1,:,:)
              conv0_phi_global(maxl,:,:)=conv0_phi_global(maxl-1,:,:)

              conv0_phi_global(:,0,:)=conv0_phi_global(:,maxm-1,:)
              conv0_phi_global(:,maxm,:)=conv0_phi_global(:,1,:)
              conv0_phi_global(:,:,maxn)=conv0_phi_global(:,:,maxn-1)
              conv0_phi_global(:,:,0)=conv0_phi_global(:,:,1)
           end if

           ! ***********************************************
           ! gather local phase field convective derivatives to master node

           If (is_master(my_id)) Then
              conv1_phi_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       conv1_phi_global(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)=conv1_phi(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr,&
                            'Failed to call mpi_cart_rank during gather local phase field convective derivatives',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,&
                            'Failed to call mpi_recv on sy_f, during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,& 
                            'Failed to call mpi_recv on ey_f during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,& 
                            'Failed to call mpi_recv on sx_f during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,& 
                            'Failed to call mpi_recv on ex_f during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,& 
                            'Failed to call mpi_recv on n_local_x during gather local phase field convective derivatives',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,& 
                            'Failed to call mpi_recv on n_local_y during gather local phase field convective derivatives',my_id)

                       Call mpi_recv(tempr_phi,n_local_x*n_local_y*n_local_z_p,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,& 
                            'Failed to call mpi_recv on tempr_phi during gather local phase field convective derivatives',my_id)

                       do k=1,n_local_z_p
                          do j=1,n_local_y
                             do i=1,n_local_x
                                conv1_phi_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_phi(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,& 
                   'Failed to call mpi_ssend on sy during gather local phase field convective derivatives ',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,& 
                   'Failed to call mpi_ssend on ey during gather local phase field convective derivatives ',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,&
                   'Failed to call mpi_ssend on sx during gather local phase field convective derivatives',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,&
                   'Failed to call mpi_ssend on ex during gather local phase field convective derivatives ',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,& 
                   'Failed to call mpi_ssend on n_local_x during gather local phase field convective derivatives ',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on  n_local_yduring gather local phase field convective derivatives',my_id)

              Call mpi_ssend(conv1_phi(sx:ex,sy:ey,sz_p:ez_p),n_local_x*n_local_y*n_local_z_p,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call  mpi_ssend on conv1_phi during gather local phase field convective derivatives',my_id)

           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              conv1_phi_global(0,:,:)=conv1_phi_global(1,:,:)
              conv1_phi_global(maxl,:,:)=conv1_phi_global(maxl-1,:,:)

              conv1_phi_global(:,0,:)=conv1_phi_global(:,maxm-1,:)
              conv1_phi_global(:,maxm,:)=conv1_phi_global(:,1,:)
              conv1_phi_global(:,:,maxn)=conv1_phi_global(:,:,maxn-1)
              conv1_phi_global(:,:,0)=conv1_phi_global(:,:,1)
           end if

           ! ***********************************************
           ! gather local viscosities to master node

           If (is_master(my_id)) Then
              visc_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       visc_global(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)=viscosity2(sx_o:ex_o,sy_o:ey_o,sz_p:ez_p)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local viscosities',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sy_f during gather local viscosities',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ey_f during gather local viscosities ',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather local viscosities ',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ex_f during gather local viscosities ',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_x during gather local viscosities ',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_y during gather local viscosities ',my_id)
                       Call mpi_recv(tempr_visc,n_local_x*n_local_y*n_local_z_p,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on tempr_visc during gather local viscosities ',my_id)

                       do k=1,n_local_z_p
                          do j=1,n_local_y
                             do i=1,n_local_x
                                visc_global(sx_f-1+i,sy_f-1+j,sz_p-1+k) =tempr_visc(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during gather local viscosities',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ey during gather local viscosities ',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sx during gather local viscosities ',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call  mpi_ssend on ex during gather local viscosities',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call  mpi_ssend on n_local_x during gather local viscosities',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local viscosities ',my_id)
              Call mpi_ssend(viscosity2(sx:ex,sy:ey,sz_p:ez_p),n_local_x*n_local_y*n_local_z_p,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on viscosity2 during gather local viscosities ',my_id)
           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              visc_global(0,:,:)=visc_global(1,:,:)
              visc_global(maxl,:,:)=visc_global(maxl-1,:,:)

              visc_global(:,0,:)=visc_global(:,maxm-1,:)
              visc_global(:,maxm,:)=visc_global(:,1,:)
              visc_global(:,:,maxn)=visc_global(:,:,maxn-1)
              visc_global(:,:,0)=visc_global(:,:,1)
           end if

           ! ***********************************************
           ! gather local u-velocities to master node

           output_u=u2

           If (is_master(my_id)) Then
              u_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       u_global(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)=output_u(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local u-velocities',my_id)

                       Call mpi_recv(sy_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on sy_f during gather local u-velocities',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ey_f during gather local u-velocities ',my_id)
                       Call mpi_recv(sx_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call  mpi_recv on sx_f during gather local u-velocities',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call  mpi_recv on ex_f during gather local u-velocities',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call  mpi_recv on n_local_x during gather local u-velocities',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on n_local_y during gather local u-velocities ',my_id)
                       Call mpi_recv(tempr_u,n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on tempr_u during gather local u-velocities ',my_id)

                       do k=1,n_local_z_uv
                          do j=1,n_local_y
                             do i=1,n_local_x
                                u_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) =tempr_u(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else

              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during gather local u-velocities',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on ey during gather local u-velocities',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on sx during gather local u-velocities',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on ex during gather local u-velocities',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on n_local_x during gather local u-velocities',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on n_local_y during gather local u-velocities',my_id)

              Call mpi_ssend(output_u(sx:ex,sy:ey,sz_uv:ez_uv),n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on output_u during gather local u-velocities',my_id)
           End If

           ! enforcing  boundary conditions on ghost cells

           if (is_master(my_id)) then
              do k=0,maxn-2
                 u_global(0,:,k)=   u_inlet(k)
              end do
              u_global(maxl-1,:,:)=u_global(maxl-2,:,:)
              u_global(maxl,:,:)=u_global(maxl-2,:,:)

              u_global(:,0,:)=u_global(:,maxm-1,:)
              u_global(:,maxm,:)=u_global(:,1,:)
           end if

           ! ***********************************************
           ! gather local u-velocity convective derivatives to master node

           If (is_master(my_id)) Then
              conv0_u_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       conv0_u_global(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)=conv0_u(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local u-velocity convective derivatives',my_id)

                       Call mpi_recv(sy_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sy_f during gather local u-velocity convective derivatives ',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ey_f during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(sx_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ex_f during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,&
                            'Failed to call mpi_recv on n_local_x during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_y during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(tempr_u,n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on tempr_u during gather local u-velocity convective derivatives',my_id)

                       do k=1,n_local_z_uv
                          do j=1,n_local_y
                             do i=1,n_local_x
                                conv0_u_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) =tempr_u(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else

              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ey during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,&
                   'Failed to call mpi_ssend on sx during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ex during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_x during gather local u-velocity convective derivatives ',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(conv0_u(sx:ex,sy:ey,sz_uv:ez_uv),n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on conv0_u during gather local u-velocity convective derivatives',my_id)
           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              conv0_u_global(0,:,:)=   conv0_u_global(1,:,:)
              conv0_u_global(maxl,:,:)=conv0_u_global(maxl-1,:,:)

              conv0_u_global(:,0,:)=conv0_u_global(:,maxm-1,:)
              conv0_u_global(:,maxm,:)=conv0_u_global(:,1,:)
           end if

           ! ***********************************************
           ! gather local u-velocity convective derivatives to master node

           If (is_master(my_id)) Then
              conv1_u_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       conv1_u_global(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)=conv1_u(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local u-velocity convective derivatives ',my_id)

                       Call mpi_recv(sy_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sy_f during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,&
                            'Failed to call mpi_recv on ey_f during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(sx_f,1,  mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ex_f during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_x during gather local u-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_y during gather local u-velocity convective derivatives',my_id)

                       Call mpi_recv(tempr_u,n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on tempr_u during gather local u-velocity convective derivatives',my_id)

                       do k=1,n_local_z_uv
                          do j=1,n_local_y
                             do i=1,n_local_x
                                conv1_u_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) =tempr_u(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else

              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,&
                   'Failed to call mpi_ssend on ey during gather local u-velocity convective derivatives ',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sx during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ex during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_x during gather local u-velocity convective derivatives',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local u-velocity convective derivatives',my_id)

              Call mpi_ssend(conv1_u(sx:ex,sy:ey,sz_uv:ez_uv),n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on conv1_u during gather local u-velocity convective derivatives ',my_id)
           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              conv1_u_global(0,:,:)=   conv1_u_global(1,:,:)
              conv1_u_global(maxl,:,:)=conv1_u_global(maxl-1,:,:)

              conv1_u_global(:,0,:)=conv1_u_global(:,maxm-1,:)
              conv1_u_global(:,maxm,:)=conv1_u_global(:,1,:)
           end if

           ! ***********************************************
           ! gather local v-velocities to master node

           If (is_master(my_id)) Then
              v_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       v_global(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)=v2(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local v-velocities ',my_id)
                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sy_f during gather local v-velocities',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ey_f during gather local v-velocities ',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather local v-velocities',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ex_f during gather local v-velocities ',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_x during gather local v-velocities',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_y during gather local v-velocities',my_id)

                       Call mpi_recv(tempr_v,n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on tempr_v during gather local v-velocities',my_id)

                       do k=1,n_local_z_uv
                          do j=1,n_local_y
                             do i=1,n_local_x
                                v_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) =tempr_v(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during gather local v-velocities',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,& 
                   'Failed to call mpi_ssend on ey during gather local v-velocities',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sx during gather local v-velocities',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on ex during gather local v-velocities',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_x during gather local v-velocities',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local v-velocities',my_id)

              Call mpi_ssend(v2(sx:ex,sy:ey,sz_uv:ez_uv),n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_ssend on v2 during gather local v-velocities',my_id)
           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              v_global(0,:,:)=   0.d0
              v_global(maxl,:,:)=v_global(maxl-1,:,:)
              v_global(:,0,:)=v_global(:,maxm-1,:)
              v_global(:,maxm,:)=v_global(:,1,:)
           end if

           ! ***********************************************
           ! gather local v-velocity convective derivatives to master node

           If (is_master(my_id)) Then
              conv0_v_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       conv0_v_global(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)=conv0_v(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local v-velocity convective derivatives',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on sy_f during gather local v-velocity convective derivatives ',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on ey_f during gather local v-velocity convective derivatives ',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on sx_f during gather local v-velocity convective derivatives',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on ex_f during gather local v-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_x during gather local v-velocity convective derivatives ',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on n_local_y during gather local v-velocity convective derivatives',my_id)
                       Call mpi_recv(tempr_v,n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on tempr_v during gather local v-velocity convective derivatives',my_id)

                       do k=1,n_local_z_uv
                          do j=1,n_local_y
                             do i=1,n_local_x
                                conv0_v_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) =tempr_v(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on sy during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on ey during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on sx during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on ex during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on n_local_x during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local v-velocity convective derivatives',my_id)

              Call mpi_ssend(conv0_v(sx:ex,sy:ey,sz_uv:ez_uv),n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on conv0_v during gather local v-velocity convective derivatives',my_id)
           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              conv0_v_global(0,:,:)=   conv0_v_global(1,:,:)
              conv0_v_global(maxl,:,:)=conv0_v_global(maxl-1,:,:)

              conv0_v_global(:,0,:)=conv0_v_global(:,maxm-1,:)
              conv0_v_global(:,maxm,:)=conv0_v_global(:,1,:)
           end if

           ! ***********************************************
           ! gather local v-velocity convective derivatives to master node

           If (is_master(my_id)) Then
              conv1_v_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       conv1_v_global(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)=conv1_v(sx_o:ex_o,sy_o:ey_o,sz_uv:ez_uv)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local v-velocity convective derivatives',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sy_f during gather local v-velocity convective derivatives ',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ey_f during gather local v-velocity convective derivatives',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather local v-velocity convective derivatives',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ex_f during gather local v-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_x during gather local v-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on n_local_y during gather local v-velocity convective derivatives',my_id)

                       Call mpi_recv(tempr_v,n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on tempr_v during gather local v-velocity convective derivatives',my_id)

                       do k=1,n_local_z_uv
                          do j=1,n_local_y
                             do i=1,n_local_x
                                conv1_v_global(sx_f-1+i,sy_f-1+j,sz_uv-1+k) =tempr_v(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on sy during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to call mpi_ssend on ey during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to call mpi_ssend on sx during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to call mpi_ssend on ex during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to call mpi_ssend on n_local_x during gather local v-velocity convective derivatives',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local v-velocity convective derivatives',my_id)

              Call mpi_ssend(conv1_v(sx:ex,sy:ey,sz_uv:ez_uv),n_local_x*n_local_y*n_local_z_uv,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,&
                   'Failed to call mpi_ssend on conv1_v during gather local v-velocity convective derivatives',my_id)
           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              conv1_v_global(0,:,:)=   conv1_v_global(1,:,:)
              conv1_v_global(maxl,:,:)=conv1_v_global(maxl-1,:,:)

              conv1_v_global(:,0,:)=conv1_v_global(:,maxm-1,:)
              conv1_v_global(:,maxm,:)=conv1_v_global(:,1,:)
           end if

           ! ***********************************************
           ! gather local w-velocities to master node

           output_w=w2

           If (is_master(my_id)) Then
              w_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       w_global(sx_o:ex_o,sy_o:ey_o,sz_w:ez_w)=output_w(sx_o:ex_o,sy_o:ey_o,sz_w:ez_w)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            ' Failed to call mpi_cart_rank during gather local w-velocities',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, & 
                            'Failed to call mpi_recv on sy_f during gather local w-velocities',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on ey_f during gather local w-velocities',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on sx_f during gather local w-velocities',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, & 
                            'Failed to call mpi_recv on ex_f during gather local w-velocities',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, & 
                            'Failed to call mpi_recv on n_local_x during gather local w-velocities',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on n_local_y during gather local w-velocities ',my_id)

                       Call mpi_recv(tempr_w,n_local_x*n_local_y*n_local_z_w,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on tempr_w during gather local w-velocities',my_id)

                       do k=1,n_local_z_w
                          do j=1,n_local_y
                             do i=1,n_local_x
                                w_global(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_w(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to call mpi_ssend on sy during gather local w-velocities',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to callmpi_ssend on ey during gather local w-velocities ',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to call mpi_ssend on sx during gather local w-velocities',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to call mpi_ssend on ex during gather local w-velocities',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, & 
                   'Failed to call mpi_ssend on n_local_x during gather local w-velocities',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local w-velocities',my_id)

              Call mpi_ssend(output_w(sx:ex,sy:ey,sz_w:ez_w),n_local_x*n_local_y*n_local_z_w,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on output_w during gather local w-velocities',my_id)
           End If

           ! enforcing conditions on ghost cells

           if (is_master(my_id)) then
              w_global(:,:,0)=0.d0
              w_global(:,:,maxn-1)=0.d0

              w_global(0,:,:)=   0.d0
              w_global(maxl,:,:)=0.d0

              w_global(:,0,:)=w_global(:,maxm-1,:)
              w_global(:,maxm,:)=w_global(:,1,:)
           end if

           ! ***********************************************
           ! gather local w-velocity convective derivatives to master node

           If (is_master(my_id)) Then
              conv0_w_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       conv0_w_global(sx_o:ex_o,sy_o:ey_o,sz_w:ez_w)=conv0_w(sx_o:ex_o,sy_o:ey_o,sz_w:ez_w)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local w-velocity convective derivatives',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sy_f during gather local w-velocity convective derivatives',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to callmpi_recv on ey_f during gather local w-velocity convective derivatives ',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather local w-velocity convective derivatives',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on ex_f during gather local w-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on n_local_x during gather local w-velocity convective derivatives',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on n_local_y during gather local w-velocity convective derivatives',my_id)
                       Call mpi_recv(tempr_w,n_local_x*n_local_y*n_local_z_w,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on tempr_w during gather local w-velocity convective derivatives',my_id)

                       do k=1,n_local_z_w
                          do j=1,n_local_y
                             do i=1,n_local_x
                                conv0_w_global(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_w(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on sy during gather local w-velocity convective derivatives',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on ey during gather local w-velocity convective derivatives',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on sx during gather local w-velocity convective derivatives',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to callmpi_ssend on ex during gather local w-velocity convective derivatives ',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on n_local_x during gather local w-velocity convective derivatives',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr, &
                   'Failed to call mpi_ssend on n_local_y during gather local w-velocity convective derivatives',my_id)

              Call mpi_ssend(conv0_w(sx:ex,sy:ey,sz_w:ez_w),n_local_x*n_local_y*n_local_z_w,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on conv0_w during gather local w-velocity convective derivatives',my_id)
           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              conv0_w_global(0,:,:)=   conv0_w_global(1,:,:)
              conv0_w_global(maxl,:,:)=conv0_w_global(maxl-1,:,:)

              conv0_w_global(:,0,:)=conv0_w_global(:,maxm-1,:)
              conv0_w_global(:,maxm,:)=conv0_w_global(:,1,:)
              conv0_w_global(:,:,0)=0.d0
              conv0_w_global(:,:,maxn-1)=0.d0
           end if

           ! ***********************************************
           ! gather local w-velocity (convective derivative) to master node

           If (is_master(my_id)) Then
              conv1_w_global=0.d0
              pz=0
              do py=0,dims(2)-1
                 do px=0,dims(1)-1

                    if((pz==0).and.(py==0).and.(px==0))then
                       conv1_w_global(sx_o:ex_o,sy_o:ey_o,sz_w:ez_w)=conv1_w(sx_o:ex_o,sy_o:ey_o,sz_w:ez_w)    
                    else

                       coords(1)=px
                       coords(2)=py
                       coords(3)=pz

                       !! Find the sender's rank from the senders virtual Cartesian coords.
                       Call mpi_cart_rank(comm2d_quasiperiodic,coords,sender_id,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_cart_rank during gather local w-velocity (convective derivative)',my_id)

                       Call mpi_recv(sy_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on sy_f during gather local w-velocity (convective derivative)',my_id)
                       Call mpi_recv(ey_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on ey_f during gather local w-velocity (convective derivative)',my_id)
                       Call mpi_recv(sx_f,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr, &
                            'Failed to call mpi_recv on sx_f during gather local w-velocity (convective derivative)',my_id)
                       Call mpi_recv(ex_f,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on ex_f during gather local w-velocity (convective derivative)',my_id)
                       Call mpi_recv(n_local_x,1,mpi_integer,  sender_id,stag,comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on n_local_x during gather local w-velocity (convective derivative)',my_id)
                       Call mpi_recv(n_local_y,1,  mpi_integer,  sender_id,etag,comm2d_quasiperiodic,status,ierr)

                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on n_local_y during gather local w-velocity (convective derivative)',my_id)
                       Call mpi_recv(tempr_w,n_local_x*n_local_y*n_local_z_w,mpi_double_precision,sender_id,0, &
                            comm2d_quasiperiodic,status,ierr)
                       if(ierr /= 0 ) call abort(ierr,  &
                            'Failed to call mpi_recv on tempr_w during gather local w-velocity (convective derivative)',my_id)

                       do k=1,n_local_z_w
                          do j=1,n_local_y
                             do i=1,n_local_x
                                conv1_w_global(sx_f-1+i,sy_f-1+j,sz_w-1+k) =tempr_w(i,j,k)
                             end do
                          end do
                       end do

                    end if

                 end do
              end do
           Else
              Call mpi_ssend(sy,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on sy during gather local w-velocity (convective derivative)',my_id)
              Call mpi_ssend(ey,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on ey during gather local w-velocity (convective derivative)',my_id)
              Call mpi_ssend(sx,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on sx during gather local w-velocity (convective derivative) ',my_id)
              Call mpi_ssend(ex,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on ex during gather local w-velocity (convective derivative)',my_id)
              Call mpi_ssend(n_local_x,1,mpi_integer,master_id,stag,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on n_local_x during gather local w-velocity (convective derivative)',my_id)
              Call mpi_ssend(n_local_y,1,  mpi_integer,master_id,etag,comm2d_quasiperiodic,ierr)

              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on n_local_y during gather local w-velocity (convective derivative)',my_id)
              Call mpi_ssend(conv1_w(sx:ex,sy:ey,sz_w:ez_w),n_local_x*n_local_y*n_local_z_w,mpi_double_precision,   &
                   master_id,0,comm2d_quasiperiodic,ierr)
              if(ierr /= 0 ) call abort(ierr,  &
                   'Failed to call mpi_ssend on conv1_w during gather local w-velocity (convective derivative)',my_id)
           End If

           ! enforcing boundary conditions on ghost cells

           if (is_master(my_id)) then
              conv1_w_global(0,:,:)=   conv1_w_global(1,:,:)
              conv1_w_global(maxl,:,:)=conv1_w_global(maxl-1,:,:)

              conv1_w_global(:,0,:)=conv1_w_global(:,maxm-1,:)
              conv1_w_global(:,maxm,:)=conv1_w_global(:,1,:)
              conv1_w_global(:,:,0)=0.d0
              conv1_w_global(:,:,maxn-1)=0.d0
           end if

           ! ***********************************************
           ! Write to files
           io_t3 = mpi_wtime()
           t_text_comms = io_t3-io_t2

           if (is_master(my_id)) then

              write(filename,'(A,I9,A)')'phi_channel_',iteration,'.dat'
              call channel_output_phi(phi_global,maxl,maxm,maxn,dx,dy,dz,filename)

              if(mod(iteration,uvw_dat_frequency).eq.0)then 
                 write(filename,'(A,I9,A)')'uvw_channel_',iteration,'.dat'
                 call channel_output_uvw(u_global,v_global,w_global,pres_global, &
                      maxl,maxm,maxn,dx,dy,dz,filename)
              end if
              io_t4 = mpi_wtime()
              t_text_reg = io_t4-io_t3

              backup_counter=backup_counter+1
              if(mod(backup_counter,2).eq.0)then

                 write(filename,'(A,I9,A)')'fieldbackup0.dat'
                 call backup_channel(u_global,v_global,w_global, &
                      conv0_u_global,conv1_u_global, &
                      conv0_v_global,conv1_v_global, &
                      conv0_w_global,conv1_w_global, &
                      pres_old_global,pres_global,   &
                      phi_global,conv0_phi_global,conv1_phi_global,maxl,maxm,maxn,filename)
              else
                 write(filename,'(A,I9,A)')'fieldbackup1.dat'
                 call backup_channel(u_global,v_global,w_global, &
                      conv0_u_global,conv1_u_global, &
                      conv0_v_global,conv1_v_global, &
                      conv0_w_global,conv1_w_global, &
                      pres_old_global,pres_global,   &
                      phi_global,conv0_phi_global,conv1_phi_global,maxl,maxm,maxn,filename)
              end if
              io_t5 = mpi_wtime()
              t_text_bak = io_t5-io_t4
              t_text_tot = io_t5-io_t2

           end if

        end if
     end if !end if for text backup option

     ! *************************************************************************************************

     !$omp parallel workshare
     conv0_u=conv1_u
     conv0_v=conv1_v
     conv0_w=conv1_w

     conv1_u=conv2_u
     conv1_v=conv2_v
     conv1_w=conv2_w

     csf_u0=csf_u1
     csf_v0=csf_v1
     csf_w0=csf_w1

     csf_u1=csf_u2
     csf_v1=csf_v2
     csf_w1=csf_w2

     conv0_phi=conv1_phi
     conv1_phi=conv2_phi

     pres_old=pres
     viscosity1=viscosity2

     u2=u3
     v2=v3
     w2=w3

     phi2=phi3
     phi_old = phi2
     !$omp end parallel workshare

     ! Timers
     if(backuptext) then
        text_back_time = text_back_time+t_text_bak
        text_reg_time = text_reg_time+t_text_reg
        text_comms_time = text_comms_time+t_text_comms
        total_text_time = total_text_time+t_text_tot
     end if
     if(backuphdf5) then
        hdf5_back_time = hdf5_back_time+t_netcdf_bak
        hdf5_reg_time = hdf5_reg_time+t_netcdf_reg
        total_hdf5_time = total_hdf5_time+t_netcdf_tot
     end if

     ! *********************************************************************
     ! Step out of time iteration loop 
  end do

  ! ********************************************************************

  t2 = mpi_wtime()

  if(is_master(my_id))then
     write(*,*) 'finalising and deallocating'
  end if

  deallocate(petsc_y_ranges)
  deallocate(u_global, v_global, w_global)
  deallocate(pres_global)
  deallocate(phi_global)
  deallocate(visc_global)
  deallocate(conv0_u_global, conv1_u_global)
  deallocate(conv0_v_global, conv1_v_global)
  deallocate(conv0_w_global, conv1_w_global)
  deallocate(conv0_phi_global, conv1_phi_global)
  deallocate(pres_old_global)

  deallocate(u2,u3,u3_old,conv0_u,conv1_u,conv2_u,diffusion_u,RHS_u,tempr_u,output_u)
  deallocate(v2,v3,v3_old,conv0_v,conv1_v,conv2_v,diffusion_v,RHS_v,tempr_v)
  deallocate(w2,w3,w3_old,conv0_w,conv1_w,conv2_w,diffusion_w,RHS_w,tempr_w,output_w)
  deallocate(pres,pres_old,RHS_p,tempr_p,output_p)
  deallocate(u2_aug,v2_aug,w2_aug,viscosity1,viscosity2,viscosity2_aug,tempr_visc)
  deallocate(phi2,phi3,conv2_phi,conv1_phi,conv0_phi,ConvPrevPhi)
  deallocate(csf_u0,csf_u1,csf_u2,csf_v0,csf_v1,csf_v2,csf_w0,csf_w1,csf_w2)
  deallocate(fx_csf,fy_csf,fz_csf)
  deallocate(u_inlet)

  call mpi_comm_free (comm2d_quasiperiodic, ierr)
  if(ierr /= 0 ) call abort(ierr, 'Failed to call mpi_comm_free',my_id)

  call PetscFinalize(ierr)
  if(ierr /= 0 ) call abort(ierr, 'Failed to call PetscFinalize',my_id)

  t3 = mpi_wtime() 
  t_total = mpi_wtime()-t0
  t_init = t1-t0
  t_iter = t2-t1
  t_iter_rate = t_iter/n_timesteps
  t_fina = t3-t2

  !Print out timings
  if(is_master(my_id))then    
     if(backuptext) then
        write(*,*) 'Timing: text_backup: ',text_back_time
        write(*,*) 'Timing: text_regular_output: ',text_reg_time
        write(*,*) 'Timing: text_communications_time: ',text_comms_time
        write(*,*) 'Timing: total_text_io_time: ',total_text_time
     end if
     if(backuphdf5) then
        write(*,*) 'Timing: netcdf_backup: ',hdf5_back_time
        write(*,*) 'Timing: netcdf_regular_output: ',hdf5_reg_time
        write(*,*) 'Timing: total_netcdf_io_time: ',total_hdf5_time
     end if
     write(*,*) 'Timing: setup: ',t_init
     write(*,*) 'Timing: iteration_loop: ',t_iter
     write(*,*) 'Timing: per_iteration: ',t_iter_rate
     write(*,*) 'Timing: finalisation: ',t_fina
     write(*,*) 'Timing: total_execution: ',t_total
  end if

end program tpls_program
