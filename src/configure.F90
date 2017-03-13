!> Subroutines to get and validate TPLS configuration options.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune, 
!! Toni Collis, David Scott, Peter Spelt, Mike Jackson.
!! @version $Revision: 252 $
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt, The University of
!! Edinburgh, all rights reserved.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module tpls_configure

  use grid_utils
  use option_names
  use options_utils
  use tpls_io

  implicit none

  !> Configuration file not found error.
  integer :: config_file_not_found_err
  parameter(config_file_not_found_err=200)
  !> Configuration directory not found error.
  integer :: config_dir_not_found_err
  parameter(config_dir_not_found_err=201)

contains

  !> Get configuration options and configuration directory from
  !! command-line.
  !! - Get configuration filename from command-line <tt>-f</tt> flag. If
  !!   none then use default_filename.
  !! - Get configuration directory from command-line <tt>-d</tt> flag. If
  !!   none then use default_directory.
  !! - Load configuration options from configuration file.
  !! - If directory does not exist or file does not exist then an
  !!   error code and message are set.
  subroutine get_configuration(default_filename,&
       default_config_dir,config_dir,options,num_options,ierr,message)

    character(len=*),  intent(in)  :: default_filename
    character(len=*),  intent(in)  :: default_config_dir
    character(len=*),  intent(out) :: config_dir
    character(len=80),dimension(:,:),allocatable,intent(out) :: options
    integer,           intent(out) :: num_options
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message

    integer :: num_cli_options
    logical :: file_exists
    character(len=80), dimension(:,:), allocatable :: cli_options
    character(len=256) :: filename

    call get_command_line_options(cli_options,num_cli_options)
    if (has_option(cli_options,num_cli_options,opt_file)) then
       call get_typed_option(cli_options,num_cli_options,opt_file,filename,ierr)
    else
       filename=default_filename
    end if
    if (has_option(cli_options,num_cli_options,opt_config_dir)) then
       call get_typed_option(cli_options,num_cli_options,opt_config_dir,config_dir,ierr)
    else
       config_dir=default_config_dir
    end if
    deallocate(cli_options)

    inquire(file=config_dir,exist=file_exists)
    if (.not. file_exists) then
       ierr = config_dir_not_found_err
       write(message,'(A,A)') 'Error: Configuration directory not found: ',trim(config_dir)
       return
    end if

    write(*,'(A,A)') 'Loading options from options file: ',trim(filename)
    call load_options(filename,options,num_options,ierr)
    if (ierr==option_file_not_found_err) then
       write(message,'(A,A)') 'Error: Options file not found: ',trim(filename)
       return
    else if (ierr/=0) then
       write(message,'(A,A)') 'Error: Loading options file: ',trim(filename)
       return
    end if
  end subroutine get_configuration


  !> Load a data file from a configuration directory.
  !! If the file does not exist or the data dimensions do not match
  !! then an error code and message are set. 
  subroutine load_configuration_data(config_dir,filename,maxl,maxm,maxn,data,ierr,message)

    character(*),      intent(in)  :: config_dir
    character(*),      intent(in)  :: filename ! Local file name in config_dir.
    integer,           intent(in)  :: maxl
    integer,           intent(in)  :: maxm
    integer,           intent(in)  :: maxn
    double precision,allocatable,dimension(:,:,:),intent(inout) :: data
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message

    character(len=60)  :: filepath
    logical            :: file_exists

    call get_filepath(config_dir,filename,filepath)
    inquire(file=filepath,exist=file_exists)
    if (.not. file_exists) then
       ierr = config_file_not_found_err
       write(message,'(A,A)') 'Error: Configuration file not found: ',trim(filepath)
       return
    end if
    call load_data_into_array(filepath,data,0,maxl,0,maxm,0,maxn,ierr)
    if (ierr/=0) then
       deallocate(data)
       write(message,'(A,A,A,I6,I6,I6)') 'Error: ',trim(filepath),&
            ' data dimensions do not match ',maxl,maxm,maxn
       return
    end if
  end subroutine load_configuration_data


  ! Subroutines used to parse both options file used as input to
  ! initial conditions program and also to parse initial conditions
  ! options used as input to TPLS itself.


  !> Get domain configuration.
  !! There are assumed to be options:
  !! - maxl
  !! - maxm
  !! - maxn
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_grid_configuration(options,num_options,&
       maxl,maxm,maxn,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,           intent(in)  :: num_options
    integer,           intent(out) :: maxl
    integer,           intent(out) :: maxm
    integer,           intent(out) :: maxn
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options,num_options,opt_maxl,maxl,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_maxl)
       return
    else if (maxl < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(opt_maxl,maxl,'Value must be >= 1')
       return
    end if
    call get_typed_option(options,num_options,opt_maxm,maxm,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_maxm)
       return
    else if (maxm < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(opt_maxm,maxm,'Value must be >= 1')
       return
    end if
    call get_typed_option(options,num_options,opt_maxn,maxn,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_maxn)
       return
    else if (maxn < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(opt_maxn,maxn,'Value must be >= 1')
       return
    end if
    if (show) then
       write(*,'(A)') 'Grid dimensions: '
       write(*,'(A,I6)') ' L: ', maxl
       write(*,'(A,I6)') ' M: ', maxm
       write(*,'(A,I6)') ' N: ', maxn
    end if
  end subroutine get_grid_configuration


  !> Get fluid configuration.
  !! There are assumed to be options:
  !! - re
  !! - dpdl
  !! - mu_minus
  !! - mu_plus
  !! - height
  !! - scap
  !! - omega_max
  !! - smooth_width_scale OR smooth_width
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_fluid_configuration(options,num_options,&
       Re,dpdl,mu_minus,mu_plus,height,scap,omega_max,smooth_width,&
       is_smooth_scaling,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,           intent(in)  :: num_options
    double precision,  intent(out) :: Re
    double precision,  intent(out) :: dpdl
    double precision,  intent(out) :: mu_minus
    double precision,  intent(out) :: mu_plus
    double precision,  intent(out) :: height
    double precision,  intent(out) :: scap
    double precision,  intent(out) :: omega_max
    double precision,  intent(out) :: smooth_width
    logical,           intent(in)  :: &
         is_smooth_scaling ! If true then treat smooth width scaling factor.
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options,num_options,opt_re,re,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_re)
       return
    end if
    call get_typed_option(options,num_options,opt_mu_minus,mu_minus,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_mu_minus)
       return
    end if
    call get_typed_option(options,num_options,opt_mu_plus,mu_plus,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_mu_plus)
       return
    end if
    call get_typed_option(options,num_options,opt_height,height,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_height)
       return
    else if ((height < 0) .or. (height > 1)) then
       ierr = option_invalid_err
       message = invalid_option_message(opt_height,height,'Value must be 0 .. 1')
       return
    end if
    call get_typed_option(options,num_options,opt_dpdl,dpdl,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_dpdl)
       return
    end if
    call get_typed_option(options,num_options,opt_scap,scap,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_scap)
       return
    end if
    call get_typed_option(options,num_options,opt_omega_max,omega_max,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_omega_max)
       return
    end if
    if (is_smooth_scaling) then
       call get_typed_option(options,num_options,opt_smooth_width_scale,smooth_width,ierr)
       if (ierr/=0) then
          message = missing_option_message(opt_smooth_width_scale)
          return
       end if
    else
       call get_typed_option(options,num_options,opt_smooth_width,smooth_width,ierr)
       if (ierr/=0) then
          message = missing_option_message(opt_smooth_width_scale)
          return
       end if
    end if
    if (show) then
       write(*,'(A,F20.10)') 'Reynolds number (Re): ', re
       write(*,'(A)') 'Viscosity of fluid: '
       write(*,'(A,F20.10)') ' Left (mu_minus): ', mu_minus
       write(*,'(A,F20.10)') ' Right (mu_plus): ', mu_plus
       write(*,'(A,F20.10)') 'Interface height: ', height
       write(*,'(A,F20.10)') 'Pressure gradient (dP/dL): ',dpdl
       write(*,'(A,F20.10)') 'Surface tension scaling parameter (scap): ', scap
       write(*,'(A,F20.10)') &
            'Maximum frequency of perturbation of interface near inlet (omega_max): ', &
            omega_max
       if (is_smooth_scaling) then
          write(*,'(A,F20.10)') 'Smooth width scale factor: ',smooth_width
       else
          write(*,'(A,F20.10)') 'Smooth width scale: ',smooth_width
       end if
    end if
  end subroutine get_fluid_configuration


  !> Get timestep.
  !! There are assumed to be options:
  !! - dt
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_timestep(options,num_options,dt,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,           intent(in)  :: num_options
    double precision,  intent(out) :: dt
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options,num_options,opt_dt,dt,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_dt)
       return
    else if (dt <= 0) then
       ierr = option_invalid_err
       message = invalid_option_message(opt_dt,dt,'Value must be > 0')
       return
    end if
    if (show) then
       write(*,'(A,F20.10)') 'Time step (dt): ', dt
    end if
  end subroutine get_timestep


  !> Get interface detection method.
  !! There are assumed to be options:
  !! - idm
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_interface_detection_method(options,num_options,&
       idm,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,           intent(in)  :: num_options
    character(len=5),  intent(out) :: idm
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options,num_options,opt_idm,idm,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_idm)
       return
    end if
    select case(idm)
    case (level_set_method)
    case (diffuse_interface_method)
    case default
       ierr = option_invalid_err
       message = invalid_option_message(opt_idm,idm,'Value must lsm or dim')
       return
    end select
    if (show) then
       write(*,'(A,A)') 'Interface detection method: ', trim(idm)
    end if
  end subroutine get_interface_detection_method


  ! Subroutines used to parse initial conditions options used as
  ! input to TPLS itself. 


  !> Get d-grid configuration.
  !! There are assumed to be options:
  !! - dx
  !! - dy
  !! - dz
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_dgrid_configuration(options,num_options,&
       dx,dy,dz,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,           intent(in)  :: num_options
    double precision,  intent(out) :: dx
    double precision,  intent(out) :: dy
    double precision,  intent(out) :: dz
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options, num_options, opt_dx, dx, ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_dx)
       return
    end if
    call get_typed_option(options, num_options, opt_dy, dy, ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_dy)
       return
    end if
    call get_typed_option(options, num_options, opt_dz, dz, ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_dz)
       return
    end if
    if (show) then
       write(*,'(A,F20.10)') 'dx: ', dx
       write(*,'(A,F20.10)') 'dy: ', dy
       write(*,'(A,F20.10)') 'dz: ', dz
    end if
  end subroutine get_dgrid_configuration


  !> Get DIM configuration.
  !! There are assumed to be options:
  !! * pe
  !! * epn
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_dim_configuration(options,num_options,Pe,epn,&
       ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,           intent(in)  :: num_options
    double precision,  intent(out) :: Pe 
    double precision,  intent(out) :: epn
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options, num_options, opt_pe, pe, ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_pe)
       return
    end if
    call get_typed_option(options, num_options, opt_epn, epn, ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_epn)
       return
    end if
    if (show) then
       write(*,'(A,F20.10)') 'DIM configuration: '
       write(*,'(A,F20.10)') ' Pe: ', PE
       write(*,'(A,F20.10)') ' epn: ', epn
    end if
  end subroutine get_dim_configuration


  !> Get DIM solver configuration.
  !! There are assumed to be options:
  !! - max_iteration_dim
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_dim_solver_configuration(options,num_options,&
    max_iteration_dim,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,           intent(in)  :: num_options
    integer,           intent(out) :: max_iteration_dim
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options, num_options, opt_max_iteration_dim,&
      max_iteration_dim, ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_max_iteration_dim)
       return
    else if (max_iteration_dim< 1) then
       ierr = option_invalid_err
       message = invalid_option_message(&
            opt_max_iteration_dim,max_iteration_dim,'Value must be >= 1')
       return
    end if
    if (show) then
       write(*,'(A)') 'DIM solver configuration: '
       write(*,'(A,I6)') ' Iterations: ', max_iteration_dim
    end if
  end subroutine get_dim_solver_configuration


  !> Get miscellaneous configuration.
  !! There is assumed to be an option:
  !! - max_u
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_misc_configuration(options,num_options,maxu,&
       ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,intent(in) :: num_options
    double precision,  intent(out) :: maxu
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options,num_options,opt_maxu,maxu,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_maxu)
       return
    end if
    if (show) then
       write(*,'(A,F20.10)') 'Maxu: ',maxu
    end if
  end subroutine get_misc_configuration


  !> Get process configuration.
  !! - There is assumed to be an option
  !!  - auto_decomp true|false.
  !! - If there is no such option then execution aborts.
  !! - If auto_decomp has value true then
  !!  - Values for num_procs_x and num_procs_y will be automatically 
  !! determined based on maxl and maxm. 
  !!  - If the difference between the ratios num_procs_x:num_procs_y 
  !! and maxl:maxm is greater than 4 then a warning is displayed
  !! to the user.
  !! - If auto_decomp has value false then there are assumed to be
  !! options: 
  !!  - num_procs_x >= 1 and divides maxl-1
  !!  - num_procs_y >= 1 and divides maxm-1
  !!  - where num_procs = num_procs_x * num_procs_y * num_procs_z
  !!  - If there are no such options then execution aborts.
  subroutine get_process_configuration(options,num_options,&
       auto_decomp,num_procs,num_procs_x,num_procs_y,num_procs_z,&
       maxl,maxm,maxn,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,intent(in)  :: num_options
    logical,intent(out) :: auto_decomp
    integer,intent(in)  :: num_procs
    integer,intent(out) :: num_procs_x
    integer,intent(out) :: num_procs_y
    integer,intent(in)  :: num_procs_z
    integer,intent(in)  :: maxl
    integer,intent(in)  :: maxm
    integer,intent(in)  :: maxn
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.
    double precision               :: mton, xtoy

    call get_typed_option(options,num_options,opt_auto_decomp,auto_decomp,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_auto_decomp)
       return
    end if
    if (auto_decomp) then
       call get_process_grid(num_procs,maxl-1,maxm-1,num_procs_x,num_procs_y)
       if (num_procs_x==0) then
          ierr = option_invalid_err
          write (message,'(A,I6,A,I6,I6,I6)') 'Error: Cannot map ',num_procs,' processes to domain ',maxl,maxm,maxn
          return
       end if
    else
       call get_typed_option(options,num_options,opt_num_procs_x,num_procs_x,ierr)
       if (ierr/=0) then
          message = missing_option_message(opt_num_procs_x)
          return
       else if (num_procs_x < 1) then
          ierr = option_invalid_err
          message = invalid_option_message(opt_num_procs_x,num_procs_x,'Value must be >= 1')
          return
       end if
       call get_typed_option(options,num_options,opt_num_procs_y,num_procs_y,ierr)
       if (ierr/=0) then
          message = missing_option_message(opt_num_procs_y)
          return
       else if (num_procs_y < 1) then
            ierr = option_invalid_err
          message = invalid_option_message(opt_num_procs_y,num_procs_y,'Value must be >= 1')
          return
       end if
       if (num_procs_x*num_procs_y*num_procs_z/=num_procs) then
          ierr = option_invalid_err
          write(message,'(A,I6,A,I6)') &
               'Error: num_procs_x*num_procs_y*num_procs_z ', &
               num_procs_x*num_procs_y*num_procs_z, &
               ' does not equal number of processes ', num_procs
          return
       else if (num_procs_z/=1) then
          ierr = option_invalid_err
          message = invalid_option_message(opt_num_procs_z,num_procs_z,'Value must be == 1')
          return
       else if (mod(maxl-1,num_procs_x)/=0) then
          ierr = option_invalid_err
          write(message,'(A,I6,A,I6)') &
               'Error: num_procs_x must divide maxl - 1. ',&
               num_procs_x,' does not divide ',maxl-1
          return
       else if (mod(maxm-1,num_procs_y)/=0) then
          ierr = option_invalid_err
          write(message,'(A,I6,A,I6)') &
               'Error: num_procs_y must divide maxm - 1. ',&
               num_procs_y,' does not divide ',maxm-1
          return
       end if
    end if ! auto_decomp
    if (show) then
       write(*,'(A)') 'Process grid: '
       write(*,'(A,I6)') ' Total: ',num_procs
       write(*,'(A,I6)') ' X: ',num_procs_x
       write(*,'(A,I6)') ' Y: ',num_procs_y
       write(*,'(A,I6)') ' Z: ',num_procs_z
       if (auto_decomp) then
          write (*,'(A)') ' TPLS selected process decomposition.'
          mton = dble(maxl) / maxm
          xtoy = dble(num_procs_x) / num_procs_y
          if (abs(mton - xtoy) > 4) then
            write (*,'(A)') ' WARNING: inefficient process decomposition!' 
            write (*,'(A,F8.2)') '  Process grid ratio: ', xtoy
            write (*,'(A,F8.2)') '  Domain grid ratio: ', mton
            write (*,'(A)') ' You may want to change the number of processes' 
            write (*,'(A)') ' or the size of the domain grid!' 
          end if
       end if
    end if
  end subroutine get_process_configuration


  !> Get solver configurations.
  !! There are assumed to be options:
  !! - max_iteration_mom_u >= 1
  !! - max_iteration_mom_v >= 1
  !! - max_iteration_mom_w >= 1
  !! - max_iteration_levelset >= 1
  !! - n_timesteps >= 1
  !! - If there are no such options then execution aborts.
  subroutine get_solver_configuration(options,num_options,&
       max_iteration_mom_u,max_iteration_mom_v,max_iteration_mom_w, &
       max_iteration_levelset,n_timesteps,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,intent(in)  :: num_options
    integer,intent(out) :: max_iteration_mom_u
    integer,intent(out) :: max_iteration_mom_v
    integer,intent(out) :: max_iteration_mom_w
    integer,intent(out) :: max_iteration_levelset
    integer,intent(out) :: n_timesteps
    integer,           intent(out)  :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options,num_options,opt_max_iteration_mom_u,max_iteration_mom_u,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_max_iteration_mom_u)
       return
    else if (max_iteration_mom_u < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(&
            opt_max_iteration_mom_u,max_iteration_mom_u,'Value must be >= 1')
       return
    end if
    call get_typed_option(options,num_options,opt_max_iteration_mom_v,max_iteration_mom_v,ierr)
    if (ierr/=0)  then
       message = missing_option_message(opt_max_iteration_mom_v)
       return
    else if (max_iteration_mom_v < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(&
            opt_max_iteration_mom_v,max_iteration_mom_v,'Value must be >= 1')
       return
    end if
    call get_typed_option(options,num_options,opt_max_iteration_mom_w,max_iteration_mom_w,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_max_iteration_mom_w)
       return
    else if (max_iteration_mom_w < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(&
            opt_max_iteration_mom_w,max_iteration_mom_w,'Value must be >= 1')
       return
    end if
    call get_typed_option(options,num_options,opt_max_iteration_levelset,max_iteration_levelset,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_max_iteration_levelset)
       return
    else if (max_iteration_levelset < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(&
            opt_max_iteration_levelset,max_iteration_levelset,'Value must be >= 1')
       return
    end if
    call get_typed_option(options,num_options,opt_n_timesteps,n_timesteps,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_n_timesteps)
       return
    else if (n_timesteps < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(&
            opt_n_timesteps,n_timesteps,'Value must be >= 1')
    end if
    if (show) then
       write(*,'(A)') 'Momentum equation solver iterations: '
       write(*,'(A,I6)') ' U: ', max_iteration_mom_u
       write(*,'(A,I6)') ' V: ', max_iteration_mom_v
       write(*,'(A,I6)') ' W: ', max_iteration_mom_w
       write(*,'(A,I6)') 'Level set equation solver iterations: ',&
            max_iteration_levelset
       write(*,'(A,I6)') 'Timesteps: ',n_timesteps
    end if
  end subroutine get_solver_configuration


  !> Get file output frequencies.
  !! There are assumed to be options:
  !! - phi_dat_frequency >=1
  !! - uvw_dat_frequency >=1
  !! - If any option is missing or has an invalid error than an error
  !! code and message are set. 
  subroutine get_output_configuration(options,num_options,backuptext, backuphdf5,&
       phi_dat_frequency,uvw_dat_frequency,backup_frequency,ierr,message,show)

    character(len=80),dimension(:,:),allocatable,intent(in) :: options
    integer,intent(in)  :: num_options
    logical,intent(out) :: backuptext
    logical,intent(out) :: backuphdf5
    integer,intent(out) :: phi_dat_frequency
    integer,intent(out) :: uvw_dat_frequency
    integer,intent(out) :: backup_frequency
    integer,           intent(out) :: ierr
    character(len=256),intent(out) :: message
    logical,           intent(in)  :: show ! Print the values.

    call get_typed_option(options,num_options, &
         opt_phi_dat_frequency,phi_dat_frequency,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_phi_dat_frequency)
       return
    else if (phi_dat_frequency < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(opt_phi_dat_frequency,phi_dat_frequency,'Value must be >= 1')
       return
    end if
    call get_typed_option(options,num_options, &
         backup_as_text, backuptext,ierr)
    if (ierr/=0) then
       message = missing_option_message(backup_as_text)
       return
    end if
    call get_typed_option(options,num_options, &
         backup_as_hdf5, backuphdf5,ierr)
    if (ierr/=0) then
       message = missing_option_message(backup_as_hdf5)
       return
    end if
    call get_typed_option(options,num_options,opt_uvw_dat_frequency,uvw_dat_frequency,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_uvw_dat_frequency)
       return
    else if (uvw_dat_frequency < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(opt_uvw_dat_frequency,uvw_dat_frequency,'Value must be >= 1')
       return
    end if
    call get_typed_option(options,num_options,opt_backup_dat_frequency,&
         backup_frequency,ierr)
    if (ierr/=0) then
       message = missing_option_message(opt_backup_dat_frequency)
       return
    else if (backup_frequency < 1) then
       ierr = option_invalid_err
       message = invalid_option_message(opt_backup_dat_frequency,backup_frequency,'Value must be >= 1')
       return
    end if
    if (show) then
       if(backuptext .and. backuphdf5) then
          write(*,*) "*************** WARNING ***************"
          write(*,*) "Printing out both netCDF and text formatted output files"
          write(*,*) "Consider choosing just one output type to improve performance"
       end if
       write(*,'(A)') 'File output frequencies: '
       write(*,'(A,I6)') ' PHI channel .dat: ',phi_dat_frequency
       write(*,'(A,I6)') ' UVW channel .dat: ',uvw_dat_frequency
       write(*,'(A,I6)') ' Backup channel .dat: ',backup_frequency
       write(*,'(A)') 'Backup file types:'
       write(*,'(A,L1)') 'netCDF HDF5 output: ',backuphdf5
       write(*,'(A,L1)') 'Text output: ',backuptext
    end if
  end subroutine get_output_configuration

end module tpls_configure
