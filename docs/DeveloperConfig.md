
Extending TPLS to support your own configuration and data
=========================================================

You can extend TPLS - both the programs to create initial conditions, `create_initial_conditions`, and to run simulations, `twophase.x` - to use configuration options and data specific to your research.

This page assumes you are familiar with [Configuring TPLS](./ConfiguringTpls.md).

How TPLS is configured
----------------------

By default TPLS is configured as follows:

`create_initial_conditions` reads:

* `initial_config.opt`, or another file specified using the `-f` command-line flag - key-value configuration options file written by the user.

`create_initial_conditions` writes:

* `initial_u.dat` - 3D data file.
* `initial_v.dat` - 3D data file.
* `initial_w.dat` - 3D data file.
* `initial_viscosity.dat` - 3D data file.
* `initial_phi.dat` - 3D data file.
* `initial_pressure.dat` - 3D data file.
* `initial_config.dat` - key-value configuration options file, auto-generated by `create_initial_conditions`.

If the `-d` command-line flag is used to specify a directory, then the files are written to that directory else they are written to the current directory.

`twophase.x` reads:

* The 7 files output by `create_initial_conditions`.
* `tpls_config.opt`, or another file specified using the `-f` command-line flag - key-value configuration options file written by the user.

If the `-d` command-line flag is used to specify a directory, then the files are assumed to be in that directory.

All these files are plain-text.

Key-value configuration options
-------------------------------

In TPLS, configuration options take the form of key-value pairs. For example:

    auto_decomp false
    maxu 10.0
    dt 0.0001
    maxn 153

Each key can be of length up to 20 characters and each value up to 60 characters.

How to read your options from `initial_config.opt` or `tpls_config.opt`
-----------------------------------------------------------------------

If you want to use your own configuration options then you can do this as follows.

For `create_initial_conditions` and `initial_config.opt`:

* Add key-value entries to the TPLS initial conditions configuration file, `initial_config.opt`.
* These will be automatically read in when `create_initial_conditions` is run.
* In `create_initial_conditions_program.F90`, they are in the `options` array.

For `twophase.x` and `tpls_config.opt`:

* Add key-value entries to the TPLS configuration options file, `tpls_config.opt`.
* These will be automatically read in when `twophase.x` is run.
* In `main_ns_hybrid.F90`, they are in the `options` array.

The subroutines and functions of `options_utils.F90` allow you to check that the values exists and to read these into variables of the desired type. For example, if you add a key-value pair:

    max_density 1.23

And want to read this in as a `double precision`, then you can do:

    integer          :: ierr
    double precision :: max_density

    call get_typed_option(options,num_options,'max_density',max_density,ierr)

For examples of reading configuration options, see the subroutines in `configure.F90`.

How to check options exist and validate options
-----------------------------------------------

`options_utils,F90` has subroutines to create error messages if there is a missing option or it has an invalid value. For example:

    if (ierr/=0) then
       write(*,'(A)') missing_option_message('max_density')
    else if ((max_density < 0) .or. (max_density > 1)) then
       write(*,'(A)') invalid_option_message('max_density',max_density,'Value must be 0 .. 1')
    end if

For examples of validating configuration options, see the subroutines in `configure.F90`.

How to load your own options file
---------------------------------

You can use your own configuration options file e.g. `my_config.opt` with your own key-value entries. You can then load this into an array using subroutines of `options_utils.F90`. For example:

    character(len=80),dimension(:,:),allocatable :: my_options
    integer                                      :: num_my_options
    integer,                                     :: ierr
    logical                                      :: file_exists
    character(len=256)                           :: filename

    filename='my_config.opt'
    write(*,'(A,A)') 'Loading options from options file: ',filename
    call load_options(filename,my_options,num_my_options,ierr)
    if (ierr==option_file_not_found_err) then
       write(*,'(A,A)') 'Error: Options file not found: ',filename
    else if (ierr/=0) then
       write(*,'(A,A)') 'Error: Loading options file: ',filename
    end if

How to modify `create_initial_conditions` to save your options
--------------------------------------------------------------

`create_initial_conditions_program.F90` saves the options used to create the initial conditions, and other values it has created, in a file called `initial_config.dat`. `initial_config.dat` is one of the files that `twophase.x` loads in by default.

`initial_config.dat` is populated as follows:

    ! Create initial_conditions with values created by get_initial_twophase.

    allocate(initial_conditions(20,2))
    num_initial_conditions=0
    call put_typed_option(initial_conditions,num_initial_conditions,opt_maxl,maxl)
    call put_typed_option(initial_conditions,num_initial_conditions,opt_maxm,maxm)
    call put_typed_option(initial_conditions,num_initial_conditions,opt_maxn,maxn)
    ...

You can save additional options into `initial_config.dat` by just adding additional calls to `put_typed_option`. For example:
    
    call put_typed_option(initial_conditions,num_initial_conditions,'density_tolerance',0.1249)

You should also edit the `allocate` line to ensure the array is big enough to hold all the options put into it. So, if we add one option we'd change the array size from 20 to 21:

    allocate(initial_conditions(21,2))

How to save your own options file
---------------------------------

You can save your own options file using subroutines in `option_utils.F90`, and in a similar way to that used by `create_initial_conditions_program.F90`. For example:

    character(len=256),dimension(:),allocatable   :: header
    integer                                       :: num_header_lines
    character(len=80), dimension(:,:),allocatable :: conditions
    integer                                       :: num_conditions

    allocate(conditions(5,2))
    num_conditions=0
    call put_typed_option(conditions,num_conditions,'density_tolerance',0.1249)
    call put_typed_option(conditions,num_conditions,'density_scaled',.false.)
    call put_typed_option(conditions,num_conditions,'density_x',0)
    call put_typed_option(conditions,num_conditions,'density_y',1)
    call put_typed_option(conditions,num_conditions,'density_z',2)

    num_header_lines=1
    allocate(header(1))
    header(1)='My initial conditions configuration'

    call get_filepath(config_dir,'density_config.dat',filename))
    call save_options(filename,conditions,num_conditions,header,num_header_lines)

`options_utils.F90` also has an `append_options` subroutine which allows options to be appended to an existing file.

In both `create_initial_conditions_program.F90` and `main_ns_hybrid.F90`, the `config_dir` variable holds the name of the directory into which files are saved. `io.F90` has a `get_filepath` subroutine that will construct a file path given a directory name and a file name.

How to modify `create_initial_conditions` to save your own data files
---------------------------------------------------------------------

If you have additional data files representing your initial conditions, then how you save these is up to you.

The `config_dir` variable holds the name of the directory into which files are saved. Whether you use this is up to you.

`io.F90` has subroutines that can be used to save 3D data which you may find of use.

In `create_initial_conditions_program.F90`, the `config_dir` variable holds the name of the directory into which files are saved. `io.F90` has a `get_filepath` subroutine that will construct a file path given a directory name and a file name.

How to modify `twophase.x` to read your options from `initial_config.opt`
-------------------------------------------------------------------------

If you modify `create_initial_conditions_program.F90` to save your own configuration options in `initial_config.opt`, then these will be automatically read in when `twophase.x` is run. In `main_ns_hybrid.F90`, they are in the `initial_conditions` array.

The subroutines and functions of `options_utils.F90` allow you to check that the values exists and to read these into variables of the desired type, as described in "How to read your options from `initial_config.opt` or `tpls_config.opt`" above.

For examples of reading configuration options and validating them, see `main_ns_hybrid.F90` and the subroutines in `configure.F90`.

How to modify `twophase.x` to load your own data files
------------------------------------------------------

If you have additional data files, then how you load these is up to you.

The `config_dir` variable holds the name of the directory from which files are loaded. Whether you use this is up to you.

`io.F90` has subroutines that can be used to load 3D data which you may find of use.

In `main_ns_hybrid.F90`, the `config_dir` variable holds the name of the directory into which files are saved. `io.F90` has a `get_filepath` subroutine that will construct a file path given a directory name and a file name.

How to read command-line options
--------------------------------

`options_utils.F90` has a `get_command_line_options` subroutine to read in command-line options provided to both `create_initial_conditions` and `twophase.x`. This subroutine returns a 2D array of key-values that can be used as for options read from options files. For example:

    character(len=80),dimension(:,:),allocatable :: options
    integer                                      :: num_options

    call get_command_line_options(options,num_options)
    if (has_option(options,num_options,'-accelerate')) then
      ...
    else
      ...
    end if

Why you should write your own `option_names`-style modules
----------------------------------------------------------

If you are creating many of your own options, then you should create your own module similar to that of `option_names.F90`, with constants for each option. For example: `density_option_names.F90` might hold:

    character(*),parameter :: opt_max_density='max_density'
    character(*),parameter :: opt_density_scaled='density_scaled'

You can then use these elsewhere in your code. For example:

    call get_typed_option(options,num_options,opt_max_density,max_density,ierr)

    call put_typed_option(my_conditions,num_my_conditions,opt_density_scaled,.false.)

If you decide to change the name of an option, you only need to change it in your module where the constant is defined rather than in every location in your code where you use it.
