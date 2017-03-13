!> TPLS configuration and command-line option names.
!!
!! @author Mike Jackson.
!! @version $Revision: 328 $
!! @copyright (c) 2015, The University of Edinburgh, all rights
!! reserved. 
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module option_names

  implicit none

  !> Default TPLS initial u velocity file name.
  character(*),parameter :: initial_u_file='initial_u.dat'
  !> Default TPLS initial velocity file name.
  character(*),parameter :: initial_v_file='initial_v.dat'
  !> Default TPLS initial w velocity file name.
  character(*),parameter :: initial_w_file='initial_w.dat'
  !> Default TPLS initial pressure file name.
  character(*),parameter :: initial_pressure_file='initial_pressure.dat'
  !> Default TPLS level-set file name.
  character(*),parameter :: initial_phi_file='initial_phi.dat'
  !> Default TPLS viscosity file name.
  character(*),parameter :: initial_viscosity_file='initial_viscosity.dat'
  !> Default TPLS settings file name.
  character(*),parameter :: initial_config_file='initial_config.dat'

  !> Command-line flag for configuration files directory.
  character(*),parameter :: opt_config_dir="-d"
  !> Command-line flag for options file.
  character(*),parameter :: opt_file="-f"

  !> Method names and indices - level set method.
  character(len=5),parameter :: level_set_method="lsm"
  !> Method names and indices - diffuse interface method (DIM).
  character(len=5),parameter :: diffuse_interface_method="dim"

  !> Domain grid - key for number of grid points in X (l) dimension.
  character(*),parameter :: opt_maxl="maxl"
  !> Domain grid - key for number of grid points in Y (m) dimension.
  character(*),parameter :: opt_maxm="maxm"
  !> Domain grid - key for number of grid points in Z (n) dimension.
  character(*),parameter :: opt_maxn="maxn"

  !> Process decomposition - key for number of processes in X dimension.
  character(*),parameter :: opt_num_procs_x="num_procs_x"
  !> Process decomposition - key for number of processes in Y dimension.
  character(*),parameter :: opt_num_procs_y="num_procs_y"
  !> Process decomposition - key for number of processes in Z dimension.
  character(*),parameter :: opt_num_procs_z="num_procs_z"
  !> Process decomposition - key for auto-selection of process decomposition.
  character(*),parameter :: opt_auto_decomp="auto_decomp"

  !> Momentum equation solver - key for number of u velocity iterations.
  character(*),parameter :: opt_max_iteration_mom_u="mom_u"
  !> Momentum equation solver - key for number of v velocity iterations.
  character(*),parameter :: opt_max_iteration_mom_v="mom_v"
  !> Momentum equation solver - key for number of w velocity iterations.
  character(*),parameter :: opt_max_iteration_mom_w="mom_w"
  !> Level-set equation solver - key for number of level-set iterations.
  character(*),parameter :: opt_max_iteration_levelset  ="levelset"

  !> Interface detection method.
  character(*),parameter :: opt_idm="idm"

  !> Fluid flow - key for Reynolds number.
  character(*),parameter :: opt_re="re"
  !> Fluid flow - key for viscosity of the lower fluid.
  character(*),parameter :: opt_mu_minus="mu_minus"
  !> Fluid flow - key for viscosity of the upper fluid.
  character(*),parameter :: opt_mu_plus="mu_plus"
  !> Fluid flow - key for interface height.
  character(*),parameter :: opt_height="height"
  !> Fluid flow - key for pressure gradient.
  character(*),parameter :: opt_dpdl="dpdl"
  !> Fluid flow - key for surface tension scaling parameter.
  character(*),parameter :: opt_scap="scap"
  !> Fluid flow - key for maximum frequency of perturbation of interface
  !! near inlet.
  character(*),parameter :: opt_omega_max="omega_max"
  !> Fluid flow - key for smooth width scale factor.
  character(*),parameter :: opt_smooth_width_scale="smooth_width_scale"
  !> Fluid flow - key for smooth width.
  character(*),parameter :: opt_smooth_width="smooth_width"
  !> Fluid flow - key for maximum u.
  character(*),parameter :: opt_maxu="maxu"

  !> Grid - key for dx.
  character(*),parameter :: opt_dx="dx"
  !> Grid - key for dy.
  character(*),parameter :: opt_dy="dy"
  !> Grid - key for dz.
  character(*),parameter :: opt_dz="dz"
  !> Grid - key for time step.
  character(*),parameter :: opt_dt="dt"

  !> DIM - key for Pe.
  character(*),parameter :: opt_pe="pe"
  !> DIM - key for epn.
  character(*),parameter :: opt_epn="epn"
  !> DIM - key for DIM interations.
  character(*),parameter :: opt_max_iteration_dim="max_iteration_dim"

  !> TPLS operation - key for PHI channel .dat file output frequency.
  character(*),parameter :: opt_phi_dat_frequency="phi_dat_frequency"
  !> TPLS operation - key for UVW channel .dat file output frequency.
  character(*),parameter :: opt_uvw_dat_frequency="uvw_dat_frequency"
  !> TPLS operation - key for backup channel .dat file output frequency.
  character(*),parameter :: opt_backup_dat_frequency="backup_frequency"
  !> TPLS operation - key for requiring text backup format.
  character(*),parameter :: backup_as_text="backup_text_format"
  !> TPLS operation - key for requiring netCDF-hdf5 backup format.
  character(*),parameter :: backup_as_hdf5="backup_hdf5_format"
  !> TPLS operation - key for number of timesteps (> 0).
  character(*),parameter :: opt_n_timesteps="num_timesteps"

end module option_names
