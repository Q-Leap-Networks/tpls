
TPLS - Linux user guide
=======================

TPLS can run on an MPI-enabled Linux cluster or network.

Prerequisites
-------------

You will need:

* A FORTRAN compiler. At present the GNU Fortran and Cray compilers can be used.
* An MPI library.
* PETSc 3.5.2.1. PETSc can be downloaded from the [PETSc web site](http://www.mcs.anl.gov/petsc/).
    * Due to API changes within PETSc, earlier versions will **not** work.
    * The `PETSC_DIR` environment variable should be set to point to the PETSc 3.4.3 distribution. If this has been installed via a package manager, then this variable should be set to the parent of the `conf` directory holding the `petscvariables` script (this is used by the TPLS Makefile). Common locations of this directory are `/usr/` or `/usr/local`.
    * If using shared libraries, the `LD_LIBRARY_PATH` environment variable should point to the `lib` directory containing the PETSc libraries.
* NetCDF 4.3.2, NetCDF-Fortran 4.2, HDF5 1.8.12
    * Note that more recent versions of HDF5 are not compatible with NetCDF 4.3.2
    * The `NETCDF_FORTRAN_DIR` environment variaible should be set to point to the directory where NetCDF Fortran is installed (should contain `/lib` and `/include` directories).
    * The `NETCDF_DIR` environment variaible should be set to point to the directory where NetCDF is installed (should contain `/lib` and `/include` directories).
    * The `HDF5_DIR` environment variaible should be set to point to the directory where HDF5 is installed (should contain `/lib` and `/include` directories).
* A network which has been configured to run MPI applications.

Download TPLS
-------------

If you have not already done so, download TPLS from [SourceForge](http://sourceforge.net/projects/tpls/). 

UnZIP TPLS:

    $ tar -zxvf tpls.tar.gz 
    $ cd src/

Build TPLS
----------

To build TPLS, run:

    $ make linux

This builds:

* TPLS initial conditions tool, `create_initial_conditions`.
* TPLS executable, `twophase.x`.
* MPI decomposition helper tool, `grids`.
* Test driver program, `run_tests` 

For troubleshooting build problems, see the [FAQ](./Faq.md).

Create initial conditions
-------------------------

To create a set of initial conditions, run:

    $ ./create_initial_conditions

For more on configuring TPLS, see [Configuring TPLS](./ConfiguringTpls.md).

Run TPLS
--------

To run TPLS, use the MPI application execution tool for your MPI library.

For example, using `mpiexec.hydra` with 128 nodes, each with 4 cores (a total of 512 cores) you could create a hosts file, `hosts.txt`, with content:

    node000:4
    node001:4
    node002:4
    node003:4
    ...
    node127:4

and then run the following, which requests 512 processes:

    $ mpiexec.hydra -np 512 -f hosts.txt ./twophase.x

When setting the number of processes, the following conditions must be respected:

* The process grid XxYxZ (default 32x16x1) must be such that X * Y * Z = number of processes.
* Given a domain grid LxMxN (default 257x145x153) then X must divide L-1 exactly and Y must divide M-1 exactly.

So, for example, to use 6 nodes, (24 cores), with the default domain, you could create a hosts file, `hosts6x4`, with content:

    node000:4
    node001:4
    node002:4
    node003:4
    node004:4
    node005:4

and, you would set, in `tpls_config.opt`:

    num_procs_x 8
    num_procs_y 3

and then run:

    $ mpiexec.hydra -np 24 -f hosts6x4.txt ./twophase.x

For more on configuring TPLS, see [Configuring TPLS](./ConfiguringTpls.md).

For troubleshooting runtime problems, see the [FAQ](./Faq.md).
