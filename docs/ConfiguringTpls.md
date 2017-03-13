Configuring TPLS
================

TPLS is configured and run as follows:

* Create initial conditions files
    * EITHER
        * Edit initial conditions configuration file.
        * Run TPLS initial conditions tool, `create_initial_conditions`.
    * OR
        * Create valid initial conditions files by hand.
* Edit TPLS configuration file.
* Run TPLS, `twophase.x`.

Initial conditions
==================

The TPLS executable, `twophase.x` reads in data files represent the domain at time t=0, the initial conditions.

These files, and their default names, are as follows (assuming an LxMxN domain):

* `initial_u.dat` - initial L-dimension velocities for points in range (0..L),(0..M),(0..N-2).
* `initial_v.dat` - initial M-dimension velocities for points in range (0..L),(0..M),(0..N-2).
* `initial_w.dat` - initial N-dimension velocities for points in range (0..L),(0..M),(0..N-1).
* `initial_viscosity.dat` - initial viscosities for points in range (0..L),(0..M),(0..N).
* `initial_phi.dat` - initial phi (level-set) function values for points in range (0..L),(0..M),(0..N).
* `initial_pressure.dat` - initial pressures for points in range (0..L),(0..M),(0..N).
* `initial_config.dat` - configuration values used to generate the foregoing data files:
    * Domain dimensions (L,M,N).
    * Grid point widths in each dimension.
    * Reynolds number.
    * Pressure gradient.
    * Viscosity of fluid one grid cell to left and right.
    * Interface height.
    * Smooth width.

Initial conditions file formats
-------------------------------

Each file is assumed to be in plain-text format.

To see the expected format of the initial conditions files, run `create_initial_conditions` (described below) and look at the output data files.

Initial conditions file location
--------------------------------

At present, `twophase.x` assumes the initial conditions files are located in the same directory as itself.

TPLS initial conditions tool
============================

The TPLS initial conditions tool, `create_initial_conditions`, can auto-generate initial conditions files. 

`create_initial_conditions` needs an initial conditions options file to run. This file specifies the domain dimensions, fluid flow, pressures etc.

Create initial conditions
-------------------------

Edit initial conditions options file, `initial_config.opt`.

Run `create_initial_conditions` to create initial conditions files:

    $ ./create_initial_conditions

By default, `create_initial_conditions` assumes `initial_config.opt` is located in the same directory. If you want to use an initial conditions options file with a different name and/or located elsewhere then use the `-f` flag. For example:

    $ ./create_initial_conditions -f /home/user/expr345_initial_config.opt

By default, the initial conditions files are written to the current directory. If you want you output the initial conditions files to another directory then use the `-d` flag. For example:

    $ ./create_initial_conditions -f /home/user/expr345_initial_config.opt -d /home/user/tpls_conditions

The directory must already exist.

Available initial conditions options
------------------------------------

`initial_config.opt` provides a full set of initial conditions options and default values. The file is commented to explain what each option is.

The default values match the default values hard-coded into `create_initial_conditions`.

Using your own initial conditions files
---------------------------------------

TPLS does not require you to use `create_initial_conditions` to create your initial conditions files. You can create your files any way you wish so long as they are of the required format.

TPLS configuration
==================

Asides from the initial conditions files, the TPLS executable, `twophase.x` needs a TPLS options file to run. This file specifies the process decomposition, number of iterations, and other options for controlling TPLS.

Available TPLS options
----------------------

`tpls_config.opt`, provides a full set of the TPLS options and default values. The file is commented to explain what each option is.

The default values match the default values hard-coded into `twophase.x`.

Configure TPLS
--------------

Create initial conditions files as described above.

Edit the TPLS options file, `tpls_config.opt`.

Run `twophase.x`. For example, if using a PBS file on ARCHER:

    aprun -n 24 twophase.x

For example, if using `mpiexec.hydra`, on a Linux MPI-enabled network:

    $ mpiexec.hydra -np 24 -f hosts6x4.txt ./twophase.x

By default, `twophase.x` assumes `tpls_config.opt` is located in the same directory. If you want to use a TPLS options file with a different name and/or located elsewhere then use the `-f` flag. For example:

    aprun -n 24 twophase.x -f /home/user/tpls.opt

    $ mpiexec.hydra -np 24 -f hosts6x4.txt ./twophase.x -f /home/user/tpls.opt

By default, the initial conditions files are read from the current directory. If you want you input the initial conditions files from another directory then use the `-d` flag. For example:

    aprun -n 24 twophase.x -f /home/user/tpls.opt -d /home/user/tpls_conditions

    $ mpiexec.hydra -np 24 -f hosts6x4.txt ./twophase.x -f /home/user/tpls.opt -d /home/user/tpls_conditions

Interface detection method
==========================

TPLS supports two interface detection methods:

* Level-set method (default).
* Diffuse interface method (DIM).

The method to use is set in the initial conditions options file. 

To use the level-set method, set the option:

    idm lsm

To use the diffuse inferface method (DIM), set the option:

    idm dim

Domain and process grids
========================

TPLS represents its domain by an LxMxN grid. This has a default size of 257x145x153.

This grid is mapped onto a process grid, with dimensions XxYxZ, where Z is 1. This has a default size of 32x16x1. 

The domain grid and process grid must be defined such that:

* X must divide L-1 exactly. For example, 32 divides (257 - 1).
* Y must divide M-1 exactly. For example, 16 divides (145 - 1).

The process grid must be defined such that:

* X * Y * Z = number of processes requested to run TPLS. For example, for the default 32x16x1, 512 processes must be requested in total.

Process grids and ARCHER
------------------------

For ARCHER, each node has 24 cores, so the number of nodes requested in the PBS file must be at least number of processes / 24 (rounded up to the nearest integer).

Examples
--------

For example, if using 24 processes then a valid process grid for the default domain grid is 8x3, and this can be configured by setting, within the TPLS options file:

    num_procs_x 8
    num_procs_y 3

For example, if using 8 processes then a valid process grid for a small domain grid of 9x5x5 is 4x2, and this can be configured via:

    num_procs_x 4
    num_procs_y 2

**Important note:** TPLS cannot be run with just a single process at this time - a 1x1x1 process grid will not work.

TPLS and auto-selection of process grids
========================================

Using the number of processes available and the L and M dimensions of the domain grid (as specified by the `maxl` and `maxm` options), TPLS can make a best guess on a suitable process grid.

To request that TPLS make a best guess, set the `auto_decomp` option in the TPLS options file to `true`:

    auto_decomp true

If `auto_decomp` is set to `true` then TPLS will ignore any values for `num_procs_x` and `num_procs_y` and TPLS will display the following message when run:

    TPLS selected process decomposition.

If TPLS cannot automatically guess a valid process decomposition, then it will exit with a message like:

    Error: Cannot map 32 processes to domain 255 143 153

To see what process grid TPLS will choose use the `grids` tool with `ptodgrid` (see below).

Select process and grids using `grids`
======================================

TPLS comes with a helper tool, `grids`, which can help you to identify a valid process and domain grid combination based on the number of processes you want to use.

To see the commands supported by `grids`, run:

    $ ./grids

Select a valid process grid
---------------------------

Given an LxM domain grid and N processors, to see the valid process grids, run:

    $ ./grids ptodgrids <N> <L> <M>

For example, valid grids for 512 processes and the default TPLS domain include:

    $ ./grids ptodgrids 512 257 145
    Process grids for placing a    257 by    145 domain grid onto    512 processes:
       256     2
       128     4
        64     8
        32    16

For example, valid grids for 24 processes and the default TPLS domain include:

    $ ./grids ptodgrids 24 257 145
    Process grids for placing a    257 by    145 domain grid onto     24 processes:
         1    24
         2    12
         8     3
         4     6

The N dimension of an LxMxN domain grid is not referred to since in process grids XxYxZ, TPLS always assumes Z=1, so N can be any value >= 1.

Select the best process grid
----------------------------

Given an LxM domain grid and N processors, to see the best process grid, run:

    $ ./grids ptodgrid <N> <L> <M>

The best grid, XxY is chosen such that, of all the grids possible, the ratio X:Y is closest to the ratio L:M.

For example, the best grid for 512 processes and the default TPLS domain is:

    $ ./grids ptodgrid 512 257 145
    Best process grid for placing a    257 by    145 domain grid onto    512 processes:
        32    16

For example, the best grid for 24 processes and the default TPLS domain is:

    $ ./grids ptodgrid 24 257 145
    Best process grid for placing a    257 by    145 domain grid onto     24 processes:
         8     3

The N dimension of an LxMxN domain grid is not referred to since in process grids XxYxZ, TPLS always assumes Z=1, so N can be any value >= 1.

Select a valid domain grid
--------------------------

Instead of matching process grids to domain grids, you can choose a domain grid based on a process grid.

To see all possible XxYx1 process grids given N processes, run:

    $ ./grids pgrids <N>

For example:

    $ ./grids pgrids 512
    Process grids for    512 processes:
       512     1
       256     2
       128     4
        64     8
        32    16
         1   512
         2   256
         4   128
         8    64
        16    32

Given an XxY process grid, to see N valid domain grids LxM that TPLS would accept, run:

    $ ./grids dgrids <X> <Y> <N>

For example, 10 domain grids that would map onto a 32x16 process grid are:

    $ ./grids dgrids 32 16 10
    Valid domain grids for a     32 by     16 process grid:
        33    17 =    561
        65    33 =   2145
        97    49 =   4753
       129    65 =   8385
       161    81 =  13041
       193    97 =  18721
       225   113 =  25425
       257   129 =  33153
       289   145 =  41905
       321   161 =  51681

And, to validate, using the `ptodgrids` option:

    $ ./grids ptodgrids 512 321 161
    Process grids for placing a    321 by    161 domain grid onto    512 processes:
        64     8
        32    16
        16    32

we can see that 32x16 is offered as an option.
