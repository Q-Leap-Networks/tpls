
TPLS - Frequently Asked Questions
=================================

1 General
=========

1.1 Can I get the latest version of TPLS?
-----------------------------------------

TPLS is an open source project. Its source code is available from TPLS's Subversion repository on SourceForge. You can do an anonymous checkout of the latest code as follows:  

    $ svn checkout https://svn.code.sf.net/p/tpls/code/ tpls-svn

**Important note:** versions of TPLS checked out in this way are not supported in any way by the TPLS project.

1.2. What does `#` mean in options and data files?
--------------------------------------------------

Any line starting with `#` in options or data files is treated by TPLS as a comment.

2 Build problems - ARCHER
=========================

2.1 /conf/petscvariables: No such file or directory
---------------------------------------------------

If you get an error like:

    Makefile:11: /conf/petscvariables: No such file or directory
    make: *** No rule to make target `/conf/petscvariables'.  Stop.

Then you have not set `PETSC_DIR` with the location of PETSc. 

Running the module setup script:

    > source scripts/setup_modules.sh

will configure `PETSC_DIR`.

3 Build problems - Linux
========================

3.1 /conf/petscvariables: No such file or directory
---------------------------------------------------

If you get an error like:

    Makefile:11: /conf/petscvariables: No such file or directory
    make: *** No rule to make target `/conf/petscvariables'.  Stop.

Then you have not set `PETSC_DIR` with the location of PETSc. 

You should set the location of PETSc e.g.

    $ export PETSC_DIR=/home/user/petsc-3.4.3

3.2 Undefined references
------------------------

If you get errors of form:

    undefined reference to `dmcreatelocalvector_'
    undefined reference to `dmdacreate3d_'
    undefined reference to `dmdagetcorners_'
    undefined reference to `dmdagetghostcorners_'
    undefined reference to `dmdavecgetarrayf903_'
    ...
    undefined reference to `kspgetconvergedreason_'
    undefined reference to `kspgetdm_'
    undefined reference to `kspgetiterationnumber_'
    ...
    undefined reference to `matassemblybegin_'
    ...
    undefined reference to `pcgettype_'
    undefined reference to `petscfinalize_'
    undefined reference to `petscinitialize_'
    ...
    undefined reference to `vecdestroy_'

then check that PETSc libraries (`petsc`) are in your library path.

If you get errors of form:

    undefined reference to `XAllocColor'
    undefined reference to `XAllocNamedColor'
    undefined reference to `XChangeGC'
    undefined reference to `XCheckTypedEvent'
    ...

then check that X11 libraries (`X11`) are in your library path.

If you get errors of form:

    undefined reference to `dlclose'
    undefined reference to `dlerror'
    undefined reference to `dlopen'
    undefined reference to `dlsym'

then check that dynamic link libraries (`dl`) are in your library path.

If you get errors of form:

    undefined reference to `dasum_'
    undefined reference to `daxpy_'
    undefined reference to `dcopy_'
    undefined reference to `ddot_'
    undefined reference to `dgeev_'
    ...

then check that LAPACK libraries (`lapack`) are in your library path.

4 Runtime problems
==================

4.1 MPI_Cart_create failed
--------------------------

If you get an error like:

    MPIR_Cart_create(58).: Size of the communicator (24) is smaller than the size of the Cartesian topology (512)

then the number of processes in the process grid is more than the number of processes available. In the above, the total number of processes in the process grid is 512, but only 24 processes have been provided.

You should do one of:

* Change the number of processes available to the application (e.g. for the above, increase the number of processes available to 512).
* Change the total number of processes in the process grid to match the number of processes available (e.g. for the above, change the process grid to 12x2x1 or 8x3x1).

4.2 Error: Process grid processes does not equal number of processes available
------------------------------------------------------------------------------

If you get an error like:

    Error: num_procs_x*num_procs_y*num_procs_z 4 does not equal number of processes available 8

then the number of processes in the process grid is less than the number of processes available. In the above, the total number of processes in the process grid is 4, but 8 processes have been provided.

You should do one of:

* Change the number of processes available to the application (e.g. for the above, reduce the number of processes available to 8).
* Change the total number of processes in the process grid to match the number of processes available (e.g. for the above, change the process grid to 4x2x1).

4.3 Error: Number of X processes must divide number of L grid points - 1
------------------------------------------------------------------------

If you get an error like:

    Error: num_procs_x must divide maxl - 1. 12 does not divide 256.

then the process grid is not suited to the domain grid. For example, this specific error was caused by a process grid of 12x4x1 was used with a domain grid of 257x145x153, and 12 does not divide 256.

You should do one of:

* Change the process grid to divide the domain grid (e.g. for the above, change the process grid to 16x3x1 since 16 divides 256).
* Change the domain grid to be divisible by the process grid (e.g. for the above, change the domain grid to 253x145x153 since 12 divides 252).

Remember when changing process grid to ensure that it still suits the number of processes you are using (otherwise 4.1 or 4.2 above may arise).

4.4 Error: Number of Y processes must divide number of M grid points - 1
------------------------------------------------------------------------

If you get an error like:

    Error: num_procs_y must divide maxm - 1. 8 does not divide 146.

then the process grid is not suited to the domain grid. For example, this specific error was caused by a process grid of 8x8x1 was used with a domain grid of 257x149x153, and 8 does not divide 148.

You should do one of:

* Change the process grid to divide the domain grid (e.g. for the above, change the process grid to 16x4x1 since 4 divides 148).
* Change the domain grid to be divisible by the process grid (e.g. for the above, change the domain grid to 257x153x153 since 8 divides 152).

Remember when changing process grid to ensure that it still suits the number of processes you are using (otherwise 4.1 or 4.2 above may arise).

5 Runtime problems - ARCHER
===========================

5.1 apsched: claim exceeds reservation's node-count
---------------------------------------------------

If you get an error like:

    apsched: claim exceeds reservation's node-count

then you have requested more processes than the number of nodes you have requested can support e.g. if your PBS file specified:

    #PBS -l select=2

    aprun -n 64 twophase.x

then this error would arise as you have requested 2 nodes which can only provide 2*24 = 48 processes. 

You should do one of:

* Increase the number of requested nodes (e.g. for the above, request 3 nodes).
* Reduce the number of processes (e.g. for the above, reduce the number of processes to 48). 

6 Runtime problems - Linux
==========================

6.1 error while loading shared libraries: libpetsc.so
-----------------------------------------------------

If you get an error like:

    error while loading shared libraries: libpetsc.so: cannot open shared object file: No such file or directory

then the shared PETSc libraries cannot be found on one or more of your nodes. 

You should check that:

* Each node has `LD_LIBRARY_PATH` set.
* `LD_LIBRARY_PATH` includes the full path to the directory with the PETSc shared libraries.

6.2 Segmentation fault 
----------------------

If you get a:

    Segmentation fault 

then this may be due to a lack of stack space available to nodes.

For each node, edit `~/.bash_profile` and `~/.bashrc` files and add:

    ulimit -s
    ulimit -s unlimited

6.3 Unable to change wdir (No such file or directory)
-----------------------------------------------------

If you get an error like:

    [proxy@ssi-trac.epcc.ed.ac.uk] HYD_pmcd_pmi_proxy_launch_procs (./pm/pmiserv/pmi_proxy_utils.c:658): unable to change wdir (No such file or directory)
    [proxy@ssi-trac.epcc.ed.ac.uk] HYD_pmcd_pmi_proxy_control_cmd_cb (./pm/pmiserv/pmi_proxy_cb.c:111): HYD_pmcd_pmi_proxy_launch_procs returned error
    [proxy@ssi-trac.epcc.ed.ac.uk] HYDT_dmx_wait_for_event (./tools/demux/demux.c:168): callback returned error status
    [proxy@ssi-trac.epcc.ed.ac.uk] wait_for_procs_to_finish (./pm/pmiserv/pmi_proxy.c:22): demux engine error waiting for event
    [proxy@ssi-trac.epcc.ed.ac.uk] main (./pm/pmiserv/pmi_proxy.c:120): error waiting for processes to finish
    Unknown signal 120 (signal 120)

then this may be due to inconsistent paths to the working directories on your nodes, for example:

    /disk/node000/home/user/tpls-svn/trunk/src
    /disk/node001/home/user/tpls-svn/trunk/src

You should ensure that the nodes have the same paths to their working directories.
