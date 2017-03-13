
TPLS - ARCHER user guide
========================

TPLS has been developed to run on ARCHER. ARCHER has all the packages required to run TPLS.

Change to your work directory
-----------------------------

TPLS should be built and run in your work directory. To get into your work directory, once you have logged in to ARCHER, do:

    > cd /work/PROJECTCODE/PROJECTCODE/USER/

where `USER` is your user name and `PROJECTCODE` is your project code.

Download TPLS
-------------

If you have not already done so, download TPLS from [SourceForge](http://sourceforge.net/projects/tpls/). 

UnZIP TPLS:

    > tar -zxvf tpls.tar.gz 
    > cd src/

Set up modules
--------------

TPLS uses the GNU FORTRAN compiler and PETSc 3.5.2.1. To set up your environment to use these modules, run:

    > source scripts/setup_modules.sh

This does:

    module swap PrgEnv-cray PrgEnv-gnu
    module load cray-petsc/3.5.2.1
    module load cray-netcdf-hdf5parallel/4.3.2

This only has to be done once when you log in.

Build TPLS
----------

To build TPLS, run:

    > make archer

This builds:

* TPLS initial conditions tool, `create_initial_conditions`.
* TPLS executable, `twophase.x`.
* MPI decomposition helper tool, `grids`.
* Test driver program, `run_tests`

For troubleshooting build problems, see the [FAQ](./Faq.md).

Create initial conditions
-------------------------

To create a set of initial conditions, run:

    > ./create_initial_conditions

For more on configuring TPLS, see [Configuring TPLS](./ConfiguringTpls.md).

Create a batch job script
-------------------------

`scripts/tpls.pbs` contains a PBS batch script. You can edit this to create a batch job to run TPLS on ARCHER.

Copy `scripts/tpls.pbs`:

    > cp scripts/tpls.pbs .

Open `tpls.pbs` in an editor.

Provide a name for your job. Edit the line:

    #PBS -N tpls

Set the maximum wall-time for your job (the time your job will be allowed to run before it's cancelled). The default time is 10 minutes. Edit the line:

    #PBS -l walltime=0:10:00

Set the project code for your job. This is the budget to which the job is charged. Edit the following line and replace `z0` with your project code:

    #PBS -A z01

Set the total number of nodes for your job. Edit the following line:

    #PBS -l select=22

Set the number of processes for your job. Edit the following line and replace `512` with the number of processes:

    aprun -n 512 twophase.x

When setting the number of nodes and processes, the following conditions must be respected:

* Each node on ARCHER has 24 cores so number of processes must be <= (number of nodes * 24).
* The process grid XxYxZ (default 32x16x1) must be such that X * Y * Z = number of processes.
* Given a domain grid LxMxN (default 257x145x153) then X must divide L-1 exactly and Y must divide M-1 exactly.

So, for example, to use 1 node, with the default domain, you would set, in `tpls_config.opt`:

    num_procs_x 8
    num_procs_y 3

and in `tpls.pbs`:

    #PBS -l select=1

    aprun -n 24 twophase.x

For more on configuring TPLS, see [Configuring TPLS](./ConfiguringTpls.md).

Submit batch job
----------------

To submit your batch job, run:

    > qsub tpls.pbs

Your job's number will be displayed.

Check the status of your job
----------------------------

To query your job, run the following, where `USER` is your user name:

    > qstat -u USER

The job queue shows your job status in an `S` column where the status will be one of:

* `S`tatus:
* `Q`ueued
* `R`unning
* `E`xiting
* `F`inished

To find out more details, run the following, where `JOBNUMBER` is the job name (e.g. `227632.sdb`):

    > qstat -f JOBNUMBER.sdb

To find out the estimated start time, run the following, where `JOBNUMBER` is the job name (e.g. `227632.sdb`):

    > qstat -T JOBNUMBER.sdb

When your job completes
-----------------------

When your job completes you will have two status files in your directory, where `JOBNAME` is the name of your job as specified in your PBS file and `JOBNUMBER` is the job number:

* `JOBNAME.oJOBNUMBER` - output file. For example, `tpls.o227632`.
* `JOBNAME.eJOBNUMBER` - error file. For example, `tpls.e227632`.

Troubleshooting
---------------

For troubleshooting runtime problems, see the [FAQ](./Faq.md).
