
TPLS - Scientific Linux 6 example
=================================

The following is an example of how a user set up TPLS to run in a Scientific Linux 6.4 environment.

    $ cat /etc/redhat-release 
    Scientific Linux release 6.4 (Carbon)

The user is familiar with installing packages under Linux and has read the [TPLS - Linux user guide](./UserGuideLinux.md).

Install GFortran and MPICH2
---------------------------

To install GFortran and MPICH2, run:

    $ sudo yum install gcc-gfortran
    $ sudo yum install mpich2
    $ sudo yum list mpich2-devel

`mpich2-devel` is needed for the `mpif90` tool.

Alternatively, just run:

    $ sudo yum install mpich2-devel

This will install the following packages:

* mpich2-devel (1.2.1-2.3.el6)
* environment-modules (3.2.9c-4.el6)
* gcc-gfortran (4.4.7-3.el6)
* mpich2 (1.2.1-2.3.el6)

To check the installation has succeeded, run:

    $ gfortran -v
    gcc version 4.4.7 20120313 (Red Hat 4.4.7-3) (GCC) 
    $ mpif90 -v
    mpif90 for MPICH2 version 1.2.1
    $ cat /usr/bin/mpif90 
    sysconfdir=/etc/mpich2-x86_64
    includedir=/usr/include/mpich2-x86_64
    modincdir=/usr/include/mpich2-x86_64
    libdir=/usr/lib64/mpich2/lib
    opalibdir=/usr/lib64/mpich2/lib
    $ mpiexec.hydra -info
    HYDRA build details:
        Version:                                 1.2.1
        CC:                                      gcc
        CXX:                                     c++
        F77:                                     gfortran
        F90:                                     gfortran
    $ mpicc -v
    mpicc for MPICH2 version 1.2.1

Install PETSc 3.4.3
-------------------

Download and build PETSc rather than using yum install so that it is configured correctly for the current FORTRAN and MPI environment.  

Install BLAS and LAPACK:

    $ sudo yum install lapack
    $ sudo yum install lapack-devel

This will install the following packages:

* lapack (3.2.1-4.el6)
* blas (3.2.1-4.el6)
* lapack-devel (3.2.1-4.el6)
* blas-devel (3.2.1-4.el6)

Download PETSc 3.4.3:

    $ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.3.tar.gz
    $ tar -zxvf petsc-3.4.3.tar.gz 
    $ cd petsc-3.4.3

Configure and build in home directory:

    $ ./configure --prefix=/home/user
    $ make
    $ make install

Create new terminal window:

    $ xterm &

Run `mpd` in this window:

    $ echo "MPD_SECRETWORD=mr45-j9z" > ~/.mpd.conf
    $ chmod 600 ~/.mpd.conf
    $ mpd

In the original terminal window, run tests:

    $ make test

PETSc content should now be in the following directories in home directory:

    $ ls bin/ conf/ include/ lib/ share/

Edit `~/.bash_profile` and `~/.bashrc` files to set PETSc and library paths on login:

    export PETSC_DIR=~/petsc-3.4.3
    export LD_LIBRARY_PATH=~/lib

Now set library path:

    $ source .bash_profile

Build TPLS
----------

To build TPLS's executables, run:

    $ make linux

Ensure nodes have enough stack
------------------------------

Edit `~/.bash_profile` and `~/.bashrc` files and add:

    ulimit -s
    ulimit -s unlimited

Set up nodes to allow password-free login
-----------------------------------------

Create a key-pair:

    $ ssh-keygen -t rsa

**Always be cautious when setting up password free logins!**

    $ cp ~/.ssh/id_rsa.pub ~/.ssh/authorized_keys
    $ cat ~/.ssh/id_rsa.pub 

Log in to each node in turn and edit `~/.ssh/authorized_keys` and paste in `id_rsa.pub` content. If you want the current node to be a processing node then do this for the current node too.

Check you can now log into each other node in turn without being prompted for a password:

    $ slogin nodeNNNN

Create initial conditions for a small domain
--------------------------------------------

Set, in `initial_config.opt`:

    maxl 9
    maxm 5
    maxn 5

Run:

    $ ./create_initial_conditions

Run TPLS on a multi-processor machine
-------------------------------------

To check that it's a multi-processor machine, run:

    $ nproc
    4

The value shown is the number of processors the machine has.

Ping:

    $ mpiexec.hydra -np 4 ls
    $ mpiexec.hydra -np 4 hostname

Set, in `tpls_config.opt`:

    num_procs_x 2
    num_procs_y 2

Run TPLS:

    $ mpiexec.hydra -np 4 ./twophase.x

Run TPLS across multiple machines
---------------------------------

Write a `hosts.txt` file with the names of each node:

    node0000
    node0001
    node0002
    node0003

Ping hosts:

    $ mpiexec.hydra -np 4 -f hosts.txt hostname
    $ mpiexec.hydra -np 4 -f hosts.txt ls

Set, in `tpls_config.opt`:

    num_procs_x 2
    num_procs_y 2

Run TPLS:

    $ mpiexec.hydra -np 4 -f hosts.txt ./twophase.x
