
Auto-generated documentation and Doxygen
========================================

Auto-generated documentation is documentation automatically generated from source code and also from comments specified using a tool-specific mark-up. Auto-generated documentation provides a way for both users and developers to understand code and understand how to use components (e.g. those bundled within libraries) or to implement components conforming to an application programming interface (API) without slogging through the source code.

Doxygen
-------

[Doxygen](http://www.stack.nl/~dimitri/doxygen/) is a popular documentation generator which can auto-generate, from Fortran source code, a hyperlinked set of HTML pages allowing functions, types, modules and programs to be browsed. It also produces graphs showing the relationship between caller and called functions, for example. If comments are written in a special format, then Doxygen can parse these too and include them in the output.

For more information see:

 * [Special commands](http://www.stack.nl/~dimitri/doxygen/manual/commands.html)
 * [Comment blocks in Fortran](http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks)

Using Doxygen on ARCHER
-----------------------

Doxygen is already available. You can load it as a module and then check it is loaded as follows:

    > module load doxygen
    > doxygen -v

Using Doxygen on Linux
----------------------

If you are using Linux, then it is recommended to build Doxygen from source to get the latest version.

Check you have the prerequisite tools needed for Doxygen:

    $ flex -V
    flex 2.5.35
    $ bison -V
    bison (GNU Bison) 2.4.1
    $ dot -V
    dot - graphviz version 2.26.0 (20091210.2329)

If you do not have these then install them, for example:

    $ sudo yum install flex
    $ sudo yum install bison
    $ sudo yum install graphviz

Download and unpack Doxygen:

    $ wget http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.6.src.tar.gz
    $ tar -zxvf doxygen-1.8.6.src.tar 
    $ cd doxygen-1.8.6/

Build and install in your home directory;

    $ ./configure --prefix ~/
    $ make
    $ make install
    $ ls ~/bin
    doxygen

Add `~/bin` to your `PATH` if it is not already there:

    $ export PATH=~/bin:$PATH

Run Doxygen:

    $ doxygen --version
    1.8.6

Run Doxygen from the command-line
---------------------------------

After installing Doxygen, it can be run out-of-the-box, as follows:

    $ doxygen --version 
    1.8.6
    $ cd tpls_1.0/src/

Create a default Doxygen configuration file:

    $ doxygen -g DefaultDoxyfile

Run Doxygen using the default configuration.

    $ doxygen DefaultDoxyfile

This outputs two directories `html/` and `latex/`.

Doxygen and TPLS
----------------

TPLS code is documented using Doxygen. The Makefile contains a target that will create Doxygen for TPLS.

If you have Doxygen available then you can run:

    $ make apidoc

This outputs two directories `html/` and `latex/` with the Doxygen documentation.

Note that on ARCHER you will get errors like:

    If you installed Graphviz/dot after a previous failing run, 
    try deleting the output directory and rerun doxygen.
    error: problems opening map file 
    /fs4/e174/e174/mjje714/tpls-svn/trunk/src/html/classpressure__solver_a9908c6a6402c1ecbd9420db1742c5c7f_icgraph.map 
    for inclusion in the docs!

To disable these, edit `Doxyfile` and set the `HAVE_DOT` option value to `NO`.

TPLS-specific Doxygen settings
------------------------------

To see the difference between the default Doxygen settings and those for TPLS, run:

    $ doxygen -g DefaultDoxyfile
    $ diff Doxyfile DefaultDoxyfile

If the `dot` program has been installed then the following three settings create call graphs:
 
    HAVE_DOT               = YES
    CALL_GRAPH             = YES
    CALLER_GRAPH           = YES

The following setting:

    JAVADOC_AUTOBRIEF      = YES

means that the first sentence of any Doxygen comment is considered to be a one-line summary of the associated entity (e.g. class or function), which can then be displayed in the various summary pages created by Doxygen.
