TPLS - Two-phase Level-Set
==========================

Welcome to TPLS 2.0 (revision $Revision: 328 $).

Getting Started
---------------

Linux users, please see [Linux User Guide](./docs/UserGuideLinux.md). 

ARCHER users, please see [ARCHER User Guide](./docs/UserGuideArcher.md).

Other pages:

* [Configuring TPLS](./docs/ConfiguringTpls.md) - information on how to configure TPLS.
* TPLS options files:
    * [initial_config.opt](./initial_config.opt) - example of initial conditions options file that can be used to create initial conditions for TPLS. 
    * [tpls_config.opt](./tpls_config.opt) - example of TPLS options file that can be used to configure TPLS for running.
    * Each file documents the available options and provides the default values used.
* [Extending TPLS to support your own configuration and data](./docs/DeveloperConfig.md)
* [FAQ](./docs/Faq.md) - frequently asked questions and troubleshooting.
* [Auto-generated documentation and Doxygen](./docs/Doxygen.md)

Copyright and licence
---------------------

TPLS is Copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, Iain Bethune, Toni Collis, David Scott, Peter Spelt, The University of Edinburgh.

This program is distributed under the BSD License

See the [license](./LICENSE.txt) for details.

Research based on TPLS should cite:

* Lennon O'Naraigh, Prashant Valluri, David Scott, Iain Bethune, and Peter D. M. Spelt (2013) "Linear instability, nonlinear instability, and ligament dynamics in three-dimensional laminar two-layer liquid/liquid flows". Accepted for publication in the Journal of Fluid Mechanics.
* PETSc. As TPLS uses [PETSc](http://www.mcs.anl.gov/petsc/) you should also see [How to cite PETSc in your publications](http://www.mcs.anl.gov/petsc/documentation/referencing.html).

Third-party code
----------------

TPLS distributions contain code written by third-parties. This code and its licences is in `thirdparty/`. The code has been imported as-is without any modifications.

The third-party code is as follows:

* FRUIT - FORTRAN Unit Test Framework
    * Directory: `thirdparty/fruit/`
    * License: BDS license
    * Web site: http://sourceforge.net/projects/fortranxunit/
    * Version: 3.3.4, downloaded 02/06/14
    * Used by: test code in `test/`

Contact Us
----------

TPLS is an open source project hosted on SourceForge. Visit our site for source code and our issue tracker:

    http//sourceforge.net/projects/tpls
