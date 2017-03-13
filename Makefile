# TPLS Makefile
#
# $Date:: 2015-02-26 22:16:56 +0000 (Thu, 26 Feb 2015)     $
# $Author:: ibethune                                       $
# $Revision:: 322                                          $

FRUIT_SRC=thirdparty/fruit/fruit_util.f90 thirdparty/fruit/fruit.f90

SRC_FILES=mpi.f90 mpi_error_check.F90 weno.f90 advect_phi.f90 levelset.f90 two_phase_levelset.f90 cahn_hilliard_solver.f90 momentum_allflux.f90 jacobi_iteration_allflux.f90 sor_iteration_allflux.f90 pressure.f90 state.F90 pressure_solver.F90 option_names.F90 options_utils.F90 io.F90 grid_utils.F90 configure.F90 configure_mpi_petsc.F90 main_ns_hybrid.F90
INITIAL_CONDITIONS_SRC_FILES=option_names.F90 options_utils.F90 error_check.F90 io.F90 grid_utils.F90 configure.F90 twophase_initialisation_wave.F90 create_initial_conditions_program.F90
GRIDS_SRC_FILES=grid_utils.F90 grids_program.F90

# Source for which tests exist, or which are needed by tests.
TEST_SRC_FILES=option_names.F90 options_utils.F90 io.F90 grid_utils.F90 io.F90
TESTS=test/tpls_fruit_utils.F90 test/grid_utils_test.f90 test/options_utils_test.f90 test/io_test.f90 test/tpls_test_driver.F90

SRC=$(addprefix src/, $(SRC_FILES))
GRIDS_SRC=$(addprefix src/, $(GRIDS_SRC_FILES))
INITIAL_CONDITIONS_SRC=$(addprefix src/, $(INITIAL_CONDITIONS_SRC_FILES))
TEST_SRC=$(addprefix src/, $(TEST_SRC_FILES))

TPLS_EXE=twophase.x 
INITIAL_CONDITIONS_EXE=create_initial_conditions
GRIDS_EXE=grids
TEST_EXE=run_tests

include $(PETSC_DIR)/conf/petscvariables

NETCDF_INCLUDES=-I$(NETCDF_FORTRAN_DIR)/include
NETCDF_LIB=-L$(NETCDF_FORTRAN_DIR)/lib -lnetcdff -L$(NETCDF_DIR)/lib -lnetcdf -L$(HDF5_DIR)/lib -lhdf5_hl -lhdf5 -lz

archer : $(addsuffix .archer, $(GRIDS_EXE) $(TPLS_EXE) $(INITIAL_CONDITIONS_EXE) $(TEST_EXE))

twophase.x.archer : $(SRC)
	ftn -craype-verbose -O3 -std=f2003 -fopenmp $^ -o $(TPLS_EXE)

create_initial_conditions.archer : $(INITIAL_CONDITIONS_SRC)
	ftn -O3 -std=f2003 $^ -o $(INITIAL_CONDITIONS_EXE)

grids.archer : $(GRIDS_SRC)
	ftn -O3 -std=f2003 $^ -o $(GRIDS_EXE)

run_tests.archer : $(TEST_SRC) $(TESTS)
	ftn -O3 -std=f2003 $(FRUIT_SRC) $^ -o $(TEST_EXE)

linux : $(addsuffix .linux, $(GRIDS_EXE) $(TPLS_EXE) $(INITIAL_CONDITIONS_EXE) $(TEST_EXE))

twophase.x.linux : $(SRC)
	mpif90 -O3 -std=f2003 -fopenmp $^ $(PETSC_FC_INCLUDES) $(NETCDF_INCLUDES) -o $(TPLS_EXE) $(PETSC_LIB) $(NETCDF_LIB)
#	mpif90 -O3 -Wall -std=f2003 -fopenmp $^ $(PETSC_FC_INCLUDES) -o $(TPLS_EXE) $(PETSC_LIB)

twophase.x.linux.debug0 : $(SRC)
	mpif90 -g -O0 -fopenmp $^ $(PETSC_FC_INCLUDES) $(NETCDF_INCLUDES) -o $(TPLS_EXE) $(PETSC_LIB) $(NETCDF_LIB)

create_initial_conditions.linux : $(INITIAL_CONDITIONS_SRC) 
	mpif90 -O3 -std=f2003 -pedantic $^ $(PETSC_FC_INCLUDES) $(NETCDF_INCLUDES) -o $(INITIAL_CONDITIONS_EXE) $(PETSC_LIB) $(NETCDF_LIB)
#	mpif90 -O3 -Wall -std=f2003 -pedantic $^ -o $(INITIAL_CONDITIONS_EXE)

grids.linux : $(GRIDS_SRC)
	mpif90 -O3 -Wall -std=f2003 -pedantic $^ -o $(GRIDS_EXE)

run_tests.linux : $(TEST_SRC) $(TESTS)
	mpif90 -O3 -Wall -std=f2003 -pedantic $(FRUIT_SRC) $^ $(PETSC_FC_INCLUDES) $(NETCDF_INCLUDES) -o $(TEST_EXE) $(PETSC_LIB) $(NETCDF_LIB)

apidoc : Doxyfile src/*.f90 src/*.F90
	doxygen

clean :
	rm -f $(GRIDS_EXE)
	rm -f $(TPLS_EXE)
	rm -f $(INITIAL_CONDITIONS_EXE)
	rm -f $(TEST_EXE)
	rm -f result*.xml
	rm -f *.mod 
	rm -f *.o
	rm -f *~
	rm -rf html
	rm -rf latex
	rm -f *.dat
