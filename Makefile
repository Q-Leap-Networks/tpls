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

BASE_DIR=/apps/local/opse

INCLUDES=-I$(BASE_DIR)/include
# Uncomment for petsc 3.7.x
# INCLUDES=$(INCLUDES) -I$(BASE_DIR)/include/petsc
LIBS= -lpetsc -lnetcdff -lnetcdf -lhdf5_hl -lhdf5

all: $(TPLS_EXE) $(INITIAL_CONDITIONS_EXE) $(GRIDS_EXE) $(TEST_EXE)

twophase.x : $(SRC)
	mpif90.openmpi-gcc-2.0 -O3 -fopenmp -std=f2003 $^ -o $(TPLS_EXE) \
	  $(INCLUDES) $(LIBS)

create_initial_conditions : $(INITIAL_CONDITIONS_SRC) 
	mpif90.openmpi-gcc-2.0 -O3 -std=f2003 -pedantic $^ \
	  -o $(INITIAL_CONDITIONS_EXE) $(INCLUDES) $(LIBS)

grids : $(GRIDS_SRC)
	mpif90.openmpi-gcc-2.0 -O3 -Wall -std=f2003 -pedantic $^ -o $(GRIDS_EXE)

run_tests : $(TEST_SRC) $(TESTS)
	mpif90.openmpi-gcc-2.0 -O3 -Wall -std=f2003 -pedantic $(FRUIT_SRC) $^ \
	   -o $(TEST_EXE) $(INCLUDES) $(LIBS)

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
