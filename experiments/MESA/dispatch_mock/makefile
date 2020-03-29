#include $(MESA_DIR)/star/work_standard_makefile

#OPTS = full_debug
OPTS = debug
#OPTS =
#OPTS = optimized

TOP            = ../../..
EXEC           = run
SOLVER         = ramses/hydro

include $(TOP)/config/Makefile
	
all: test.x
	cd ../mesa_template; pwd; ./mk; cd -;

run: test.x
test.x: mesa_mod.o
mesa_mod.o: io_unit_mod.o mpi_mod.o