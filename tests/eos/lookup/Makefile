# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Common configuration
SOLVER = stagger2/e
EOS = eos_table_mod.o
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

test.x: mpi_mod.o $(EOS)

perflog: test.x
	mkdir -p performance
	uname -a | tee -a performance/$(HOST).log
	./test.x | tee -a performance/$(HOST).log
