# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target
PDE = radiation
EXEC= experiment.x

# Common configuration
TOP = ../..
sinclude $(TOP)/config/Makefile
TOP = ../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../../../..
sinclude $(TOP)/config/Makefile

OBJ = \
  cgs_mod.o mpi_mod.o omp_mod.o omp_timer_mod.o io_mod.o trace_mod.o \
  $(PDE)_mod.o patches_mod.o scalar_mod.o rt_mod.o microphysics_mod.o scaling_mod.o

# Dependencies
$(PDE)_mod.o      : omp_timer_mod.o io_mod.o trace_mod.o patches_mod.o scalar_mod.o rt_mod.o
patches_mod.o     : trace_mod.o io_mod.o bits_mod.o
rt_mod.o          : io_mod.o trace_mod.o patches_mod.o mpi_mod.o microphysics_mod.o
scaling_mod.o     : io_mod.o cgs_mod.o
microphysics_mod.o: patches_mod.o scaling_mod.o

$(EXEC): $(OBJ)
