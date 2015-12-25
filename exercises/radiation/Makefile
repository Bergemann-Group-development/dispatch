# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target
PDE = radiation
EXEC= experiment.x

# Common configuration
AUTO = yes
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
VPATH:=$(VPATH):solvers/

# Dependencies
$(PDE)_mod.o      : omp_timer_mod.o io_mod.o trace_mod.o patches_mod.o scalar_mod.o cgs_mod.o scaling_mod.o microphysics_mod.o rt_mod.o 
patches_mod.o     : trace_mod.o io_mod.o bits_mod.o
rt_mod.o          : io_mod.o trace_mod.o patches_mod.o feautrier.o integral.o
scaling_mod.o     : io_mod.o cgs_mod.o
planck_mod.o      : io_mod.o cgs_mod.o scaling_mod.o
microphysics_mod.o: scaling_mod.o cgs_mod.o io_mod.o planck_mod.o

# Depeding on the AUTO setting, use auto_dependencies or not
ifeq ($(AUTO),yes)
  OBJ = $(PDE)_mod.o
else
  OBJ = \
    cgs_mod.o mpi_mod.o omp_mod.o omp_timer_mod.o io_mod.o trace_mod.o \
    $(PDE)_mod.o patches_mod.o scalar_mod.o rt_mod.o microphysics_mod.o scaling_mod.o \
    integral.o feautrier.o feautrier8.o planck_mod.o
endif

# Executable name and dependency
$(EXEC): $(OBJ)
