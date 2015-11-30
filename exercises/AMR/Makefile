# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Default target
default: experiment.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

OBJ = \
  mpi_timer_mod.o upwind_mod.o central_mod.o derivatives_mod.o \
  io_mod.o patch_list_mod.o trace_mod.o omp_mod.o \
  scalar_mod.o vector_mod.o tensor_mod.o \
  hd_mod.o task_mod.o refine_mod.o timestep_mod.o

# Dependencies
task_mod.o              : trace_mod.o io_mod.o
hd_mod.o                : task_mod.o scalar_mod.o vector_mod.o tensor_mod.o
patch_list_mod.o        : io_mod.o timestep_mod.o refine_mod.o
meshes_mod.o            : trace_mod.o
derivatives_mod.o       : meshes_mod.o
timestep_mod.o          : io_mod.o hd_mod.o
refine_mod.o            : timestep_mod.o
experiment.o            : io_mod.o patch_list_mod.o

experiment.x: $(OBJ)
