
# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = pic
STAGGER         = stagger_2nd.o
MESH            = staggered
EXEC            = dispatch.x
TOP             = ../..

# PIC component selection (needs to be defined before sinclude)
PIC_STAGGER         = 2nd
INTERPOLATION   = tsc

# main Makefile -- sincludes other Makefiles
sinclude $(TOP)/config/Makefile

# To include the PPcode, PPCODE should point to the PPcode top directory
PPCODE          = $(TOP)/solvers/pic/PPcode
#PPCODE          = $(TOP)/../ppcode

 # Search paths for pic/PP module
  VPATH:=\
  $(VPATH):\
  $(PPCODE):\
  $(PPCODE)/Boundaries:\
  $(PPCODE)/Fields:\
  $(PPCODE)/Init:\
  $(PPCODE)/Interpolation:\
  $(PPCODE)/Io:\
  $(PPCODE)/Mpi:\
  $(PPCODE)/Utilities:\
  $(PPCODE)/Particles:\
  #$(PPCODE)/Experiments/DISPATCH

  # Object file list (automatic generation is not useful here)
  PICOBJ += params.o header.o debug.o grid.o units.o species.o dump.o stat.o textio.o \
  omp.o $(MPI).o distributions.o particleoperations.o io.o random.o fileoperations.o errors.o trace.o \
  pic_stagger_$(PIC_STAGGER).o mesh_sources.o maxwell_sources.o math.o laplace.o energy_filter.o \
  maxwell_solver_$(PIC_STAGGER).o load_balance.o sort.o init_ppcode.o bitoperations.o timer.o set_dt.o \
  move_particles.o celloperations.o init_pic_experiment.o interpolation_$(INTERPOLATION)_3d.o periodic.o \
  init_boundaries.o tests.o picmhd.o init_ppcode.o pic_mod.o

# List the .o files the main program depends on after the colon:
$(EXEC): cartesian_mod.o $(PICOBJ)

# Dependencies
experiment_mod.o        : solver_mod.o pic_mod.o trace_mod.o io_mod.o timestep_mod.o download_mod.o link_mod.o \
			  bits_mod.o params.o units.o species.o debug.o particleoperations.o textio.o \
                          initial_mod.o driver_mod.o boundary_mod.o
dispatch.o              : io_mod.o mpi_mod.o cartesian_mod.o
initial_mod.o           : io_mod.o trace_mod.o mesh_mod.o textio.o params.o grid.o units.o species.o debug.o particleoperations.o index_mod.o
boundary_mod.o          : io_mod.o trace_mod.o mesh_mod.o
bndry_mod.o             : io_mod.o trace_mod.o mesh_mod.o
driver_mod.o            : io_mod.o mesh_mod.o scalar_mod.o mhd_mod.o scaling_mod.o $(STAGGER) bndry_mod.o

