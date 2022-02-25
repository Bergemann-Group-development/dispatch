# Option groups
#OPTS = full_debug
#OPTS = vscode
#OPTS = debug
#OPTS =
OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
#SOLVER         = stagger
SOLVER          = ramses/hlld_eos
SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd_eos/s
#SOLVER          = stagger2/bifrost
#STAGGER         = stagger_2nd.o
STAGGER         = stagger_2nd.o
#STAGGER         = stagger_mod.o
EXEC            = dispatch.x
#EOS             = data/eos/tomida/e
TOP             = ../../..
sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC): restart_mod.o cartesian_mod.o setup_mod.o dispatcher_mod.o initial_conditions_mod.o

# Dependencies
benchmark_mod.o          : trace_mod.o
cooling_mod.o            : io_mod.o kinds_mod.o scratch_mod.o benchmark_mod.o patch_mod.o
experiment_mod.o         : io_mod.o bits_mod.o solver_mod.o trace_mod.o interpolate_mod.o
initial_conditions_mod.o : io_mod.o trace_mod.o kinds_mod.o mesh_mod.o \
                           patch_mod.o index_mod.o units_mod.o random_mod.o \
													math_mod.o
extras_mod.o             : patch_mod.o scratch_mod.o initial_conditions_mod.o \
                           gravity_mod.o viscosity_mod.o conduction_mod.o \
													cooling_mod.o
boundaries_mod.o         : eos_mod.o io_mod.o trace_mod.o kinds_mod.o bits_mod.o \
                           mesh_mod.o index_mod.o standard_boundaries_mod.o \
                           patch_mod.o benchmark_mod.o scratch_mod.o

plot: 
	@python make_plot.py
