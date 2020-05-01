
# Option groups
#OPTS = full_debug
#OPTS = debug
OPTS =
#OPTS = optimized

# Default target and methods
UPD             = adams_bashforth
SOLVER          = muscl
#SOLVER          = ramses/hydro
#SOLVER          = ramses/mhd
TOP             = ../../..

# Tomida & Hori, with opacity
#EOS             = data/eos/Tomida+Hori_2016/e
#OPACITY         = Tomida+Hori_2016

sinclude $(TOP)/config/Makefile

#VPATH := $(VPATH):$(TOP)/microphysics/opacities/stub

# List the .o files the main program depends on after the colon:
$(EXEC)                 : setup_mod.o cartesian_mod.o dispatcher_mod.o lagrange_mod.o restart_mod.o \
                          task_block_mod.o
experiment_mod.o        : io_mod.o io_unit_mod.o trace_mod.o omp_mod.o bits_mod.o solver_mod.o initial_mod.o \
                          timer_mod.o link_mod.o
scaling_mod.o           : io_mod.o units_mod.o
extras_mod.o            : forces_mod.o

# .o files the mockup2 program depends on:
mockup2.o               : muscl.o
muscl.o                 : riemann.o
riemann.o               :
riemann2.o              :
muscl2.o                : riemann2.o
mockup3.o               : muscl2.o
mockup4.o               :

mockup.x: mockup.o
	$(FC) $(LDFLAGS) $(BUILD)/mockup.o -o mockup.x

mockup2.x: mockup2.o
	$(FC) $(LDFLAGS) $(BUILD)/{mockup2,muscl,riemann}.o -o mockup2.x

mockup3.x: mockup3.o
	$(FC) $(LDFLAGS) $(BUILD)/{mockup3,muscl2,riemann2}.o -o mockup3.x

mockup4.x: mockup4.o
	$(FC) $(LDFLAGS) $(BUILD)/mockup4.o -o mockup4.x
