# Option groups
#OPTS = full_debug
#OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = optimized
OPTS = hoppe# use this when compiling for different hosts at once.

# Default target and paths
EXEC            = dispatch.x
TOP             = ../..

# Default methods
ifndef mpiifort
  MPI=
endif

ifeq ($(origin COMPILER), undefined)
  ifneq ($(shell command -v gfortran 2>/dev/null),)
    COMPILER = gfortran
  else ifneq ($(shell command -v ifort 2>/dev/null),)
    COMPILER = ifort
  else
    $(error "No Fortran compiler found. Install gfortran or ifort, or set COMPILER explicitly.")
  endif
endif

SOLVER          = None
EOS             = microphysics/None
STAGGER         = 

sinclude $(TOP)/config/Makefile
VPATH := ./SRC:./HOP:$(VPATH) # Search SRC and HOP first

# Automatically (re)generate the SRC/HOP module dependency graph by scanning
# the source text for "module"/"use" statements (config/deps.pl). 
EXTERNAL_DEP_DIRS := io mpi omp tasks lists meshes states microphysics solvers/rt
EXTERNAL_SOURCES  := $(wildcard $(addsuffix /*.f90,$(addprefix $(TOP)/,$(EXTERNAL_DEP_DIRS))))
LOCAL_DEPFILE := $(BUILD)/Makefile.dep
LOCAL_SOURCES := $(wildcard SRC/*.f90 HOP/*.f90 SRC/*.f HOP/*.f)
LOCAL_OBJECTS := $(notdir $(patsubst %.f90,%.o,$(patsubst %.f,%.o,$(LOCAL_SOURCES))))
# Keep any generated line that mentions a local object on EITHER side of the
# colon: local-depends-on-external (e.g. mpi_debug_mod.o uses mpi_mod.o) as
# well as external-depends-on-local (e.g. core tasks/gpatch_mod.f90 uses this
# experiment's own SRC/boundaries_mod.f90 -- a dependency the shared
# tasks/Makefile has no way to know about, since it's experiment-specific).
$(LOCAL_DEPFILE): $(LOCAL_SOURCES) $(EXTERNAL_SOURCES) $(TOP)/config/deps.pl
	@ls $(LOCAL_SOURCES) $(EXTERNAL_SOURCES) > $(BUILD)/files.txt
	@for o in $(LOCAL_OBJECTS); do echo "$$o"; done > $(BUILD)/local_objects.txt
	@perl $(TOP)/config/deps.pl $(BUILD)/files.txt | grep -Fwf $(BUILD)/local_objects.txt > $(LOCAL_DEPFILE) || true
sinclude $(LOCAL_DEPFILE)

ifeq ($(COMPILER), ifort)
  OPT += -diag-disable=10448
endif


CLEANDIRS = {./SRC,./HOP,$(BUILD)}

# Symbolication for sampling profilers (sample/dtrace/Instruments). -g adds DWARF
# only and does NOT change -Ofast codegen; frame pointers help stack unwinding.
# OPT += -g -fno-omit-frame-pointer

# List the .o files the main program depends on after the colon:
$(EXEC): read_external_mod.o mpi_debug_mod.o setup_mod.o