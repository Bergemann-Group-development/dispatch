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
$(EXEC): read_external_mod.o mpi_debug_mod.o setup_mod.o rt_solvers_mod.o

# rt_solvers_mod.f90 is a Fortran submodule of rt_patch_mod (SUBMODULE (rt_patch_mod)
# rt_solvers_smod): it has no "use" statement for deps.pl to see, so nothing pulls its
# object into the link automatically (hence the explicit prerequisite above), and its
# compile order relative to its ancestor needs stating explicitly too.
rt_solvers_mod.o: rt_patch_mod.o

# The commit hash and branch reach the compiler as -DHASH and -DBRANCH. io_mod.f90 prints
# them in the start banner and job_record_mod.f90 records the hash, but make cannot see a
# changed -D value, so an incremental build in a checkout that moved to another commit
# would keep the old stamp. The stamp file changes with them and both objects depend on it.
HASH_STAMP := $(BUILD)/hash.stamp
$(HASH_STAMP): FORCE
	@mkdir -p $(BUILD)
	@echo '$(HASH) $(BRANCH)' | cmp -s - $@ || echo '$(HASH) $(BRANCH)' > $@
FORCE:
io_mod.o job_record_mod.o: $(HASH_STAMP)

# Fortran INCLUDE files: the compiler's .d files only feed the link list here, so
# their includers are named explicitly
rt_solvers_mod.o bftrans_mod.o level_mod.o: HOP/exp_neg.inc
# The standalone HOP package is built in its own directory with -DHOP_STANDALONE, so its
# objects never mix with the framework-enabled ones above.  Its default target builds
# what the Python layer runs: libhop.so (tumult.hop, through the committed libhop.h) and
# the hop_prune.x line-list pruner.  Building both here keeps a plain `make` sufficient
# for every entry point -- a build that produced only one of them fails at run time in
# the other.  The package takes this build's compiler and option group.
default: hop_standalone
hop_standalone:
	@$(MAKE) --no-print-directory -C HOP/standalone COMPILER=$(COMPILER) OPTS=$(OPTS)
clean:: clean_standalone
clean_standalone:
	@$(MAKE) --no-print-directory -C HOP/standalone clean
.PHONY: hop_standalone clean_standalone

# Input data that is tracked compressed, unpacked on first build to the path HOP reads
default: input/Anstee-1802-2014.dat
input/%.dat: input/%.dat.gz
	gunzip -c $< > $@
