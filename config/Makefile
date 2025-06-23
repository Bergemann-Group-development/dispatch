#############################################################################
# Default target and shell
SHELL = /bin/bash

#############################################################################
# Generate macro with commit hash tag
BRANCH := $(shell git name-rev --name-only HEAD)
HASH := $(shell git log --pretty='format:%h' -n1 HEAD)
#FCFLAGS = $(FOPT) $(DEBUG) $(LIBS) -D'GIT_REV="$(REV)"' -cpp

#############################################################################
# Defalt compiler and directives
COMPILER ?= gfortran
OMP      ?= yes
MPI      ?= mpi
MOD      ?= mod
SOLVER   ?= VOID
OPTS     ?= devel
UPD      ?= adams_bashforth
PREC     ?= SINGLE
EXEC     ?= dispatch.x
OFFLOAD  ?= CPU

#############################################################################
# Default target
default: $(EXEC)

#############################################################################
# 1st part of Makefile macro info
info::
	@echo "Selected Makeflags and options:"
	@echo "$B  VERBOSE     $N= $(VERBOSE)"
	@echo "$B  HOST        $N= $(HOST)"
	@echo "$B  COMPILER    $N= $(COMPILER)"
	@echo "$B  MPI         $N= $(MPI)"
	@echo "$B  PATCH       $N= $(PATCH)"
	@echo "$B  VPATH       $N= $(VPATH)"
	@echo "$B  SOLVER      $N= $(SOLVER)"
	@echo "$B  OPACITY     $N= $(OPACITY)"
	@echo "$B  EOS         $N= $(EOS)"

#############################################################################
# Print hostinfo (relevant for hpc systems)
hostinfo:
	$(info $(HOSTINFO))
	@make info

#############################################################################
# sinclude files that are likely to change $(COMPILER)
CONFIG = $(TOP)/config
sinclude options.mkf
sinclude $(TOP)/options.mkf
sinclude $(CONFIG)/host/$(HOST)/Makefile
target: default

DEFS = -D$(PREC) -DDISPATCH -D'SOLVER="$(SOLVER)"' -D'BRANCH="$(BRANCH)"' -D'HASH="$(HASH)"' -D'OPTS="$(OPTS)"'

#############################################################################
# Check that the $(COMPILER) Makefile is present (typo trap)
ifneq ($(COMPILER),)
  COMPILER_MK = $(TOP)/config/compiler/$(COMPILER)/Makefile
  COMPILER_OS  = $(shell ls $(COMPILER_MK))
  ifneq ("$(COMPILER_OS)", "$(COMPILER_MK)")
    $(error $(COMPILER_MK) does not exist)
  endif
endif

#############################################################################
# Check that the $(SOLVER) Makefile is present (typo trap)
ifneq ($(SOLVER),VOID)
  SOLVER_MK = $(TOP)/solvers/$(SOLVER)/Makefile
  SOLVER_OS  = $(shell ls $(SOLVER_MK))
  ifneq ("$(SOLVER_OS)", "$(SOLVER_MK)")
    $(error $(SOLVER_MK) does not exist)
  endif
endif

#############################################################################
# Check that the $(OPTS) file is present (typo trap)
ifneq ($(OPTS),)
  OPTS_COMP = $(shell ls $(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf)
  OPTS_FILE = $(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf
  ifneq ("$(OPTS_COMP)", "$(OPTS_FILE)")
    $(error $(OPTS_FILE) does not exist)
  endif
endif

#############################################################################
# Optional separate build directories
SEPARATE_BUILDS ?= $(SEPARATE_BUILD)
ifeq ($(SEPARATE_HOSTS),)
	ifeq ($(SEPARATE_BUILDS),)
		BUILD    ?= build
	else
		BUILD    ?= build/$(SOLVER)/$(OPTS)
	endif
else
	ifeq ($(SEPARATE_BUILDS),)
		BUILD    ?= build/$(HOST)
	else
		BUILD    ?= build/$(HOST)/$(SOLVER)/$(OPTS)
	endif
endif
$(shell mkdir -p $(BUILD) >/dev/null)
$(info make HOST=$(HOST) COMPILER=$(COMPILER) OPTS=$(OPTS) BUILD=$(BUILD))

#############################################################################
# Infrastructure to automatically determine dependencies
DEPDIR := $(BUILD)/.d
$(shell mkdir -p $(DEPDIR) >/dev/null)

#############################################################################
# sinclude $(OPTS) and $(HOST) dependent config files
sinclude $(CONFIG)/compiler/$(COMPILER)/Makefile
sinclude $(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf
sinclude $(CONFIG)/compiler/$(COMPILER)/$(HOST)/Makefile
sinclude $(CONFIG)/compiler/$(COMPILER)/$(HOST)/$(OPTS).mkf

#############################################################################
# HDF5 option
ifneq ($(HDF5),)
  INCL += -I$(HDF5)/include
  LIBS += -L$(HDF5)/lib $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a $(HDF5)/lib/libhdf5hl_fortran.a $(HDF5)/lib/libhdf5_hl.a -lz -ldl -lm
  DEFS += -DHDF5
  HDF5_MOD = h5_mod.o
  ifneq ($(VERBOSE),)
    $(info $(TOP)/config/Makefile HDF5=$(HDF5) HDF5_MOD=$(HDF5_MOD))
  endif
endif
FFLAGS += $(INCL)
LDFLAGS += $(INCL)

#############################################################################
# Executable
#EXEC=dispatch.x
ifeq ($(OPTS),phi)
  EXEC := $(EXEC).mic
endif
default: $(BUILD) $(EXEC)

#############################################################################
# Sources directories are searched in this exact order (VPATH inside Makefile)
VPATH = $(TOP)
sinclude $(TOP)/Makefile
sinclude $(TOP)/omp/Makefile
sinclude $(TOP)/mpi/Makefile
ifeq ($(MPI),mpi)
  DEFS += -DMPI
endif
sinclude $(TOP)/io/Makefile
sinclude $(TOP)/solvers/Makefile
#sinclude $(TOP)/extras/Makefile
#sinclude $(TOP)/components/Makefile
sinclude $(TOP)/dispatchers/Makefile
sinclude $(TOP)/offload/Makefile
sinclude $(TOP)/lists/Makefile
sinclude $(TOP)/interpolation/Makefile
sinclude $(TOP)/tasks/Makefile
sinclude $(TOP)/meshes/Makefile
sinclude $(TOP)/microphysics/Makefile
#sinclude $(TOP)/timestep/$(UPD)/Makefile
#sinclude $(TOP)/particles/Makefile
sinclude $(TOP)/states/Makefile
#sinclude $(TOP)/boundaries/Makefile

#############################################################################
# Optional patches and extra dependencies
ifdef PATCH
  VPATH := $(TOP)/patch/$(PATCH):$(VPATH)
  sinclude $(TOP)/patch/$(PATCH)/Makefile
endif
VPATH := ./:$(VPATH):$(BUILD)
sinclude $(CONFIG)/Makefile.dep

#############################################################################
# Macro to hide compiler output
E = @

# Remove default rule for .f files.
%.o: %.f
%.o: %.mod
#############################################################################
# Pattern rules
CP = cp -p
COPY = $(BUILD)/$(EXEC)
%.o: %.f90 $(DEPDIR)/%.d
	@if [ "x$(E)" == "x@" ] ; then echo " [ FC ]" $@ ; fi
	$(E)$(F90) $(FFLAGS) $(DEFS) $(DEPFLAGS) -c $< -o $(BUILD)/$@
	$(E)$(POSTCOMPILE)
	$(eval OBJS+=$@)
%.o: %.f $(DEPDIR)/%.d
	@if [ "x$(E)" == "x@" ] ; then echo " [ FC ]" $@ ; fi
	$(E)$(F90) $(FFLAGS) $(DEFS) $(DEPFLAGS) -c $< -o $(BUILD)/$@
	$(E)$(POSTCOMPILE)
	$(eval OBJS+=$@)
%.x: %.f90
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $@ ; fi
	$(eval OBJS+=$(subst .d,.o,$(notdir $(wildcard $(DEPDIR)/*.d))))
	$(E)$(LD) $(LDFLAGS) $< $(sort $(patsubst %,$(BUILD)/%,$(OBJS))) -I$(BUILD) -o $@ $(LIBS)
	$(CP) $@ $(COPY)
%.x.mic: %.f90
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $@ ; fi
	$(eval OBJS+=$(subst .d,.o,$(notdir $(wildcard $(DEPDIR)/*.d))))
	$(E)$(LD) $(LDFLAGS) $< $(sort $(patsubst %,$(BUILD)/%,$(OBJS))) -I$(BUILD) -o $@ $(LIBS)
%.w: %.f90
	@echo $<
%.w: %.f
	@echo $<
$(DEPDIR)/%.d: ;

# Prevent deletion of dependency files after successful compilation
.PRECIOUS: $(DEPDIR)/%.d
.PHONY: showopt clean cleanlog cleandata cleandeps cleaner cleanest info bitvalidate pure

$(BUILD):
	@mkdir -p $(BUILD)

showopt:
	@echo "FC = $(F90) $(FFLAGS) $(DEFS)"
	@echo "LD = $(LD) $(LDFLAGS)"

#############################################################################
## Clean up
CLEANDIRS = {.,$(BUILD)}

clean:: cleandeps
	rm -f $(CLEANDIRS)/*.o $(CLEANDIRS)/*.$(MOD) *$(CLEANDIRS)/.stackdump
cleanlog:
	rm -f *.log log*
cleandata:
	rm -f *.dat
cleanfort:
	rm -f fort.*
cleandirs:
	rm -rf data $(DEPDIR) $(BUILD)
cleandeps:
	rm -f $(DEPDIR)/*.d
cleaner: clean cleanfort
	rm -f *.x *.bak *.optrpt *.exe *.pyc
	rm -rf build_*
cleanest: cleaner cleandata cleanlog cleandeps cleandirs
pure: cleanest
	rm -rf build*

#############################################################################
# Generate list of source files (requires successful build)
FILES = $(subst .d,.o,$(notdir $(wildcard $(DEPDIR)/*.d)))
SRC   = $(subst .d,.w,$(notdir $(wildcard $(DEPDIR)/*.d)))
files: $(FILES)
	@echo $(FILES)
src: $(SRC)
bitvalidate:
	csh $(TOP)/utilities/scripts/bitvalidate.csh $(SOLVER)
	OMP_NUM_THREADS=1 ./$(EXEC) $(SOLVER)/bitvalidate.nml

# Generate list of dependencies
deps: files
	@make files > files.txt
	@perl $(CONFIG)/deps.pl files.txt

#############################################################################
# update documentation via doxygen
doc: doxygen
doxygen: ../doc/doxyfile
	@doxygen ../doc/doxyfile >& /dev/null

#############################################################################
# 2nd part of Makefile macro info
info::
	@echo "$B  OPTS        $N= $(OPTS)"
	@echo "$B  OPT         $N= $(OPT)"
	@echo "$B  HW          $N= $(HW)"
	@echo "$B  DEBUG       $N= $(DEBUG)"
	@echo "$B  PAR         $N= $(PAR)"
	@echo "$B  F90         $N= $(F90)"
	@echo "$B  LD          $N= $(LD)"
	@echo "$B  FFLAGS      $N= $(FFLAGS)"
	@echo "$B  LDFLAGS     $N= $(LDFLAGS)"
	@echo "$B  BUILD       $N= $(BUILD)"
	@echo "$B  EXEC        $N= $(EXEC)"
	@echo "$B  DEFS        $N= $(DEFS)"
	@echo "$B  LIBS        $N= $(LIBS)"
	@echo "$B  BRANCH      $N= $(BRANCH)"
	@echo "$B  HASH        $N= $(HASH)"

compiler:
	@echo "$(COMPILER)"

error:
	@echo "$B  OPTS_FILE   $N= :$(OPTS_FILE): DOES NOT EXIST!"

#############################################################################
# Check where options are set
PATTERN = OPT
showinclude:
	@echo ""; echo "searching for PATTERN=$(PATTERN), use 'make showinclude PATTERN=xxx' to search for xxx"; echo ""
	@-for f in \
Makefile \
$(CONFIG)/host/$(HOST)/Makefile \
$(CONFIG)/compiler/$(COMPILER)/Makefile \
$(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf \
$(CONFIG)/compiler/$(COMPILER)/$(HOST)/Makefile \
$(CONFIG)/compiler/$(COMPILER)/$(HOST)/$(OPTS).mkf \
$(PATCH)/Makefile \
	; do if [ -e $$f ]; then echo "------ $$f ------"; grep -w \\$(PATTERN) $$f; fi; done
