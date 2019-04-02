#############################################################################
# Default target and shell
SHELL = /bin/bash

#############################################################################
# Defalt compiler and directives
COMPILER = gfortran
OMP  = yes
MPI  = mpi
MOD  = mod
DEFS = -D$(PREC) -DDISPATCH
SOLVER  ?= VOID
PREC    ?= SINGLE
RT_PREC ?= RT_SINGLE
BUILD   ?= build

#############################################################################
# Check for options
CONFIG = $(TOP)/config
# sinclude files that are likely to change $(COMPILER)
sinclude options.mkf
sinclude $(CONFIG)/host/$(HOST)/Makefile
# check that the $(OPTS) file is present (typo trap)
ifeq ($(OPTS),)
  target: default
else
  OPTS_COMP = $(shell ls $(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf) 
  OPTS_FILE = $(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf
  ifeq ("$(OPTS_COMP)","$(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf ")
    target: default
  else
    target: error
  endif
endif
# sinclude $(OPTS) and $(HOST) dependent config files
sinclude $(CONFIG)/compiler/$(COMPILER)/Makefile
sinclude $(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf
sinclude $(CONFIG)/compiler/$(COMPILER)/$(HOST)/Makefile
sinclude $(CONFIG)/compiler/$(COMPILER)/$(HOST)/$(OPTS).mkf

#############################################################################
# Executable
#EXEC=dispatch.x
ifeq ($(OPTS),phi)
  EXEC := $(EXEC).mic
endif
default: build $(EXEC)

#############################################################################
# Sources directories are searched in this exact order (VPATH inside Makefile)
VPATH = $(TOP)
sinclude $(TOP)/omp/Makefile
sinclude $(TOP)/mpi/Makefile
ifeq ($(MPI),mpi)
  DEFS += -DMPI
endif
sinclude $(TOP)/io/Makefile
#sinclude $(TOP)/operators/Makefile
sinclude $(TOP)/solvers/Makefile
sinclude $(TOP)/extras/Makefile
sinclude $(TOP)/components/Makefile
sinclude $(TOP)/dispatchers/Makefile
sinclude $(TOP)/lists/Makefile
sinclude $(TOP)/interpolation/Makefile
sinclude $(TOP)/tasks/Makefile
sinclude $(TOP)/meshes/Makefile
sinclude $(TOP)/microphysics/Makefile
sinclude $(TOP)/timestep/$(UPD)/Makefile
sinclude $(TOP)/particles/Makefile
sinclude $(TOP)/states/Makefile
sinclude $(TOP)/boundaries/Makefile

#############################################################################
# HDF5 option
ifneq ($(HDF5),)
  INCL += -I$(HDF5)/include
  LIBS += -L$(HDF5)/lib $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a $(HDF5)/lib/libhdf5hl_fortran.a $(HDF5)/lib/libhdf5_hl.a -lz -ldl -lm
endif
FFLAGS += $(INCL)
LDFLAGS += $(INCL)

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

#############################################################################
# Infrastructure to automatically determine dependencies
DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)

# Remove default rule for .f files.
%.o: %.f
%.o: %.mod
#############################################################################
# Pattern rules
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
.PHONY: showopt clean cleanlog cleandata clean-deps cleaner cleanest info bitvalidate

$(BUILD):
	@mkdir -p $(BUILD)

showopt:
	@echo "FC = $(F90) $(FFLAGS) $(DEFS)"
	@echo "LD = $(LD) $(LDFLAGS)"

#############################################################################
## Clean up
CLEANDIRS = {.,$(BUILD)}

clean:: clean-deps
	rm -f $(CLEANDIRS)/*.o $(CLEANDIRS)/*.$(MOD) *$(CLEANDIRS)/.stackdump
cleanlog:
	rm -f *.log log*
cleandata:
	rm -f *.dat fort.*
cleandirs:
	rm -rf data .d $(BUILD)
clean-deps:
	rm -f .d/*.*d
cleaner: clean cleandata
	rm -f *.x *.bak *.optrpt *.exe *.pyc
cleanest: cleaner cleandata cleanlog clean-deps cleandirs

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
# Output info about compiler options
info:
	@echo "Selected Makeflags and options:"
	@echo "$B  HOST        $N= $(HOST)"
	@echo "$B  COMPILER    $N= $(COMPILER)"
	@echo "$B  MPI         $N= $(MPI)"
	@echo "$B  PATCH       $N= $(PATCH)"
	@echo "$B  VPATH       $N= $(VPATH)"
	@echo "$B  SOLVER      $N= $(SOLVER)"
	@echo "$B  OPACITY     $N= $(OPACITY)"
	@echo "$B  EOS         $N= $(EOS)"
	@echo "$B  PROG        $N= $(PROG).o"
	@echo "$B  OPTS        $N= $(OPTS)"
	@echo "$B  OPT         $N= $(OPT)"
	@echo "$B  HW          $N= $(HW)"
	@echo "$B  DEBUG       $N= $(DEBUG)"
	@echo "$B  PAR         $N= $(PAR)"
	@echo "$B  F90         $N= $(F90)"
	@echo "$B  LD          $N= $(LD)"
	@echo "$B  FFLAGS      $N= $(FFLAGS)"
	@echo "$B  LDFLAGS     $N= $(LDFLAGS)"
	@echo "$B  DEFS        $N= $(DEFS)"
	@echo "$B  LIBS        $N= $(LIBS)"
	@echo "$B  EXEC        $N= $(EXEC)"
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
