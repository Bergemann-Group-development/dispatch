#############################################################################
# Default target and shell
target: default
SHELL = /bin/bash

#############################################################################
# Defalt compiler and directives
OMP=yes
COMPILER = gfortran
MPI = mpi
MOD = mod
AUTO = yes

#############################################################################
# Check for options
CONFIG = $(TOP)/config
sinclude options.mkf
sinclude $(CONFIG)/host/$(HOST)/Makefile
sinclude $(CONFIG)/compiler/$(COMPILER)/Makefile
sinclude $(CONFIG)/compiler/$(COMPILER)/$(HOST).mkf
sinclude $(CONFIG)/compiler/$(HOST)/$(OPTS).mkf
sinclude $(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf
sinclude $(CONFIG)/compiler/$(COMPILER)/$(OPTS)_$(HOST).mkf

#############################################################################
# Executable
#EXEC=dispatch.x
ifeq ($(OPTS),phi)
  EXEC := $(EXEC).mic
endif
default: $(EXEC)

#############################################################################
# Sources directories are searched in this exact order (VPATH inside Makefile)
VPATH = $(TOP)
sinclude $(TOP)/par/omp/Makefile
ifeq ($(MPI),mpi)
  sinclude $(TOP)/par/mpi/Makefile
else
  sinclude $(TOP)/par/nompi/Makefile
endif
sinclude $(TOP)/io/Makefile
#sinclude $(TOP)/operators/Makefile
sinclude $(TOP)/solvers/$(SOLVER)/Makefile
sinclude $(TOP)/solvers/rt/$(RT_SOLVER)/Makefile
sinclude $(TOP)/components/Makefile
sinclude $(TOP)/lists/Makefile
sinclude $(TOP)/tasks/Makefile
sinclude $(TOP)/meshes/$(MESH)/Makefile
sinclude $(TOP)/microphysics/Makefile
sinclude $(TOP)/patches/$(MESH)/Makefile
sinclude $(TOP)/timestep/$(UPD)/Makefile
sinclude $(TOP)/data/opacities/$(OPACITIES)/Makefile
sinclude $(TOP)/particles/Makefile
sinclude $(TOP)/states/Makefile

#############################################################################
# Optional patches and extra dependencies
ifdef PATCH
  VPATH := $(PATCH):$(VPATH)
  sinclude $(PATCH)/Makefile
endif
VPATH := ./:$(VPATH)
sinclude $(CONFIG)/Makefile.dep

#############################################################################
# Macro to hide compiler output
E = @

ifeq ($(AUTO),yes)
#############################################################################
# Infrastructure to automatically determine dependencies
DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
ifeq ($(COMPILER),gfortran)
  DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
endif
ifeq ($(COMPILER),ifort)
  DEPFLAGS = -gen-depformat=make -gen-dep=$(DEPDIR)/$*.Td
endif
POSTCOMPILE= mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d
FILES = $(subst .d,.o,$(notdir $(wildcard $(DEPDIR)/*.d)))

# Remove default rule for .f files.
%.o: %.f
%.o: %.mod
#############################################################################
# Pattern rules to get .o files from .f90 files
%.o: %.f90 $(DEPDIR)/%.d
	@if [ "x$(E)" == "x@" ] ; then echo " [ FC ]" $@ ; fi
	$(E)$(F90) $(FFLAGS) $(DEPFLAGS) -c $< -o $@
	$(E)$(POSTCOMPILE)
	$(eval OBJS+=$@)
%.o: %.f $(DEPDIR)/%.d
	@if [ "x$(E)" == "x@" ] ; then echo " [ FC ]" $@ ; fi
	$(E)$(F90) $(FFLAGS) $(DEPFLAGS) -c $< -o $@
	$(E)$(POSTCOMPILE)
	$(eval OBJS+=$@)
%.x: %.f90
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $@ ; fi
	$(eval OBJS+=$(subst .d,.o,$(notdir $(wildcard $(DEPDIR)/*.d))))
	$(E)$(LD) $(LDFLAGS) $< $(sort $(OBJS)) -o $@ $(LIBS)
%.x.mic: %.f90
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $@ ; fi
	$(eval OBJS+=$(subst .d,.o,$(notdir $(wildcard $(DEPDIR)/*.d))))
	$(E)$(LD) $(LDFLAGS) $< $(sort $(OBJS)) -o $@ $(LIBS)
%.w: %.f90
	@echo $<
%.w: %.f
	@echo $<
$(DEPDIR)/%.d: ;

# Prevent deletion of dependency files after successful compilation
.PRECIOUS: $(DEPDIR)/%.d

.PHONY: showopt clean cleanlog cleandata clean-deps cleaner cleanest info
else
#############################################################################
# Rule to get .o files from .f90 files
.SUFFIXES:
.SUFFIXES: .x.mic .x .o .f90 .f .w
.f90.x.mic: $(OBJ)
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $@ ; fi
	$(E)$(LD) $(LDFLAGS) $(OBJ) $< -o $@ $(LIBS)
.f90.x: $(OBJ)
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $@ ; fi
	$(E)$(LD) $(LDFLAGS) $(OBJ) $< -o $@ $(LIBS)
.f90.o:
	@if [ "x$(E)" == "x@" ] ; then echo " [ FC ]" $@ ; fi
	$(E)$(F90) $(FFLAGS) -c $< -o $@
.f.o:
	@if [ "x$(E)" == "x@" ] ; then echo " [ FC ]" $@ ; fi
	$(E)$(F90) $(FFLAGS) -c $< -o $@
.f90.w:
	@echo $<
.f.w:
	@echo $<
$(EXEC): $(OBJ)
endif

showopt:
	@echo "FC = $(F90) $(FFLAGS)"
	@echo "LD = $(LD) $(LDFLAGS)"

#############################################################################
## Clean up
clean: clean-deps
	rm -f *.o *.$(MOD) *.stackdump
cleanlog:
	rm -f *.log log*
cleandata:
	rm -f *.dat fort.*
cleandirs:
	rm -rf data .d
clean-deps:
	rm -f .d/*.*d
cleaner: clean cleandata
	rm -f *.x *.bak *.optrpt *.exe *.pyc
cleanest: cleaner cleandata cleanlog clean-deps cleandirs

#############################################################################
# Generate list of source files
files: 
	@echo "$(FILES)"
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
	@echo "$B  LIBS        $N= $(LIBS)"
	@echo "$B  EXEC        $N= $(EXEC)"

#############################################################################
# Check where options are set
PATTERN = OPT
showinclude:
	@echo ""; echo "searching for PATTERN=$(PATTERN), use 'make showinclude PATTERN=xxx' to search for xxx"; echo ""
	@-for f in \
Makefile \
$(CONFIG)host/$(HOST)/Makefile \
$(CONFIG)/compiler/$(COMPILER)/Makefile \
$(CONFIG)/compiler/$(COMPILER)/$(HOST).mkf \
$(CONFIG)/compiler/$(HOST)/$(OPTS).mkf \
$(CONFIG)/compiler/$(COMPILER)/$(OPTS).mkf \
$(CONFIG)/compiler/$(COMPILER)/$(OPTS)_$(HOST).mkf \
$(PATCH)/Makefile \
	; do if [ -e $$f ]; then echo "------ $$f ------"; grep -w \\$(PATTERN) $$f; fi; done
