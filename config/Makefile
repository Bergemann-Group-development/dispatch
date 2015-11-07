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
EXEC=dispatch.x
ifeq ($(OPTS),phi)
  EXEC=dispatch.x.mic
endif
default: $(EXEC)

#############################################################################
# Sources directories are searched in this exact order (VPATH inside Makefile)
VPATH = $(TOP)/modules
sinclude $(TOP)/modules/Makefile
sinclude $(TOP)/mpi/Makefile
#sinclude $(TOP)/$(MPI)/Makefile

#############################################################################
# Optional patches
ifdef PATCH
  VPATH := $(PATCH):$(VPATH)
  sinclude $(PATCH)/Makefile
endif
VPATH := ./:$(VPATH)
sinclude $(CONFIG)/Makefile.dep

#############################################################################
# Macro to hide compiler output
E = @

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
.f90.w:
	@echo $<
showopt:
	@echo "FC = $(F90) $(FFLAGS)"
	@echo "LD = $(LD) $(LDFLAGS)"

#############################################################################
# Recipe for linking
$(EXEC): $(OBJ)

#############################################################################
clean:
	rm -f *.o *.$(MOD) *.stackdump
cleandata:
	rm -f *.dat
cleaner: clean
	rm -f *.x *.bak *.optrpt
cleanest: cleaner cleandata

#############################################################################
# Generate list of source files
FILES = $(subst .o,.w,$(OBJ))
files: $(FILES)
	@echo xx >/dev/null

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
	@echo "$B  OBJ         $N= $(OBJ)"
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
