# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/$(SOLVER)

include $(TOP)/solvers/zeus3d/Makefile

PREC=DOUBLE

# **Don't forget to update the LD_LIBRARY_PATH environment variable when you run the code!**
ifeq ($(SOLVER),zeus3d/zeustw)
  LIBS+=-L/home/jramsey/repos/ZeusTWCore -lzeustwcore
  mhd_mod.o:      FFLAGS+=-I/home/jramsey/repos/ZeusTWCore
endif
