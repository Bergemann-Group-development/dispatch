# Append to VPATH
ifeq ($(PREC),SINGLE)
  DEFS := $(DEFS) -DPREC=4
else
  DEFS := $(DEFS) -DPREC=8
endif

VPATH := $(VPATH):$(TOP)/solvers/stagger2/pic

include $(TOP)/solvers/stagger2/Makefile
#sinclude $(TOP)/solvers/pic/Makefile.dep
