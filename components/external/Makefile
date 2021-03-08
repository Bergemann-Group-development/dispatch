
ifneq ($(EXTERNAL),)
  VPATH := $(TOP)/components/external/$(EXTERNAL):$(VPATH)
  sinclude $(TOP)/components/external/$(EXTERNAL)/Makefile
endif

read_external_mod.o  : task_list_mod.o lagrange_mod.o experiment_mod.o scaling_mod.o
