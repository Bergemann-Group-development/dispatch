VPATH := $(VPATH):$(TOP)/microphysics/eos/$(EOS_TYPE)/:$(TOP)/data/eos/$(EOS_READER)/

ifeq ($(EOS_TYPE),)
else
  sinclude $(TOP)/microphysics/eos/$(EOS_TYPE)/Makefile
endif

ifeq ($(EOS_READER),)
else
  sinclude $(TOP)/data/eos/$(EOS_READER)/Makefile
endif


