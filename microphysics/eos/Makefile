
ifeq ($(EOS),)
  EOS = ideal
endif
sinclude $(TOP)/microphysics/eos/$(EOS)/Makefile

ifeq ($(EOS_DATA),)
else
  sinclude $(TOP)/data/eos/$(EOS_DATA)/Makefile
endif

VPATH := $(VPATH):$(TOP)/microphysics/eos/$(EOS)/:$(TOP)/data/eos/$(EOS_DATA)/
