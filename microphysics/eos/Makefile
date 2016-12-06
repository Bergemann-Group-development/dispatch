
ifeq ($(EOS),)
  EOS = ideal
endif
VPATH := $(VPATH):$(TOP)/microphysics/eos/$(EOS)/
sinclude          $(TOP)/microphysics/eos/$(EOS)/Makefile

ifeq ($(EOS_DATA),)
else
  VPATH := $(VPATH):$(TOP)/data/eos/$(EOS_DATA)/
  sinclude          $(TOP)/data/eos/$(EOS_DATA)/Makefile
endif

