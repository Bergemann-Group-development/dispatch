VPATH := $(VPATH):$(TOP)/data/opacities/$(OPACITY)/

opacity_mod.o           : scaling_mod.o opacity
opacity_reader_mod.o    : io_mod.o trace_mod.o opacity_mod.o
