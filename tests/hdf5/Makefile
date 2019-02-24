# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# The location of HDF5 may vary.  At steno the default
# location is given in config/host/steno/Makefile, as:
#HDF5 = /software/astro/hdf5-1.10.4

VERBOSE ?= 1
default: test.x
	./test.x $(VERBOSE)
	ls -l hdf5.dat
ifneq ($(VERBOSE),0)
	$(HDF5)/bin/h5ls.exe -lr hdf5.dat
endif

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

test.x: h5_mod.o
