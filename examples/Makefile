# Template makefile for example directories.   If they only require 
# compilation of a test.f90 (which should be the default name), then
# a copy of this Makefile will suffice.  Available make commands:

#     make                      # compile
#     make info                 # show Makefile macros
#     make clean                # remove *.o and *.mod files, plus some
#     make cleaner              # also remove *.x
#     make cleanest             # also remove *.dat

# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Optional extra object files
OBJ =
test.x: $(OBJ)

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

