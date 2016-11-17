# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Common configuration
EOS_READER = simple
EOS_TYPE   = E
TOP        = ../../..
EXEC       = test.x
include $(TOP)/config/Makefile

test.x: io_mod.o trace_mod.o eos_reader_mod.o eos_mod.o
