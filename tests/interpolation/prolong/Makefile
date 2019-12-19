# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

default: test.x

# Common configuration
TOP = ../../..
include $(TOP)/config/Makefile

test.x         : prolong_mod.o
prolong_mod.o  : io_mod.o
