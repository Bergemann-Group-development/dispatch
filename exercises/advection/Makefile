# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target
default: run

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

# Dependencies
run: classical.x module_based.x object_oriented.x
	./classical.x       < input.nml | tee classical.log
	./module_based.x    < input.nml | tee module_based.log
	./object_oriented.x < input.nml | tee object_oriented.log
