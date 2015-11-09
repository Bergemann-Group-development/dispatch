# Default target
default: run

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

# Dependencies

run: advection_test.x advection_mod.x advection_type.x
	./advection_test.x < input.nml | tee advection_test.log
	./advection_mod.x  < input.nml | tee advection_mod.log
	./advection_type.x < input.nml | tee advection_type.log
