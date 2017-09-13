# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

default: test.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

test.x: particles_mod.o

particle_history_mod.o 	: particle_mod.o
particles_mod.o 	: particle_history_mod.o
