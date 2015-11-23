# This Makefile is a template for experiments, and may only need to be changed in the Dependencies section

# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized
#OPTS =

# Problem specification
MESH = centered
PROBLEM = advection
TIMESTEP = adams_bashforth

# Local object files
OBJ += envelope_mod.o planets_mod.o pebbles_mod.o atmosphere_mod.o experiment_mod.o component_mod.o

# Common configuration
EXEC = dispatch.x
TOP = ../../..
include $(TOP)/config/Makefile

# Dependencies
atmosphere_mod.o:				planet_mod.o patch_mod.o task_list_mod.o
planets_mod.o:          planet_mod.o task_list_mod.o io_mod.o trace_mod.o atmosphere_mod.o
pebble_mod.o:           particle_mod.o task_mod.o io_mod.o trace_mod.o
pebbles_mod.o:          pebble_mod.o task_list_mod.o task_mod.o io_mod.o trace_mod.o
envelope_mod.o:         mhd_mod.o task_list_mod.o task_mod.o io_mod.o trace_mod.o
experiment_mod.o:       task_list_mod.o envelope_mod.o star_mod.o planets_mod.o pebbles_mod.o trace_mod.o
