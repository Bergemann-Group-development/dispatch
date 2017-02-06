# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Common configuration
UPD=adams_bashforth
SOLVER=zeus3d
EXEC = test.x
TOP = ../..
include $(TOP)/config/Makefile

LDFLAGS += -cpp
LIBS+=-L/home/jramsey/repos/ZEUScore -lzeuscore

test.x: mesh_mod.o io_mod.o mhd_mod.o patch_mod.o
