VPATH := $(VPATH):./

# Option groups
#OPTS = debug
#OPTS = full_debug
#OPTS = optimized
OPTS = 

OBJ += \
omp_mod.o \
mpi_basic.o \
mpi_timer.o \
io_mod.o \
trace_mod.o \
mesh_mod.o \
list_mod.o \
item_mod.o \
task_mod.o \
patch_mod.o \
mhd_mod.o \
particle_mod.o \
dust_mod.o \
star_mod.o \
planet_mod.o \
message_mod.o \
mesg_list_mod.o \
experiment_mod.o

default: dispatch.x

# Common configuration
TOP = ../..
include $(TOP)/config/Makefile

test.x: $(OBJ)
