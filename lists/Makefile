# Append to VPATH
VPATH := $(VPATH):$(TOP)/lists

# Local object files to include (hardwired later)
ifeq ($(EXEC),dispatch.x)
	OBJ += \
	list_mod.o \
	mesg_list_mod.o \
	task_list_mod.o
endif

# Dependencies
list_mod.o        : task_mod.o particle_mod.o planet_mod.o star_mod.o pebble_mod.o
task_list_mod.o   : list_mod.o tasks_mod.o
mesg_list_mod.o	  : task_list_mod.o
nbor_list_mod.o   : simple_list_mod.o task_list_mod.o task_mod.o
