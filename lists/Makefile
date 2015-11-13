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
list_mod.o				: item_mod.o
task_list_mod.o 	: list_mod.o patch_mod.o
mesg_list_mod.o		: task_list_mod.o
