
# Local object files to include
OBJ += \
task_mod.o \
list_mod.o \
patch_mod.o \
timer_mod.o

list_mod.o: task_mod.o

patch_mod.o:	mesh_mod.o task_mod.o
hd_mod.o:	patch_mod.o
mhd_mod.o:	hd_mod.o
