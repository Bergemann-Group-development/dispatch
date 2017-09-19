
# Option groups
#OPTS = debug
OPTS = full_debug
#OPTS = optimized

# Path to ZEUScore shared library.
# Don't forget to update the LD_LIBRARY_PATH environment variable.
LIBS+=-L/home/jramsey/repos/ZEUScore -lzeuscore
PREC = DOUBLE

# Targets and methods
UPD = adams_bashforth
SOLVER = zeus3d
EXEC = dispatch.x
TOP = ../..
sinclude $(TOP)/config/Makefile

# Dependencies
$(EXEC) : setup_mod.o
setup_mod.o:		task_list_mod.o io_mod.o trace_mod.o experiment_mod.o components_mod.o scaling_mod.o
scaling_mod.o:		io_mod.o trace_mod.o
pboundary_mod.o: 	io_mod.o mesh_mod.o ic_mod.o mhd_mod.o kinds_mod.o
experiment_mod.o:	io_mod.o trace_mod.o timestep_mod.o mesh_mod.o ic_mod.o mhd_mod.o pboundary_mod.o
components_mod.o:	io_mod.o trace_mod.o task_mod.o mhd_mod.o task_list_mod.o experiment_mod.o bits_mod.o immersed_boundary_mod.o
ic_mod.o:			io_mod.o trace_mod.o mesh_mod.o mhd_mod.o kinds_mod.o
immersed_boundary_mod.o: io_mod.o patch_mod.o mhd_mod.o
