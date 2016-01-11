# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Default target
EXEC= experiment.x

# Common configuration
AUTO = yes
TOP = ../..
sinclude $(TOP)/config/Makefile
TOP = ../../..
sinclude $(TOP)/config/Makefile
TOP = ../../../..
sinclude $(TOP)/config/Makefile

# Dependencies
particle_mod.o          : io_mod.o
star_mod.o              : io_mod.o particle_mod.o scaling_mod.o
task_mod.o		: io_mod.o trace_mod.o bits_mod.o
particles_mod.o         : io_mod.o particle_mod.o task_mod.o random_mod.o star_mod.o standard_disk_mod.o
standard_disk_mod.o     : io_mod.o trace_mod.o star_mod.o cgs_mod.o scaling_mod.o

# Executable name and dependency
OBJ = particles_mod.o star_mod.o
$(EXEC): $(OBJ)
