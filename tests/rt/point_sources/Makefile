
# Option groups
#OPTS = full_debug
#OPTS = debug
#OPTS =
OPTS = optimized

# Default target and methods	
SOLVER    = AN/hlld
RT_SOLVER = point_sources
EXEC      = test.x
TOP       = ../../..

sinclude $(TOP)/config/Makefile

# List the .o files the main program depends on after the colon:
$(EXEC) : patch_mod.o point_sources_mod.o
