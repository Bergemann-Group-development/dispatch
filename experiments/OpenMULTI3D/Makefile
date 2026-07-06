# Option groups
#OPTS = full_debug
#OPTS = vscode
#OPTS = debug
#OPTS =
#OPTS = hoppe
OPTS = optimized

# Default target and paths
EXEC            = dispatch.x
TOP             = ../..

# Default methods
ifndef mpiifort
  MPI=
endif

OS := $(shell uname -s)

ifeq ($(OS),Darwin)
    COMPILER = gfort10
else
    COMPILER = ifort
endif

#COMPILER        = gfort10
#COMPILER        = ifort
SOLVER          = None
EOS             = microphysics/None
STAGGER         = 

sinclude $(TOP)/config/Makefile
VPATH := ./SRC:./HOP:$(VPATH) # Search SRC and HOP first
sinclude ./SRC/Makefile
sinclude ./HOP/Makefile

ifeq ($(COMPILER), ifort)
  OPT += -diag-disable=10448
endif


CLEANDIRS = {./SRC,./HOP,$(BUILD)}

# Symbolication for sampling profilers (sample/dtrace/Instruments). -g adds DWARF
# only and does NOT change -Ofast codegen; frame pointers help stack unwinding.
# OPT += -g -fno-omit-frame-pointer

# List the .o files the main program depends on after the colon:
$(EXEC): read_external_mod.o mpi_debug_mod.o setup_mod.o