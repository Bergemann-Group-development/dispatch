# Append to VPATH
VPATH := $(TOP)/solvers/rt/$(RT_SOLVER):$(VPATH)

# Dependencies
radiation_mod.o:   io_mod.o trace_mod.o mesh_mod.o link_mod.o \
                   cooling_Dalgarno.f90 \
                   cooling_Rosen_etal.f90 \
                   cooling_Casiana.f90 \
                   cooling_Gnedin.f90
