# Append to VPATH
VPATH := $(VPATH):$(TOP)/solvers/particles:$(TOP)/solvers/particles/$(PARTICLE_SOLVER)

# Dependencies
particle_history_mod.o  : io_mod.o trace_mod.o
particle_patch_mod.o    : io_mod.o trace_mod.o particle_solver_mod.o dll_mod.o particle_history_mod.o
#sinkparticle_mod.o      : io_mod.o trace_mod.o particle_solver_mod.o patch_mod.o download_mod.o data_hub_mod.o
# Particle solver components
particle_mod.o          : io_unit_mod.o
particle_list_mod.o     : trace_mod.o particle_mod.o omp_timer_mod.o dll_mod.o patch_mod.o
particle_force_mod.o    : particle_list_mod.o
particle_solver_mod.o   : particle_list_mod.o
