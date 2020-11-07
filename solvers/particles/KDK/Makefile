# dependencies:

patch_force_mod.o       : trace_mod.o dll_mod.o
patch_forces_mod.o      : patch_force_mod.o
particle_history_mod.f90: trace_mod.o io_mod.o
particle_patch_mod.o    : particle_history_mod.o omp_timer_mod.o
particle_mod.o          : patch_forces_mod.o io_unit_mod.o lagrange_mod.o
particles_mod.o         : particle_mod.o patch_mod.o hash_table_mod.o \
                          scaling_mod.o scratch_mod.o
particle_solver_mod.o   : particles_mod.o
