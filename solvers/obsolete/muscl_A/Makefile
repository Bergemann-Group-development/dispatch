
mhd_mod.o       : io_mod.o trace_mod.o riemann_mod.o omp_lock_mod.o extras_mod.o
riemann_mod.o   : io_mod.o timer_mod.o
initial_mod.o   : io_unit_mod.o trace_mod.o random_mod.o io_mod.o mesh_mod.o \
                  index_mod.o units_mod.o
solver_mod.o    : io_mod.o trace_mod.o kinds_mod.o patch_mod.o mhd_mod.o stagger_2nd.o
force_mod.o     : io_mod.o trace_mod.o mesh_mod.o random_mod.o scaling_mod.o
