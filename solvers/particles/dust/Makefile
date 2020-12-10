# Dependencies

# Search path
ifneq ($(DUST_EVOLVE),)
  VPATH := $(TOP)/solvers/particles/$(PARTICLE_SOLVER)/$(DUST_EVOLVE):$(VPATH)
  sinclude $(TOP)/solvers/particles/$(PARTICLE_SOLVER)/$(DUST_EVOLVE)/Makefile
endif

particles_mod.o         : io_mod.o trace_mod.o omp_mod.o mpi_mod.o \
                          mpi_buffer_mod.o omp_lock_mod.o
particle_solver_mod.o   : omp_timer_mod.o timer_mod.o task_mod.o link_mod.o \
                          gpatch_mod.o particles_mod.o dust_mod.o list_mod.o \
                          omp_lock_mod.o timer_mod.o dust_evolve_mod.o
dust_mod.o              : gpatch_mod.o particles_mod.o random_mod.o \
                          scaling_mod.o scratch_mod.o
dust_evolve_mod.o       : trace_mod.o io_mod.o math_mod.o scaling_mod.o \
                          units_mod.o particles_mod.o

