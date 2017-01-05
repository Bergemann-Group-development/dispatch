# If using mpi and openmp add link option for threadsafe mpi library
ifeq ($(strip $(DO_OPENMP)),yes)
ifeq ($(strip $(DO_MPI)),yes)
LDFLAGS += -mt_mpi 
endif
endif
