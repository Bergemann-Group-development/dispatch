# If using mpi and openmp add link option for threadsafe mpi library
ifeq ($(strip $(DO_OPENMP)),yes)
ifeq ($(strip $(DO_MPI)),yes)
PAR += -mt_mpi 
endif
endif
