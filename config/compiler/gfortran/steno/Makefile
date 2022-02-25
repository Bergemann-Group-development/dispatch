# Depending on MPI use mpif90 or not
HDF5=
ifeq ($(strip $(MPI)),mpi)
  FC  = mpif90 -I${BUILD} -J${BUILD}
endif
F90 = $(FC) -x f95-cpp-input -ffree-line-length-0 $(W)
LD  = $(FC)
