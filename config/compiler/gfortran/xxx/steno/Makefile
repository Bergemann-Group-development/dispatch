# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
F90 = mpif90 -x f95-cpp-input -ffree-line-length-0 $(W)
FC  = mpif90 -x f95-cpp-input -ffree-line-length-0 $(W)
LD  = mpif90
else
F90 = gfortran -x f95-cpp-input -ffree-line-length-0 $(W)
FC  = gfortran -x f95-cpp-input -ffree-line-length-0 $(W)
LD  = gfortran
endif
