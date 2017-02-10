FC = fort

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  F90 = mpi$(FC) -cpp -assume byterecl -diag-disable remark
  LD  = mpi$(FC)
else
  F90 = $(FC) -cpp -assume byterecl -diag-disable remark
  LD  = $(FC)
endif
