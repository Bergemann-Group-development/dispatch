HW = -axAVX -xSSE4.1 -diag-disable remark

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  FC  = mpif08
else
  FC  = ifort
endif
