HW = -axAVX -xSSE4.1 -diag-disable remark

# Depending on MPI use mpif90 or not
ifeq ($(strip $(MPI)),mpi)
  FC  = mpif08
  F90 = $(FC) -cpp -assume byterecl -diag-disable remark -DPREC=4
  LD  = $(FC)
else
  FC  = ifort
  F90 = $(FC) -cpp -assume byterecl -diag-disable remark -DPREC=4
  LD  = $(FC)
endif
