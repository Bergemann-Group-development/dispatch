HW      = -align array64byte -qopt-assume-safe-padding -fp-model source
PAR     = -qopenmp
LDFLAGS = 

ifeq ($(strip $(COMPILER)),intel)
	# Use intel MPI :: mpiifort
	ifeq ($(strip $(MPI)),mpi)
  	FC  = mpiifort 
	else
  	FC  = ifort
	endif
else
	# Use mpif90 if MPI is set 
	ifeq ($(strip $(MPI)),mpi)
  	FC  = mpif90
	else
  	FC  = ifort
	endif
endif

F90 = $(FC) -cpp -module $(BUILD)
LD = $(FC) -module $(BUILD) $(PAR)