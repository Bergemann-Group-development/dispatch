# Append to VPATH
VPATH := $(VPATH):$(TOP)/io

# Dependencies
trace_mod.o:            io_mod.o omp_timer_mod.o omp_mod.o
io_mod.o:               omp_mod.o omp_timer_mod.o
output_mod.o:           io_mod.o
