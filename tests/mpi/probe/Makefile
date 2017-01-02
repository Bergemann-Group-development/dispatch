# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

# Common configuration
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

ifneq ($(MPI),mpi)
  $(error "This test *requires* MPI!")
endif

test.x:
