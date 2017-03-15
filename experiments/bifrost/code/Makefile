VPATH := $(VPATH):code:$(BIFROST):$(BIFROST)/UNITS:$(BIFROST)/RADIATION:$(BIFROST)/BOUNDARIES:$(BIFROST)/EXTRAS

sinclude $(BIFROST)/EXPERIMENTS/dispatch

# filling in blanks

MAIN		 ?=
PARAMS		 ?= params
PARAMS_IO	 ?= params_split
MP		 ?= mpi
UNITS		 ?= solar
TIMESTEP	 ?=
MHD		 ?= mhd_$(MP)
MHD_VAR		 ?= mhd_var_$(MP)
PDE		 ?= pde_$(MP)
STAGGER		 ?= stagger_mesh_$(MP)
QUENCH		 ?= quench3_$(MP)
BOUNDARY	 ?= $(EXPERIMENT)_boundaries
BOUNDARY_SUB 	 ?= standard_boundaries
IO		 ?= io
IOC		 ?=
CIO		 ?=
RADIATION	 ?= radiation_GS_$(MP)
EOS		 ?= square_gas_$(MP)
DRIVING		 ?= driving_none_$(MP)
TIMING		 ?= timing_$(MP)
MATH		 ?= math_$(MP)
EXTRAS		 ?= $(EXPERIMENT)_extras_$(MP)
EXTRAS_SUB	 ?=
ADAPTION	 ?= adaption_empty
WORK		 ?= work_$(MP)
FFT		 ?= fftpack 
INFNAN		 ?= infnan

BASE 		?= $(MAIN) $(PARAMS_IO) $(PARAMS) $(IO) $(MP) $(UNITS) $(STAGGER) $(TIMESTEP) $(PDE) $(MHD) $(MHD_VAR) $(WORK) $(EOS) $(QUENCH) $(MATH) $(IOC) $(CIO) $(BOUNDARY) $(BOUNDARY_SUB) $(RADIATION) $(EXTRAS) $(EXTRAS_SUB) $(TIMING) $(ADAPTION) $(COLMASS) $(FFT) $(INFNAN)

OBJ		 ?= $(addsuffix .o,$(strip $(BASE)))

dispatch_mod.o: params.o mpi.o $(UNITS).o stagger_mesh_mpi.o $(TIMESTEP).o

sinclude $(BIFROST)/Makefile.dep
