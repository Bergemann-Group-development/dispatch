FC = gfortran
LD = gfortran
FLAGS = -std=f2008ts
FLAGS += -g -ffpe-trap=zero,invalid,overflow,underflow -fbacktrace -fbounds-check

BUILD = build
FLAGS +=  -I$(BUILD) -J$(BUILD) 

FCFLAGS = $(FLAGS) -cpp
LFLAGS  = $(FLAGS)


# Object Files for build
SRCS = const_mod prec_mod op_math_mod handyfuncs_mod level_mod ion_mod molecule_mod element_mod transition_mod osmet_mod optab_mod composition_mod main
OBJS = $(SRCS:=.o)

MAIN = main

$(MAIN) : $(BUILD) $(OBJS)
	$(FC) $(FCFLAGS) $(sort $(patsubst %,$(BUILD)/%,$(OBJS))) -I$(BUILD) -o $(MAIN)

$(BUILD):
	mkdir -p $(BUILD)

.SUFFIXES : .o .f90	

%.o: %.f90
	$(LD) $(LFLAGS) -c $<  -o $(BUILD)/$@

clean :
	$(RM) *.o *.mod $(MAIN) $(BUILD)/*.o $(BUILD)/*.mod