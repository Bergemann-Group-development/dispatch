ifeq ($(MESA_DIR),)
ifeq ($($MESA_DIR_INTENTIONALLY_EMPTY),)
  $(error MESA_DIR enviroment variable is not set)
endif
endif

include $(MESA_DIR)/star/work_standard_makefile
#include /lustre/astro/ltw906/MESA/mesa-r9793/star/work_standard_makefile

STAR = star


## ADDED SSJ:

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job
BIN_DIR = ../../../../../bin


#####
mesa_star_mod.o mesa_star_mod.mod: $(OBJS) $(WORK_DIR)/mesa_star_mod.f90
	$(LOADER) $(FCopenmp) -c $(WORK_DIR)/mesa_star_mod.f90 -I$(MESA_DIR)/include/ -I$(MESA_DIR)/star/make/ -J$(BIN_DIR) -o $(BIN_DIR)/mesa_star_mod.o

