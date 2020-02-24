include $(MESA)/work_standard_makefile

OBJS = $(MESA)/standard_work_dir/make/run_star_extras.o $(MESA)/standard_work_dir/make/run_star.o mesa_star_mod.o

# Objects that have to linked in, and are compiled within the code
HYDROOBJ += stellar_evolution.o

# External MESA dependencies
LIBS += $(OBJS) $(LOAD_MESA)
