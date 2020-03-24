#OPTS = full_debug
OPTS = debug
#OPTS =
#OPTS = optimized

TOP            = ../../..
EXEC           = run
SOLVER         = ramses/hydro

sinclude $(TOP)/config/Makefile

run: mesa_task.x mesa_code.x
	./mesa_task.x

mesa_code.x: io_unit_mod.o
mesa_task.x: io_unit_mod.o
