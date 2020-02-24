TOP = ../../..

sinclude $(TOP)/config/Makefile

default:
	$(FC) mesa_task.f90 -o mesa_task.x
	$(FC) mesa_code.f90 -o mesa_code.x
	./mesa_task.x
