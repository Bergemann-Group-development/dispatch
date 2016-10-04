# Option groups
#OPTS = debug
#OPTS = full_debug
OPTS = optimized

# Common configuration
EXEC = test.x
TOP = ../../..
include $(TOP)/config/Makefile

test: test.x
	mkdir -p performance
	( uname -a; grep "model name" /proc/cpuinfo | head -1; ./test.x ) | tee -a performance/$(HOST).log

test.x: io_unit_mod.o
