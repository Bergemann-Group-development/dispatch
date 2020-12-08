
# Option groups
OPTS = debug
#OPTS = full_debug
#OPTS = optimized

default: nml2dat.x

TOP    = ../..
EXEC   = 
sinclude $(TOP)/config/Makefile
VPATH := .:$(VPATH)
CP     = mv
COPY   = $(TOP)/bin

nml2dat.x            : nml2dat.f90
