.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
objdir = .obj
srcdir = src
libdir = src/libs

# Command-line options at make call
debug ?= 0 

## GFORTRAN OPTIONS ##
FC  = gfortran
INC_NC  = -I/opt/local/include
LIB_NC  = -L/opt/local/lib -lnetcdff -lnetcdf
INC_COORD = -I/Users/robinson/models/EURICE/coord/.obj
LIB_COORD = /Users/robinson/models/EURICE/coord/libcoordinates.a

FLAGS  = -I$(objdir) -J$(objdir) $(INC_COORD) $(INC_NC) 
LFLAGS = $(LIB_COORD) $(LIB_NC)

DFLAGS = -O3
ifeq ($(debug), 1)   # Debugging options
    DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
endif


###############################################
##							
## Rules for individual libraries or modules
##
###############################################

$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/stress.o: $(srcdir)/stress.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/viscosity.o: $(srcdir)/viscosity.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<


$(objdir)/yelmo_topography.o: $(srcdir)/yelmo_topography.f90 $(objdir)/lithosphere.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/yelmo_dynamics.o: $(srcdir)/yelmo_dynamics.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/yelmo_energy.o: $(srcdir)/yelmo_energy.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/yelmo_tracers.o: $(srcdir)/yelmo_tracers.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/yelmo_exchange.o: $(srcdir)/yelmo_exchange.f90 $(objdir)/yelmo.o $(objdir)/lithosphere.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/yelmo_io.o: $(srcdir)/yelmo_io.f90 $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/control.o: control.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/yelmo.o: $(srcdir)/yelmo.f90 $(objdir)/yelmo_topography.o \
	               $(objdir)/yelmo_dynamics.o $(objdir)/yelmo_energy.o \
	               $(objdir)/yelmo_tracers.o \
	               $(objdir)/yelmo_io.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

###############################################
##							
## List of grisli related files
## And local yelmo-grisli rules 
##
###############################################

# # Define grisli source directory 
# grisli_srcdir = GRISLI/SOURCES

# # Include grisli compilation rules
# include $(grisli_srcdir)/Makefile_grisli.mk

# grisli_common = $(objdir)/runparam_mod.o $(objdir)/3D-physique-gen_mod.o

# grisli_domain = $(objdir)/paradim-hemin40_mod.o $(objdir)/geography-hemin40_mod.o

# $(objdir)/yelmo_grisli.o: $(srcdir)/yelmo_grisli.f90 $(objdir)/yelmo.o $(grisli_domain) \
# 						  $(grisli_common)
# 	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

###############################################
##							
## List of yelmo related files
##
###############################################

yelmo_libs = $(objdir)/nml.o $(objdir)/ncio.o $(objdir)/control.o

yelmo_physics = $(objdir)/lithosphere.o $(objdir)/viscosity.o \
				$(objdir)/energy.o 

yelmo_base = $(objdir)/yelmo_topography.o  \
	   $(objdir)/yelmo_dynamics.o \
	   $(objdir)/yelmo_energy.o $(objdir)/yelmo_tracers.o \
	   $(objdir)/yelmo_exchange.o $(objdir)/yelmo.o \
	   $(objdir)/yelmo_io.o 

yelmo_grisli = $(objdir)/yelmo_grisli.o $(grisli_common) $(grisli_domain)

###############################################
##							
## Compilation of complete programs
##
###############################################

yelmo: $(yelmo_libs) $(yelmo_physics) $(yelmo_base) 
	$(FC) $(DFLAGS) $(FLAGS) -o test_yelmo.x $^ test_yelmo.f90 $(LFLAGS)
	@echo " "
	@echo "    test_yelmo.x is ready."
	@echo " "

test_energy: $(objdir)/energy.o $(objdir)/viscosity.o
	$(FC) $(DFLAGS) $(FLAGS) -o tests/test_energy.x $^ tests/test_energy.f90 $(LFLAGS)
	@echo " "
	@echo "    tests/test_energy.x is ready."
	@echo " "

test_stress: $(objdir)/stress.o $(objdir)/viscosity.o
	$(FC) $(DFLAGS) $(FLAGS) -o tests/test_stress.x $^ tests/test_stress.f90 $(LFLAGS)
	@echo " "
	@echo "    tests/test_stress.x is ready."
	@echo " "

test_indices:
	$(FC) $(DFLAGS) $(FLAGS) -o tests/test_indices.x $^ tests/test_indices.f90 $(LFLAGS)
	@echo " "
	@echo "    tests/test_indices.x is ready."
	@echo " "


yelmo-grl: $(yelmo_libs) $(yelmo_base) 
	$(FC) $(DFLAGS) $(FLAGS) -o test_yelmo-grl.x $^ test_yelmo-grl.f90 $(LFLAGS)
	@echo " "
	@echo "    test_yelmo-grl.x is ready."
	@echo " "

ygrisli-ant: $(yelmo_libs) $(yelmo_base) $(yelmo_grisli)
	$(FC) $(DFLAGS) $(FLAGS) -o test_ygrisli-ant.x $^ test_yelmo-ant.f90 $(LFLAGS)
	@echo " "
	@echo "    test_ygrisli-ant.x is ready."
	@echo " "


.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make yelmo-grl    : compiles the program test_yelmo-grl.x (ice sheet only)"
	@echo " make ygrisli-ant  : compiles the program test_ygrisli-ant.x (ice sheet only)"
	@echo " make clean        : cleans object files"
	@echo ""

clean:
	rm -f  *.x $(objdir)/*.o $(objdir)/*.mod gmon.out
	rm -rf *.x.dSYM   

