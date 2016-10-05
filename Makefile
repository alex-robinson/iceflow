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

$(objdir)/iceflow.o: $(srcdir)/iceflow.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/viscosity.o: $(srcdir)/viscosity.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

###############################################
##							
## List of files to make library
##
###############################################

iceflow_obj = $(objdir)/nml.o $(objdir)/ncio.o $(objdir)/iceflow.o

###############################################
##							
## Compilation of complete programs
##
###############################################

# coordinates static library - using subset2
libiceflow-static: $(iceflow_obj)
	ar rc libiceflow.a $^
	@echo " "
	@echo "    libiceflow.a is ready."
	@echo " "

# coordinates shared library - using subset2
libiceflow-shared: $(iceflow_obj)
	$(FC) $(DFLAGS) $(FLAGS) -shared -fPIC -o libiceflow.so $^ $(LFLAGS)
	@echo " "
	@echo "    libiceflow.so is ready."
	@echo " "

test_iceflow: $(iceflow_obj) 
	$(FC) $(DFLAGS) $(FLAGS) -o test_iceflow.x $^ $(testdir)/test_iceflow.f90 $(LFLAGS)
	@echo " "
	@echo "    test_iceflow.x is ready."
	@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make libiceflow-static : compiles the static library libiceflow.a"
	@echo " make libiceflow-shared : compiles the shared library libiceflow.so"
	@echo " make test_iceflow : compiles the program test_iceflow.x"
	@echo " make clean        : cleans object files"
	@echo ""

clean:
	rm -f  *.x $(objdir)/*.o $(objdir)/*.mod gmon.out
	rm -rf *.x.dSYM   

