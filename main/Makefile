################################################################################
#
#   Makefile to compile and link Fortran programs.
#
#  "make" compiles and links the specified main programs and modules
#  using the specified libraries (if any), and produces the executables;
#  by default it deletes also object files after producing executables.
# 
#  "make clean" removes all files apart from source codes.
#
################################################################################

all: rmxeq cmpff mkfxeq rmpmk
.PHONY: all

SHELL=/bin/bash


#############################  Fortran programs  ###############################

MAIN = eikonal2d_mpi  eikonal2d_pe
#theta_diffs
#short_corridor

GLOBAL = kind constants

CLUSTER =
#jackknife

EXTRAS = string_extras
#extras

LEARNING =
#policy_gradient

RANDOM = randomf

MODULES = $(GLOBAL) $(CLUSTER) $(LEARNING) $(EXTRAS) $(RANDOM)


##############################  set search paths  ##############################

# search path for modules
MDIR = ../modules
VPATH = $(MDIR)/global:$(MDIR)/cluster:$(MDIR)/simulations:$(MDIR)/start:\
	$(MDIR)/random:$(MDIR)/extras:$(MDIR)/learning


#########################  flags for Fortran programs  #########################

INCPATH = /usr/include ../include

LIBPATH = /usr/local/lib


#FC=gfortran
#FFLAGS= -C -g3 -O3 -DSSE4 -ffpe-summary='none' -fbacktrace -ffpe-trap=zero,overflow,underflow -fbounds-check
# -Wall -Wextra

FC = mpifort
FFLAGS = -C -g3 -O3 -DSSE4 -ffpe-summary='none' -Wunused-variable


FLIBS =
#lapack  blas 

PROGRAMS = $(MODULES) $(MAIN)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(FLIBS))


# compile source programs
$(addsuffix .o,$(PROGRAMS)): %.o: %.f90 Makefile
	@ $(FC) $(LDFLAGS) $(FFLAGS) $(INCDIRS) -c $<

# link object files and produce executables
$(MAIN): %: %.o $(OBJECTS)
	@ $(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^


############################  List of instructions  ############################

# remove old executables and old error log file
rmxeq:
	@ -rm -f $(MAIN); \
#        echo "  delete old executables"


# compile source files and produce objects
cmpff: $(addsuffix .o,$(PROGRAMS))


# link objects and produce executables
mkfxeq: $(MAIN)
#	@ echo "  link object files and produce executables"


# remove files except for executables and .mod
rmpmk:
	@ -rm -rf *.d *.o *.mod *.mod0 .tmp; \
#		echo -e "  delete object and temporary files\n"


# clean the main directory
.PHONY: clean
clean: rmpmk
	@ rm -f *.mod *.MOD *~ $(PROGRAMS); \
#		echo "  main directory cleaned"


################################################################################
