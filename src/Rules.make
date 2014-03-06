#$Id$

SHELL   = /bin/sh

# The compilation mode is obtained from $COMPILATION_MODE - 
# default production - else debug or profiling
ifndef COMPILATION_MODE
compilation=production
else
compilation=$(COMPILATION_MODE)
endif

DEFINES=-D$(FORTRAN_COMPILER) -D_FABM_F2003_

FEATURES	=
FEATURE_LIBS	=
EXTRA_LIBS	=
INCDIRS		=
LDFLAGS		=

FABMBASE = ${LIBFABM}(fabm_driver.o) ${LIBFABM}(fabm_standard_variables.o) ${LIBFABM}(fabm_properties.o) ${LIBFABM}(fabm_types.o) ${LIBFABM}(fabm_particle.o) ${LIBFABM}(fabm_expressions.o)

#
# phony targets
#
.PHONY: clean realclean distclean dummy

# Top of this version of FABM.
ifndef FABMDIR
FABMDIR := $(HOME)/FABM/fabm-git
endif
ifeq ($(wildcard $(FABMDIR)/src/fabm.F90),)
$(error the directory FABMDIR=$(FABMDIR) is not a valid FABM directory)
endif

ifndef FABMHOST
FABMHOST = gotm
endif

CPP	= /lib/cpp

FEATURES += fabm
FEATURE_LIBS += -lfabm$(buildtype)
INCDIRS	+= -I$(FABMDIR)/src/drivers/$(FABMHOST)

# Directory related settings.

ifndef BINDIR
BINDIR	= $(FABMDIR)/bin
endif

ifndef LIBDIR
LIBDIR	= $(FABMDIR)/lib/$(FABMHOST)/$(FORTRAN_COMPILER)
endif

LIBFABM = $(LIBDIR)/libfabm$(buildtype).a

ifndef MODDIR
MODDIR	= $(FABMDIR)/modules/$(FABMHOST)/$(FORTRAN_COMPILER)
endif
INCDIRS	+= -I/usr/local/include -I$(FABMDIR)/include -I$(MODDIR)

ifdef FABM_PMLERSEM
DEFINES += -DFABM_PMLERSEM
endif

ifdef FABM_F2003
$(warning Fortran 2003 support is now on by default; there is no need to set FABM_F2003)
endif

ifdef FABM_NO_F2003
$(error Fortran 2003 support is now required; use of FABM_NO_F2003 is no longer supported)
endif

# Normally this should not be changed - unless you want something very specific.

# The Fortran compiler is determined from the environment variable FORTRAN_COMPILER.
# The value of this variable must match the extension of one of the compilers/compiler.* files.
# Currently tested:
#   IFORT:     Intel Fortran compiler 12.1 and up
#   GFORTRAN:  GNU Fortran compiler 4.7 and up
#   FTN:       Cray Fortran compiler 8.2
#   PGFORTRAN: Portland Group (PGI) Fortran compiler 13.2

# Sets options for debug compilation
ifeq ($(compilation),debug)
buildtype = _debug
DEFINES += -DDEBUG $(STATIC)
FLAGS   = $(DEBUG_FLAGS)
endif

# Sets options for profiling compilation
ifeq ($(compilation),profiling)
buildtype = _prof
DEFINES += -DPROFILING $(STATIC)
FLAGS   = $(PROF_FLAGS)
endif

# Sets options for production compilation
ifeq ($(compilation),production)
buildtype = _prod
DEFINES += -DPRODUCTION $(STATIC)
FLAGS   = $(PROD_FLAGS)
endif

include $(FABMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

# For making the source code documentation.
PROTEX	= protex -b -n -s

.SUFFIXES:
.SUFFIXES: .F90

LINKDIRS	+= -L$(LIBDIR)

CPPFLAGS	= $(DEFINES) $(INCDIRS)
FFLAGS  	= $(DEFINES) $(FLAGS) $(MODULES) $(INCDIRS) $(EXTRAS)
F90FLAGS  	= $(FFLAGS)
LDFLAGS		+= $(FFLAGS) $(LINKDIRS)

#
# Common rules
#
%.o: %.F90
#	@ echo "Compiling $<"
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
