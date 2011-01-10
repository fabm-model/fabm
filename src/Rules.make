#$Id$

SHELL   = /bin/sh

# The compilation mode is obtained from $COMPILATION_MODE - 
# default production - else debug or profiling
ifndef COMPILATION_MODE
compilation=production
else
compilation=$(COMPILATION_MODE)
endif

DEFINES=-D$(FORTRAN_COMPILER)

FEATURES	=
FEATURE_LIBS	=
EXTRA_LIBS	=
INCDIRS		=
LDFLAGS		=

#
# phony targets
#
.PHONY: clean realclean distclean dummy

# Top of this version of RMBM.
ifndef RMBMDIR
RMBMDIR  := $(HOME)/rmbm
endif

CPP	= /lib/cpp

DEFINES += -DRMBM
FEATURES += rmbm
FEATURE_LIBS += -lrmbm$(buildtype)
INCDIRS	+= -I$(RMBMDIR)/src/drivers/gotm

# Directory related settings.

ifndef BINDIR
BINDIR	= $(RMBMDIR)/bin
endif

ifndef LIBDIR
LIBDIR	= $(RMBMDIR)/lib/$(FORTRAN_COMPILER)
endif

ifndef MODDIR
MODDIR	= $(RMBMDIR)/modules/$(FORTRAN_COMPILER)
endif
INCDIRS	+= -I/usr/local/include -I$(RMBMDIR)/include -I$(MODDIR)

# Normaly this should not be changed - unless you want something very specific.

# The Fortran compiler is determined from the EV FORTRAN_COMPILER - options 
# sofar NAG(linux), FUJITSU(Linux), DECF90 (OSF1 and likely Linux on alpha),
# SunOS, PGF90 - Portland Group Fortran Compiler (on Intel Linux).

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

include $(RMBMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

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
ifeq  ($(can_do_F90),true)
%.o: %.F90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
else
%.f90: %.F90
#	$(CPP) $(CPPFLAGS) $< -o $@
	$(F90_to_f90)
%.o: %.f90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
endif
