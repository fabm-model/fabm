#$Id: Rules.make,v 1.22 2009-11-20 08:16:56 kb Exp $

SHELL   = /bin/sh

# The compilation mode is obtained from $COMPILATION_MODE - 
# default production - else debug or profiling
ifndef COMPILATION_MODE
compilation=production
else
compilation=$(COMPILATION_MODE)
endif

DEFINES=-DNUDGE_VEL
DEFINES=-D$(FORTRAN_COMPILER)

DEFINES += -DLOCATIONTYPE=integer -DLOCATIONVARIABLE=ilev

# What do we include in this compilation
NetCDF=false
NetCDF=true
SEDIMENT=false
#SEDIMENT=true
SEAGRASS=false
SEAGRASS=true
BIO=false
BIO=true
RMBM=false
RMBM=true

FEATURES	=
FEATURE_LIBS	=
EXTRA_LIBS	=
INCDIRS		=
LDFLAGS		=

# If we want NetCDF - where are the include files and the library

ifeq ($(NetCDF),true)

DEFINES += -DNETCDF_FMT

ifdef NETCDFINC
INCDIRS         += -I$(NETCDFINC)
endif

ifdef NETCDFLIBDIR
LINKDIRS        += -L$(NETCDFLIBDIR)
endif

ifdef NETCDFLIBNAME
NETCDFLIB       = $(NETCDFLIBNAME)
else
NETCDFLIB       = -lnetcdf -lnetcdff
endif

ifeq ($(NETCDF_VERSION),NETCDF4)

DEFINES         += -DNETCDF4
ifdef HDF5_DIR
INCDIRS         += -I$(HDF5_DIR)/include
LINKDIRS        += -L$(HDF5_DIR)/lib
endif
HDF5LIB         = -lhdf5_hl -lhdf5 -lz

else  # NetCDF3 is default

DEFINES         += -DNETCDF3
HDF5LIB         =

endif

EXTRA_LIBS      += $(NETCDFLIB) $(HDF5LIB)

endif
# NetCDF/HDF configuration done

#
# phony targets
#
.PHONY: clean realclean distclean dummy

# Top of this version of GOTM.
ifndef GOTMDIR
GOTMDIR  := $(HOME)/gotm
endif

CPP	= /lib/cpp

# Here you can put defines for the [c|f]pp - some will also be set depending
# on compilation mode.
ifeq ($(SEDIMENT),true)
DEFINES += -DSEDIMENT
FEATURES += extras/sediment
FEATURE_LIBS += -lsediment$(buildtype)
endif
ifeq ($(SEAGRASS),true)
DEFINES += -DSEAGRASS
FEATURES += extras/seagrass
FEATURE_LIBS += -lseagrass$(buildtype)
endif
ifeq ($(BIO),true)
DEFINES += -DBIO
FEATURES += extras/bio
FEATURE_LIBS += -lbio$(buildtype)
endif
ifeq ($(NO_0D_BIO),true)
DEFINES         += -DNO_0D_BIO
endif

# Directory related settings.

ifndef BINDIR
BINDIR	= $(GOTMDIR)/bin
endif

ifndef LIBDIR
LIBDIR	= $(GOTMDIR)/lib/$(FORTRAN_COMPILER)
endif

ifndef MODDIR
MODDIR	= $(GOTMDIR)/modules/$(FORTRAN_COMPILER)
endif
INCDIRS	+= -I/usr/local/include -I$(GOTMDIR)/include -I$(MODDIR)

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

include $(GOTMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

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
