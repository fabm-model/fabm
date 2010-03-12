#$Id$

include ../../../Rules.make

CORE_LIBS	=	\
		-lairsea$(buildtype)		\
		-lmeanflow$(buildtype) 		\
		-lturbulence$(buildtype) 	\
		-lobservations$(buildtype)	\
		-loutput$(buildtype)		\
		-lutil$(buildtype)

ALL_LIBS	= $(FEATURE_LIBS) $(CORE_LIBS) $(EXTRA_LIBS)

EXEC = driver_0d$(buildtype)_$(FORTRAN_COMPILER)

SUBDIRS  = ../../..

all: $(EXEC)

$(EXEC): $(SUBDIRS)
	$(FC) -o $@ ./driver_0d.F90 ./main.F90 $(LDFLAGS) $(ALL_LIBS)

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: $(SUBDIRS)

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
