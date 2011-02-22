#$Id: Makefile 103 2010-12-15 12:38:33Z jornbr $
#
# Makefile to build the Bio Module - libbio.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC  += gotm/npzd/npzd.F90
SUBDIRS += gotm/npzd

DOCSRC  += pml/ersem/pmlersem.F90
SUBDIRS += pml/ersem

DOCSRC  += metu/mnemiopsis/mnemiopsis.F90
SUBDIRS += metu/mnemiopsis

DOCSRC  += pml/co2sys/co2sys.F90
SUBDIRS += pml/co2sys

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/fabm.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
