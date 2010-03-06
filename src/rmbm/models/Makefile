#$Id: Makefile,v 1.28 2009-11-20 08:16:56 kb Exp $
#
# Makefile to build the Bio Module - libbio.a
#

include ../../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC  += npzd/npzd.F90
SUBDIRS += npzd

DOCSRC  += jellyfish/jellyfish.F90
SUBDIRS += jellyfish

DOCSRC  += co2sys/co2sys.F90
SUBDIRS += co2sys

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/rmbm.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
