#$Id: Makefile 103 2010-12-15 12:38:33Z jornbr $
#
# Makefile to build the Bio Module - libbio.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

SUBDIRS += gotm/

SUBDIRS += pml/

SUBDIRS += metu/

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
