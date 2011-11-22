#
# Makefile to build the FABM version of the iow module
#

include ../../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	=  ergom.F90

OBJS   = \
${LIB}(ergom.o)

all: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/fabm.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
