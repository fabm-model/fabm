#
# Makefile to build the 0D version of the carbonate system module
#

include ../../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	=  carbonate.F90

OBJS   = \
${LIB}(co2_dyn.o)	\
${LIB}(carbonate.o)

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
