#
# Makefile to build the passive tracer model
#

include ../../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	= \
passive.F90

OBJS   = \
${LIB}(passive.o)

all: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/bb_passive.tex

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Jorn Bruggeman (BB)                             !
#-----------------------------------------------------------------------
