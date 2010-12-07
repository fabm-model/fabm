#
# Makefile to build the 0D version of the jellyfish bio-module
#

include ../../../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	=  jellyfish.F90

OBJS   = \
${LIB}(jellyfish.o)

all: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/rmbm.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
