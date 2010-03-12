#$Id$
#
# Makefile to build the Seagrass Module - libseagrass.a
#

include ../../Rules.make

LIB	= $(LIBDIR)/libseagrass$(buildtype).a

DOCSRC	= seagrass.F90

OBJ   = \
${LIB}(seagrass.o)

all: ${OBJ}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/seagrass.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 1999 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
