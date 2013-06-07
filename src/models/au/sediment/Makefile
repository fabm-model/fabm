#
# Makefile to build the au/daneco_sediment FABM compliant model
#

include ../../../Rules.make

DOCSRC	= sediment.F90

OBJS    = ${LIBFABM}(sediment.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/fabm.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Karsten Bolding (BB), Karen Timmerman (AU)      !
#-----------------------------------------------------------------------
