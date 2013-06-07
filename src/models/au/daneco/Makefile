#
# Makefile to build the au/daneco FABM compliant model
#

include ../../../Rules.make

DOCSRC	=  daneco.F90

OBJS    = ${LIBFABM}(daneco.o)

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
