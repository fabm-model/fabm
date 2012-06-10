#
# Makefile to build the FABM version of the iow module
#

include ../../../Rules.make

DOCSRC	=  ergom.F90

OBJS    = ${LIBFABM}(ergom.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/iow_ergom.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
