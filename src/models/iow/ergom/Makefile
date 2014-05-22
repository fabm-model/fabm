#
# Makefile to build the FABM version of the GOTM implementation of 
# the ERGOM model
#

include ../../../Rules.make

DOCSRC	=  ergom_cgt.F90

OBJS    = ${LIBFABM}(ergom_cgt.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/gotm_ergom.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
