#
# Makefile to build the FABM version of the fasham module
#

include ../../../Rules.make

DOCSRC	=  fasham.F90

OBJS    = ${LIBFABM}(fasham.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/gotm_fasham.tex 

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
