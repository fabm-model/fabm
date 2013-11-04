#
# Makefile to build the 0D version of the npzd bio-module
#

include ../../../Rules.make

DOCSRC	=  npzd.F90

OBJS    = ${LIBFABM}(npzd.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/gotm_npzd.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
