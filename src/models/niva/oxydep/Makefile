#
# Makefile to build oxydep
#

include ../../../Rules.make

DOCSRC	= oxydep.F90

OBJS    = ${LIBFABM}(oxydep.o)

all: objs

${OBJS}: $(FABMBASE) ${LIBFABM}(base.o)

${LIBFABM}(base.o): $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/niva.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Jorn Bruggeman (BB)                             !
#-----------------------------------------------------------------------
