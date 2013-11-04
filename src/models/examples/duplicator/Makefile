#
# Makefile to build the passive tracer model
#

include ../../../Rules.make

DOCSRC	= duplicator.F90

OBJS    = ${LIBFABM}(duplicator.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_duplicator.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Jorn Bruggeman (BB)                             !
#-----------------------------------------------------------------------
