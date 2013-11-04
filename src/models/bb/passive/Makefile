#
# Makefile to build the passive tracer model
#

include ../../../Rules.make

DOCSRC	= passive.F90

OBJS    = ${LIBFABM}(passive.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/bb_passive.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Jorn Bruggeman (BB)                             !
#-----------------------------------------------------------------------
