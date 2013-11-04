#
# Makefile to build the prey_predator model
#

include ../../../Rules.make

DOCSRC	= prey_predator.F90

OBJS    = ${LIBFABM}(prey_predator.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/prey_predator.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2013 - Fenjuan Hu                                      !
#-----------------------------------------------------------------------
