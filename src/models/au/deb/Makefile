#
# Makefile to build the au/deb_mussel FABM compliant model
#

include ../../../Rules.make

DOCSRC = deb.F90

OBJS   = ${LIBFABM}(deb.o)

all: objs

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc: $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/fabm.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Karsten Bolding (BB), Marie Maar (AU)           !
#-----------------------------------------------------------------------
