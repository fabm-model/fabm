#
# Makefile to build the 0D version of the carbonate system module
#

include ../../../Rules.make

DOCSRC	= co2_dyn.F90 carbonate.F90

OBJS    = ${LIBFABM}(co2_dyn.o) ${LIBFABM}(carbonate.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/pml_carbonate.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
