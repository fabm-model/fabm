#
# Makefile to build the 0D version of the npzd bio-module
#

include ../../../Rules.make

DOCSRC	= det.F90 nut.F90 phy.F90 zoo.F90

OBJS    = ${LIBFABM}(det.o) ${LIBFABM}(nut.o) ${LIBFABM}(phy.o) ${LIBFABM}(zoo.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_npzd.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Casper Boon and Karsten Bolding (BB)            !
#-----------------------------------------------------------------------
