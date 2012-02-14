#
# Makefile to build the 0D version of the npzd bio-module
#

include ../../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	= \
det.F90 \
nut.F90 \
phy.F90 \
zoo.F90

OBJS   = \
${LIB}(det.o) \
${LIB}(nut.o) \
${LIB}(phy.o) \
${LIB}(zoo.o)

all: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_npzd.tex

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Casper Boon and Karsten Bolding (BB)            !
#-----------------------------------------------------------------------
