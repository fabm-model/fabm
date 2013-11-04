#
# Makefile to build the PND model
#

include ../../../Rules.make

DOCSRC	= phy_feedback.F90

OBJS    = ${LIBFABM}(phy_feedback.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/klimacampus_phy_feedback.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Jorn Bruggeman (BB)                             !
#-----------------------------------------------------------------------
