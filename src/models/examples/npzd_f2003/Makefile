#
# Makefile to build the Fortran 2003 version of the npzd module
#

include ../../../Rules.make

DOCSRC	=  npzd_f2003.F90

OBJS    = ${LIBFABM}(npzd_f2003.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_npzd_f2003.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2011 - Jorn Bruggeman                                  !
#-----------------------------------------------------------------------
