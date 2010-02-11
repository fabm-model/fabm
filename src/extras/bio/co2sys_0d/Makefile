#
# Makefile to build the 0D version of the carbonate system module
#

include ../../../Rules.make

LIB	= $(LIBDIR)/libbio$(buildtype).a

DOCSRC	=  bio_co2_sys_0d.F90

OBJS   = \
${LIB}(co2_dyn.o)	\
${LIB}(bio_co2_sys_0d.o)

all: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/bio.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
