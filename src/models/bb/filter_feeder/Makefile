#
# Makefile to build the filter_feeder model
#

include ../../../Rules.make

DOCSRC	= filter_feeder.F90

OBJS    = ${LIBFABM}(filter_feeder.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_filter_feeder.tex

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
