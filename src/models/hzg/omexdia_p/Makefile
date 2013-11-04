#
# Makefile to build the omexdia_p model
#

include ../../../Rules.make

DOCSRC	= omexdia_p.F90

OBJS    = ${LIBFABM}(omexdia_p.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/examples_omexdia_p.tex

clean:
	$(RM) *.o *~

