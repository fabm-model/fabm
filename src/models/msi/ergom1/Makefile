#
# Makefile to build the ergom1 model
#

include ../../../Rules.make

DOCSRC	= ergom1.F90

OBJS    = ${LIBFABM}(ergom1.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/msi_ergom1.tex

clean:
	$(RM) *.o *~

