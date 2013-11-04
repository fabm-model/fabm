#
# Makefile to build the 0D version of the GETM spm module
#

include ../../../Rules.make

DOCSRC	=  spm.F90

OBJS   = \
${LIBFABM}(spm.o)

all: objs

${OBJS}: $(FABMBASE)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/fabm.tex

realclean:
	$(RM) *.o *~

