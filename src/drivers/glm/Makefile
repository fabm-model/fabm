include ../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	=

OBJ   = \
${LIB}(fabm_driver.o)

all: $(OBJ)
        $(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/fabm.tex
	touch doc

clean:
	$(RM) -f ${LIB} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
