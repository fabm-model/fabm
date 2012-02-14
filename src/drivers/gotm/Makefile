include ../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	= fabm_driver.F90

OBJ   = \
${LIB}(fabm_driver.o)

all: $(OBJ)
        $(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/driver_code_gotm.tex

clean:
	$(RM) -f ${LIB} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
