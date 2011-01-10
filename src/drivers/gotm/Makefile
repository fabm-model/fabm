include ../../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	= 

OBJ   = \
${LIB}(rmbm_driver.o)

all: $(OBJ)
        $(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/rmbm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
