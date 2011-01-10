include ../../../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	= \
gotm_rmbm.F90

OBJ   = \
${LIB}(gotm_rmbm.o)

OBJ1 = \
${LIB}(rmbm_driver.o)

all: $(OBJ1) $(OBJ)
        $(MOVE_MODULES_COMMAND)

driver: $(OBJ1)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/rmbm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/gotm.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
