include ../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	= \
gotm_rmbm.F90

MOD	= \
${LIB}(rmbm.o)

OBJ   = \
gotm_rmbm.o
rmbm_driver.o

all: $(MOD) $(OBJ)
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/rmbm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/gotm.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
