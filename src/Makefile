include Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	= \
fabm.F90  \
fabm_types.F90

FABMMOD	= ${LIB}(fabm.o)
FABMBASE = ${LIB}(fabm_types.o)

all: makedirs $(FABMMOD)
        $(MOVE_MODULES_COMMAND)

makedirs:
	mkdir -p $(MODDIR) $(LIBDIR)

driver:
	$(MAKE) -C drivers/gotm

$(FABMBASE): driver

models: driver $(FABMBASE)
	$(MAKE) -C models

$(FABMMOD): $(FABMBASE) models

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../doc/fabm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
