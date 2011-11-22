include Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	= \
fabm.F90  \
fabm_types.F90

FABMMOD	= ${LIB}(fabm.o)
FABMBASE = ${LIB}(fabm_types.o)
FABMMODELLIB = ${LIB}(fabm_library.o)

all: makedirs $(FABMMOD)
        $(MOVE_MODULES_COMMAND)

0d getm gotm:
	( export FABMHOST=$@ ; $(MAKE) )

makedirs:
	mkdir -p $(MODDIR) $(LIBDIR)

driver:
	$(MAKE) -C drivers/$(FABMHOST)

$(FABMBASE): driver

models: driver $(FABMBASE)
	$(MAKE) -C models

$(FABMMOD): $(FABMBASE) models $(FABMMODELLIB)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../doc/fabm.tex 
	touch doc

clean:
	$(RM) -f *.o *~

realclean: clean
	$(RM) -r $(FABMDIR)/modules/$(FABMHOST)/$(FORTRAN_COMPILER)
	$(RM) -r $(FABMDIR)/lib/$(FABMHOST)/$(FORTRAN_COMPILER)

distclean: realclean
	$(RM) -r doc
	$(RM) -r $(FABMDIR)/modules/ $(FABMDIR)/lib/
