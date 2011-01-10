include Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	= \
rmbm.F90  \
rmbm_types.F90

RMBMMOD	= ${LIB}(rmbm.o)
RMBMBASE = ${LIB}(rmbm_types.o)

all: gotmrmbm
        $(MOVE_MODULES_COMMAND)

driver:
	$(MAKE) -C drivers/gotm driver

$(RMBMBASE): driver

models: driver $(RMBMBASE)
	$(MAKE) -C models

$(RMBMMOD): $(RMBMBASE) models

gotmrmbm: $(RMBMMOD)
	 $(MAKE) -C drivers/gotm

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../doc/rmbm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
