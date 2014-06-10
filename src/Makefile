include Rules.make

SUBDIRS = models config

DOCSRC	= fabm.F90 fabm_types.F90

FABMMOD	 = ${LIBFABM}(fabm.o)
FABMMODELLIB = ${LIBFABM}(fabm_library.o)

all: makedirs models config
        $(MOVE_MODULES_COMMAND)

makedirs:
	mkdir -p $(MODDIR) $(LIBDIR)

models: $(FABMBASE)
	$(MAKE) -C models

config: $(FABMMOD)
	$(MAKE) -C config

$(FABMMOD): $(FABMMODELLIB)

0d getm gotm mossco python:
	( export FABMHOST=$@ ; $(MAKE) )

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../doc/fabm_interface_code.tex
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

clean:
	$(RM) -f *.o *~
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

realclean: clean
	$(RM) -r $(FABMDIR)/modules/$(FABMHOST)/$(FORTRAN_COMPILER)
	$(RM) -r $(FABMDIR)/lib/$(FABMHOST)/$(FORTRAN_COMPILER)

distclean: realclean
	$(RM) -r doc
	$(RM) -r $(FABMDIR)/modules/ $(FABMDIR)/lib/

.PHONY: all makedirs models config doc clean realclean distclean
