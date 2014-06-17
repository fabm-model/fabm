include ../../Rules.make

export FABMHOST=python

EXTRA_FFLAGS+=-fpic

all: python

fabm:
	$(MAKE) EXTRA_FFLAGS=-fPIC -C $(FABMDIR)/src all

${LIBFABM}(helper.o): fabm

python: ${LIBFABM}(helper.o)
	$(FC) -shared -fpic $(LDFLAGS) $(LIBS) fabm_python.F90 $(LIBFABM) -o fabm-python.so

clean:
	$(RM) *.o *.mod *~

realclean: clean
	$(RM) *.so

distclean: realclean
	$(MAKE) -C $(FABMDIR)/src $@

.PHONY: all fabm python clean realclean distclean
