include ../../Rules.make

export FABMHOST=python

ifeq ($(FORTRAN_COMPILER),IFORT)
SETUPARGS+=--fcompiler=intel
endif

ifeq ($(FORTRAN_COMPILER),GFORTRAN)
SETUPARGS+=--fcompiler=gnu95
endif

EXTRA_FFLAGS+=-fPIC

all: python

fabm:
	$(MAKE) EXTRA_FFLAGS=-fPIC -C $(FABMDIR)/src all

${LIBFABM}(helper.o): fabm

python: ${LIBFABM}(helper.o)
	./setup.py build_ext --inplace $(SETUPARGS)

clean:
	$(RM) *.o *.mod *~
	$(RM) -rf build dist

realclean: clean
	$(RM) -rf fabm_fortran

distclean: realclean
	$(MAKE) -C $(FABMDIR)/src $@

.PHONY: all fabm python clean realclean distclean
