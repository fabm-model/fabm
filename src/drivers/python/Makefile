include ../../Rules.make

export FABMHOST=python

ifeq ($(FORTRAN_COMPILER),IFORT)
SETUPARGS+=--fcompiler=intelem
endif

ifeq ($(FORTRAN_COMPILER),GFORTRAN)
SETUPARGS+=--fcompiler=gnu95
endif

OBJ = ${LIBFABM}(fabm_driver.o)
EXTRA_FFLAGS+=-fPIC

all: obj

obj: $(OBJ)

fabm:
	$(MAKE) EXTRA_FFLAGS=-fPIC -C $(FABMDIR)/src all

${LIBFABM}(helper.o): fabm

python: ${LIBFABM}(helper.o)
	./setup.py build_ext --inplace $(SETUPARGS)
clean:
	$(RM) *.o *.mod *~
	$(RM) -r build dist

realclean: clean
	$(RM) -r fabm_fortran

distclean: realclean
	$(MAKE) -C $(FABMDIR)/src $@

.PHONY: makedirs
