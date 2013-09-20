include ../../Rules.make

export FABMHOST=0d

OBJ = ${LIBFABM}(fabm_driver.o)

all: obj fabm0d

obj: makedirs $(OBJ)

0d makedirs:
	$(MAKE) -C $(FABMDIR)/src $@

ifndef GOTMDIR
GOTMDIR = $(HOME)/GOTM/gotm-git
endif
FFLAGS	+= -I$(GOTMDIR)/modules/$(FORTRAN_COMPILER)
LDFLAGS += -L$(FABMDIR)/lib/$(FORTRAN_COMPILER) -L$(GOTMDIR)/lib/$(FORTRAN_COMPILER)
LIBS	+= -lairsea$(buildtype) -linput$(buildtype) -lutil$(buildtype) -lfabm$(buildtype)

fabm0d: 0d fabm0d.o main.o
	$(FC) -o $@ main.o fabm0d.o $(LDFLAGS) $(LIBS)
	$(RM) fabm0d.o main.o
	cp -p $@ $(HOME)/bin

clean:
	$(RM) *.o *.mod *~

realclean: clean
	$(RM) fabm0d

distclean: realclean
	$(MAKE) -C $(FABMDIR)/src $@

.PHONY: makedirs
