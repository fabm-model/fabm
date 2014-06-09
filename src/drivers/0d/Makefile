include ../../Rules.make

export FABMHOST=0d

all: fabm0d

fabm:
	$(MAKE) -C $(FABMDIR)/src

ifndef GOTMDIR
GOTMDIR = $(HOME)/GOTM/gotm-git
endif


FFLAGS	+= -I$(GOTMDIR)/include -I$(GOTMDIR)/src/util -I$(GOTMDIR)/modules/$(FORTRAN_COMPILER)
LDFLAGS += -L$(FABMDIR)/lib/$(FORTRAN_COMPILER) -L$(GOTMDIR)/lib/$(FORTRAN_COMPILER)
LIBS	+= -lairsea$(buildtype) -linput$(buildtype) -lutil$(buildtype) -lfabm$(buildtype)

ifeq ($(NETCDF_VERSION),NETCDF4)
INCDIRS += -DNETCDF4 -I$(shell nf-config --includedir)
LDFLAGS += $(shell nf-config --flibs)
endif

fabm0d.o: fabm shared.o output.o

output.o: shared.o

main.o: fabm0d.o

fabm0d: main.o
	$(FC) -o $@ main.o fabm0d.o shared.o output.o $(LDFLAGS) $(LIBS)
	$(RM) fabm0d.o main.o
	mkdir -p $(HOME)/bin
	cp -p $@ $(HOME)/bin

clean:
	$(RM) *.o *.mod *~

realclean: clean
	$(RM) fabm0d

distclean: realclean
	$(MAKE) -C $(FABMDIR)/src $@

.PHONY: all fabm fabm0d clean realclean distclean
