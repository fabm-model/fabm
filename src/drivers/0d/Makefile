include ../../Rules.make

export FABMHOST=0d

all: fabm0d

ifndef GOTMDIR
GOTMDIR = $(HOME)/GOTM/gotm-git
endif

LIB0D	= $(FABMDIR)/lib/0d/$(FORTRAN_COMPILER)/lib0D.a

FFLAGS	+= -I$(GOTMDIR)/include -I$(GOTMDIR)/src/util -I$(GOTMDIR)/modules/$(FORTRAN_COMPILER)
LDFLAGS += -L$(FABMDIR)/lib/$(FORTRAN_COMPILER) -L$(GOTMDIR)/lib/$(FORTRAN_COMPILER)
LIBS	+= -l0D -lairsea$(buildtype) -linput$(buildtype) -lutil$(buildtype) -lfabm$(buildtype)

ifeq ($(NETCDF_VERSION),NETCDF4)
INCDIRS += -DNETCDF4 -I$(shell nf-config --includedir)
LIBS += $(shell nf-config --flibs)
endif

fabm0d: main.o
	$(FC) -o $@ main.o $(LDFLAGS) $(LIBS)
	$(RM) main.o
	mkdir -p $(HOME)/bin
	cp -p $@ $(HOME)/bin

main.o: ${LIB0D}(fabm0d.o)

${LIB0D}(fabm0d.o): fabm ${LIB0D}(output.o)

fabm:
	$(MAKE) -C $(FABMDIR)/src

${LIB0D}(output.o): ${LIB0D}(shared.o)

$(LIB0D): fabm

clean:
	$(RM) *.o *.mod *~ lib0D.a

realclean: clean
	$(RM) fabm0d

distclean: realclean
	$(MAKE) -C $(FABMDIR)/src $@

.PHONY: all clean realclean distclean
