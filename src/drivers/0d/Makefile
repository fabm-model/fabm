include ../../Rules.make

LIB = $(LIBDIR)/libfabm$(buildtype).a

OBJ = ${LIB}(fabm_driver.o)

all: $(OBJ)

ifndef GOTMDIR
GOTMDIR = $(HOME)/GOTM/gotm-git
endif
FFLAGS	+= -I$(GOTMDIR)/modules/$(FORTRAN_COMPILER)
LDFLAGS += -L$(FABMDIR)/lib/$(FORTRAN_COMPILER) -L$(GOTMDIR)/lib/$(FORTRAN_COMPILER)
LIBS	+= -lairsea$(buildtype) -lutil$(buildtype) -lfabm$(buildtype)

exec: fabm0d

fabm0d: $(OBJ) fabm0d.o main.o
	$(FC) -o $@ main.o fabm0d.o $(LDFLAGS) $(LIBS)
	$(RM) fabm0d.o main.o
	cp -p $@ $(HOME)/bin

clean:
	$(RM) fabm0d.o main.o fabm0d

realclean: clean
	$(RM) *.o *.mod *~

distclean: realclean
