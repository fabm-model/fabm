include ../../Rules.make

DOCSRC	= fabm_driver.F90

OBJ = ${LIBFABM}(fabm_driver.o)

all: makedirs obj

obj: $(OBJ)

makedirs:
	$(MAKE) -C $(FABMDIR)/src $@

doc:    $(DOCSRC)
#	$(PROTEX) $(DOCSRC) > ../../../doc/driver_code_mossco.tex

clean:
	$(RM) -f *.o *~

.PHONY: makedirs
