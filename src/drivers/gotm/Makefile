include ../../Rules.make

DOCSRC	= fabm_driver.F90

OBJ = ${LIBFABM}(fabm_driver.o)

all: obj

obj: $(OBJ) makedirs

makedirs:
	$(MAKE) -C $(FABMDIR)/src $@

all: obj

obj: $(OBJ)
        $(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/driver_code_gotm.tex

clean:
	$(RM) -f *.o *~

.PHONY: makedirs
