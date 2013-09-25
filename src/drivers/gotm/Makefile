include ../../Rules.make

ifndef GOTMDIR
GOTMDIR = $(HOME)/GOTM/gotm-git
endif
FFLAGS += -I$(GOTMDIR)/include

DOCSRC	= fabm_driver.F90

OBJ = ${LIBFABM}(fabm_driver.o)

all: obj

obj: $(OBJ) makedirs

$(OBJ): gotmversion

gotmversion:
	$(MAKE) -C $(GOTMDIR) VERSION

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
