#
# Makefile to build the FABM-ERSEM interface
#

include ../../../Rules.make

DOCSRC	=  ersem.F90

ERSEMDIR = /common/work/ERSEM/FABM
ERSEMLIBDIR = $(LIBDIR)/
ERSEMLIB=$(ERSEMLIBDIR)libersem.a
ERSEMINCDIR = $(MODDIR)/
PPERSEM = -DMASSTRACER
OBJS    = ${LIBFABM}(ersem.o)

all: objs

${OBJS}: $(FABMBASE)

ifdef FABM_PMLERSEM
objs: $(ERSEMLIB) ${OBJS}
else
objs: $(OBJS)
endif
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/pml_ersem.tex

$(ERSEMLIB):
	$(MAKE) -C $(ERSEMDIR) incremental ERSEMLIB=$(ERSEMLIB) ERSEMLIBDIR=$(ERSEMLIBDIR) ERSEMINCDIR=$(ERSEMINCDIR) CPP=$(CPP) FC=$(FC) AC="$(AR) -r" ERSEMROOT=$(ERSEMDIR) PPDEFS=$(PPERSEM)

clean:
	$(RM) *.o *~

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
