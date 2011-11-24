#
# Makefile to build the FABM-ERSEM interface
#

include ../../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	=  ersem.F90
ERSEMDIR = /common/work/ERSEM/FABM
ERSEMLIBDIR = $(LIBDIR)/
ERSEMLIB=$(ERSEMLIBDIR)libersem.a
ERSEMINCDIR = $(MODDIR)/
PPERSEM = -DMASSTRACER
OBJS   = \
${LIB}(ersem.o)

ifdef FABM_PMLERSEM
all:  $(ERSEMLIB) ${OBJS}
else
all: $(OBJS)
endif
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../../doc/fabm.tex
	touch doc

$(ERSEMLIB):
	$(MAKE) -C $(ERSEMDIR) incremental ERSEMLIB=$(ERSEMLIB) ERSEMLIBDIR=$(ERSEMLIBDIR) ERSEMINCDIR=$(ERSEMINCDIR) CPP=$(CPP) FC=$(FC) AC="$(AR) -r" ERSEMROOT=$(ERSEMDIR) PPDEFS=$(PPERSEM)

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
