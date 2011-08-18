#$Id$
#
# Makefile to build the FABM models from the Plymouth Marine Laboratory
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

DOCSRC  += carbonate/carbonate.F90 ersem/pml_ersem.F90
SUBDIRS += carbonate ersem

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/fabm.tex 

clean:
	$(RM) ${LIB} doc

