#$Id$
#
# Makefile to build the FABM models from the Middle East Technical University (METU)

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

DOCSRC  += mnemiopsis/mnemiopsis.F90
SUBDIRS += mnemiopsis

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/fabm.tex

clean:
	$(RM) ${LIB} doc

