#$Id$
#
# Makefile to build the FABM models that have been ported from GOTM
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

DOCSRC  += npzd/npzd.F90
SUBDIRS += npzd

DOCSRC  += fasham/fasham.F90
SUBDIRS += fasham

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/fabm.tex 

clean:
