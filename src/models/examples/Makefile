#$Id$
#
# Makefile to build the FABM example models
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

DOCSRC  += benthic_predator/benthic_predator.F90
SUBDIRS += benthic_predator
SUBDIRS += npzd
SUBDIRS += npzd_f2003

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/fabm.tex

clean:
