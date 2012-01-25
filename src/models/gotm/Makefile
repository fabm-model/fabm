#
# Makefile to build the FABM models that have been ported from GOTM
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

SUBDIRS += npzd

SUBDIRS += fasham

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done
#	$(PROTEX) $(DOCSRC) > ../../../doc/driver_gotm_code.tex

clean:
