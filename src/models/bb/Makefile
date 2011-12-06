#
# Makefile to build Bolding & Burchard models
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

DOCSRC  += passive/passive.F90
SUBDIRS += passive

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/fabm.tex

clean:
