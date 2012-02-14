#
# Makefile to build Bolding & Burchard models
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

SUBDIRS += passive

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

clean:
