#
# Makefile to build KlimaCampus models
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

SUBDIRS += phy_feedback

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

clean:
