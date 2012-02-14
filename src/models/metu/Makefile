#
# Makefile to build the FABM models from the Middle East Technical University (METU)

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

SUBDIRS += mnemiopsis

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

clean:
	$(RM) ${LIB} doc

