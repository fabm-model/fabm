#
# Makefile to build the FABM models from IOW
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

DOCSRC  += ergom/ergom.F90
SUBDIRS += ergom

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/fabm.tex

clean:
	$(RM) ${LIB} doc

