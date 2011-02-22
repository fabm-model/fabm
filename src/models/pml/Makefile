#$Id$
#
# Makefile to build the pml FABM models
#

include ../../Rules.make

LIB     = $(LIBDIR)/libfabm$(buildtype).a

DOCSRC  += co2sys/co2sys.F90 ersem/pmlersem.F90
SUBDIRS += co2sys ersem

all: OBJS
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/fabm.tex 

clean:
	$(RM) ${LIB} doc

