#
# Makefile to build the FABM models from the Plymouth Marine Laboratory
#
include ../../Rules.make

SUBDIRS += oxydep

all: objs

objs: ${LIBFABM}(niva_model_library.o)

doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

${LIBFABM}(niva_model_library.o): subdirobjs

subdirobjs:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir objs ; done

.PHONY: objs doc clean subdirobjs (SUBDIRS)
