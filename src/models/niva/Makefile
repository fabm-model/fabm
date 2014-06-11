#
# Makefile to build the FABM models from the Norsk Institutt for Vannforskning
#

include ../../Rules.make

SUBDIRS += oxydep

all: objs

objs: ${LIBFABM}(niva_model_library.o)

${LIBFABM}(niva_model_library.o): ${LIBFABM}(oxydep.o)

${LIBFABM}(oxydep.o): oxydep
	$(MAKE) -C oxydep objs

doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: doc clean $(SUBDIRS)
