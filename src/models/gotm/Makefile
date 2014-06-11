#
# Makefile to build the FABM models that have been ported from GOTM
#

include ../../Rules.make

SUBDIRS  = ergom
SUBDIRS += fasham
SUBDIRS += npzd

all: objs

objs: ${LIBFABM}(gotm_model_library.o)

${LIBFABM}(gotm_model_library.o): ${LIBFABM}(ergom.o) ${LIBFABM}(fasham.o) ${LIBFABM}(npzd.o)

${LIBFABM}(ergom.o): ergom
	$(MAKE) -C ergom objs

${LIBFABM}(fasham.o): fasham
	$(MAKE) -C fasham objs

${LIBFABM}(npzd.o): npzd
	$(MAKE) -C npzd objs

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: doc clean $(SUBDIRS)
