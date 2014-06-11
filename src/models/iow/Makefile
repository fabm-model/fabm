#
# Makefile to build the FABM models from IOW
#

include ../../Rules.make

SUBDIRS += age
SUBDIRS += ergom
SUBDIRS += spm

all: objs

objs: ${LIBFABM}(iow_model_library.o)

${LIBFABM}(iow_model_library.o): ${LIBFABM}(iow_age.o) ${LIBFABM}(ergom_cgt.o) ${LIBFABM}(spm.o)
#${LIBFABM}(iow_model_library.o): ${LIBFABM}(spm.o)
#${LIBFABM}(iow_model_library.o): ${LIBFABM}(ergom_cgt.o)
#${LIBFABM}(iow_model_library.o): ${LIBFABM}(iow_age.o) 

${LIBFABM}(iow_age.o): age
	$(MAKE) -C age objs

${LIBFABM}(ergom_cgt.o): ergom
	$(MAKE) -C ergom objs

${LIBFABM}(spm.o): spm
	$(MAKE) -C spm objs

doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: doc clean $(SUBDIRS)
