#
# Makefile to build the FABM example models
#

include ../../Rules.make

SUBDIRS += benthic_predator
SUBDIRS += duplicator
SUBDIRS += mean
SUBDIRS += npzd
SUBDIRS += npzd_f2003

all: objs

objs: ${LIBFABM}(examples_model_library.o)

${LIBFABM}(examples_model_library.o): ${LIBFABM}(benthic_predator.o) ${LIBFABM}(duplicator.o) ${LIBFABM}(mean.o) ${LIBFABM}(det.o) ${LIBFABM}(nut.o) ${LIBFABM}(phy.o) ${LIBFABM}(zoo.o) ${LIBFABM}(npzd_f2003.o)

${LIBFABM}(benthic_predator.o): benthic_predator
	$(MAKE) -C benthic_predator objs

${LIBFABM}(duplicator.o): duplicator
	$(MAKE) -C duplicator objs

${LIBFABM}(mean.o): mean
	$(MAKE) -C mean objs

${LIBFABM}(det.o) ${LIBFABM}(nut.o) ${LIBFABM}(phy.o) ${LIBFABM}(zoo.o): npzd
	$(MAKE) -C npzd objs

${LIBFABM}(npzd_f2003.o): npzd_f2003
	$(MAKE) -C npzd_f2003 objs

doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: doc clean $(SUBDIRS)

