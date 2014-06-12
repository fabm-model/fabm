#
# Makefile to build the FABM example models
#

include ../../Rules.make

LIBRARY  = ${LIBFABM}(examples_model_library.o)

SUBDIRS += benthic_predator
BENTHIC_PREDATOR = ${LIBFABM}(benthic_predator.o)

SUBDIRS += duplicator
DUPLICATOR = ${LIBFABM}(duplicator.o)

SUBDIRS += mean
MEAN = ${LIBFABM}(mean.o)

SUBDIRS += npzd
NPZD = ${LIBFABM}(examples_npzd.o)

SUBDIRS += npzd_f2003
NPZD_F2003 = ${LIBFABM}(npzd_f2003.o)


all: objs

objs: $(LIBRARY)

$(LIBRARY): $(BENTHIC_PREDATOR) $(DUPLICATOR) $(MEAN) $(NPZD) $(NPZD_F2003)

$(BENTHIC_PREDATOR): benthic_predator
	$(MAKE) -C benthic_predator objs

$(DUPLICATOR): duplicator
	$(MAKE) -C duplicator objs

$(MEAN): mean
	$(MAKE) -C mean objs

$(NPZD): npzd
	$(MAKE) -C npzd objs

$(NPZD_F2003): npzd_f2003
	$(MAKE) -C npzd_f2003 objs

doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: doc clean $(SUBDIRS)

