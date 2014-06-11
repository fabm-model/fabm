#
# Makefile to build Bolding & Burchard models
#

include ../../Rules.make

SUBDIRS += passive
SUBDIRS += filter_feeder

all: objs

objs: ${LIBFABM}(bb_model_library.o)

${LIBFABM}(bb_model_library.o): ${LIBFABM}(filter_feeder.o) ${LIBFABM}(passive.o)

${LIBFABM}(filter_feeder.o): filter_feeder
	$(MAKE) -C filter_feeder objs

${LIBFABM}(passive.o): passive
	$(MAKE) -C passive objs

doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: doc clean $(SUBDIRS)
