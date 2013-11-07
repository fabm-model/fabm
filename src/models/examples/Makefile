#
# Makefile to build the FABM example models
#

SUBDIRS += benthic_predator
SUBDIRS += npzd
SUBDIRS += npzd_f2003
SUBDIRS += duplicator
SUBDIRS += mean

#all: OBJS
#	$(MOVE_MODULES_COMMAND)

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)

