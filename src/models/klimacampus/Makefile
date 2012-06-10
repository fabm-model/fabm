#
# Makefile to build KlimaCampus models
#

SUBDIRS += phy_feedback

#all: objs

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)
