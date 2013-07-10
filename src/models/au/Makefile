#
# Makefile to build Bolding & Burchard models
#

SUBDIRS += prey_predator

#all: objs

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)
