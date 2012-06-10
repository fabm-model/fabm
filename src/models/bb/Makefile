#
# Makefile to build Bolding & Burchard models
#

SUBDIRS += passive

#all: objs

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)
