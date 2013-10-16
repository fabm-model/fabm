#
# Makefile to build the FABM models that have been ported from GOTM
#

SUBDIRS  = ergom
SUBDIRS += npzd
SUBDIRS += fasham

#all: objs

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)
