#
# Makefile to build the FABM models from the Plymouth Marine Laboratory
#

SUBDIRS += carbonate

#all: objs

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)
