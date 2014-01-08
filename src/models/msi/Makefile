#
# Makefile to build the FABM models from MSI
#

SUBDIRS += ergom1

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean
