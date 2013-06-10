#
# Makefile to build the au FABM models
#

SUBDIRS += daneco
SUBDIRS += deb
SUBDIRS += sediment

all: objs

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)

