#
# Makefile to build the FABM model library - libfabm.a
#

SUBDIRS += bb/
SUBDIRS += examples/
SUBDIRS += gotm/
SUBDIRS += iow/
SUBDIRS += klimacampus/
SUBDIRS += metu/
SUBDIRS += pml/

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@; done

realclean: clean

distclean: realclean

.PHONY: objs doc clean $(SUBDIRS)

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BB)          !
#-----------------------------------------------------------------------
