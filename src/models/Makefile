#
# Makefile to build the FABM model library - libfabm.a
#

SUBDIRS += au/           # University of Aarhus
SUBDIRS += bb/           # Bolding & Burchard
SUBDIRS += examples/     # Test examples
SUBDIRS += gotm/         # Converted from old GOTM/BIO models
SUBDIRS += iow/          # Institute for Ostsee Forchung
SUBDIRS += klimacampus/  # Hamburg ??
SUBDIRS += metu/         # Middle East Technical University
SUBDIRS += pml/          # Plymouth Marine Laboratory
SUBDIRS += aed/          # AED models
SUBDIRS += hzg/          # Helmholtz-Zentrum Geesthacht
SUBDIRS += msi/          # Marine Systems Institute, Tallinn University of Technology 

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@; done

realclean: clean

distclean: realclean

.PHONY: objs doc clean $(SUBDIRS)

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BB)          !
#-----------------------------------------------------------------------
