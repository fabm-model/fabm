#
# Makefile to build the FABM models from IOW
#

SUBDIRS  = spm
SUBDIRS += age
SUBDIRS += ergom

#all: objs
#	$(MOVE_MODULES_COMMAND)

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)
