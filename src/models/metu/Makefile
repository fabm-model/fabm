#
# Makefile to build the FABM models from the Middle East Technical University (METU)

SUBDIRS += mnemiopsis

#all: OBJS
#       $(MOVE_MODULES_COMMAND)

objs doc clean:
	set -e; for dir in $(SUBDIRS); do $(MAKE) -C $$dir $@ ; done

.PHONY: objs doc clean $(SUBDIRS)
