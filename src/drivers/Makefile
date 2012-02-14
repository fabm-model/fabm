
include ../Rules.make

SUBDIRS += 0d getm gotm mom4p0d mom4p1

doc:
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

