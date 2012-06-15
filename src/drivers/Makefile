
include ../Rules.make

#SUBDIRS += 0d getm gotm mom4p0d mom4p1
SUBDIRS += 0d gotm

doc clean:
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

