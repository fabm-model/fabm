ifndef GOTMDIR
FABMDIR = $(HOME)/FABM/fabm-git
endif

all: exec

exec:
	( export FABMHOST=0d ; $(MAKE) -C $(FABMDIR)/src/drivers/0d $@ )
	mv $(FABMDIR)/src/drivers/0d/fabm0d .

