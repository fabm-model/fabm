#
# Makefile to build the config modules
#

include ../Rules.make

CONFIGOBJ = ${LIBFABM}(fabm_config.o)
CONFIGTYPESOBJ = ${LIBFABM}(fabm_config_types.o)
YAMLOBJ = ${LIBFABM}(fabm_yaml.o)

OBJS = $(CONFIGTYPESOBJ) $(YAMLOBJ) $(CONFIGOBJ)

all: objs

$(YAMLOBJ): $(CONFIGTYPESOBJ)

$(CONFIGOBJ): $(YAMLOBJ) $(CONFIGTYPESOBJ)

objs: ${OBJS}
	$(MOVE_MODULES_COMMAND)

doc:

clean:
	$(RM) *.o *~

.PHONY: all objs doc clean

#-----------------------------------------------------------------------
# Copyright (C) 2013 - Jorn Bruggeman                                  !
#-----------------------------------------------------------------------
