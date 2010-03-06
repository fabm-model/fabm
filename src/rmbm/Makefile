include ../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	= \
rmbm.F90  \
rmbm_types.F90

MOD	= \
${LIB}(rmbm.o)

OBJ   = \
rmbm_types.o \
drivers/gotm/rmbm_driver.F90

SUBDIRS += models
SUBDIRS += drivers/gotm

OBJS: $(OBJ) $(SUBDIRS) 
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

all: $(MOD) OBJS
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/rmbm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/gotm.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
