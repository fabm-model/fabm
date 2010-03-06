include ../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	= \
rmbm.F90  \
rmbm_types.F90

MOD	= \
${LIB}(rmbm.o)

OBJ1   = \
drivers/gotm/rmbm_driver.F90

OBJ   = \
rmbm_types.o

SUBDIRS += models
SUBDIRS += drivers/gotm

OBJ: OBJ1

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
