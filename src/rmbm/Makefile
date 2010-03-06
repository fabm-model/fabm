include ../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	= \
rmbm.F90  \
rmbm_types.F90

MOD	= \
${LIB}(rmbm.o)

OBJ   = \
rmbm_types.o

SUBDIRS += models
SUBDIRS += drivers/gotm

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

all: $(MOD) $(OBJ) OBJS
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/rmbm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/gotm.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
