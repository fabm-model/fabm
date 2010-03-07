include ../Rules.make

LIB	= $(LIBDIR)/librmbm$(buildtype).a

DOCSRC	= \
rmbm.F90  \
rmbm_types.F90

MOD	= \
${LIB}(rmbm.o)

OBJ   = \
${LIB}(rmbm_types.o)

SUBDIRS += models
SUBDIRS += drivers/gotm

all: OBJS $(MOD)
        $(MOVE_MODULES_COMMAND)

driver:
	$(MAKE) -C drivers/gotm driver

$(OBJ): driver

$(SUBDIRS): $(OBJ)

OBJS: $(SUBDIRS) 
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/rmbm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/gotm.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean
