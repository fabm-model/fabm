#$Id: Makefile,v 1.28 2009-11-20 08:16:56 kb Exp $
#
# Makefile to build the Bio Module - libbio.a
#

include ../../Rules.make

LIB	= $(LIBDIR)/libbio$(buildtype).a

DOCSRC	=  bio.F90                   \
           bio_var.F90               \
           npzd/bio_npzd.F90         \
           fasham/bio_fasham.F90     \
           iow/bio_iow.F90           \
           photo/bio_photo.F90       \
           bio_fluxes.F90            \
           bio_alloc_info.F90        \
           bio_save.F90   

OBJ1   = \
${LIB}(bio_var.o)		\
${LIB}(bio_alloc_info.o)        \
${LIB}(bio_fluxes.o)

OBJ2   = \
${LIB}(bio_save.o)		\
${LIB}(bio.o)

#BIO_TEMPLATE=false
#BIO_NPZD=false
#BIO_CL=false
#BIO_IOW=false
#BIO_SED=false
#BIO_FASHAM=false
BIO_MAB=false
#BIO_ROLM=false
#BIO_NPZD_FE=false
#BIO_PHOTO=false

ifneq ($(BIO_TEMPLATE),false)
SUBDIRS += template
DEFINES += -DBIO_TEMPLATE
endif

ifneq ($(BIO_NPZD),false)
DOCSRC  += npzd/bio_npzd.F90
SUBDIRS += npzd
DEFINES += -DBIO_NPZD
endif

ifneq ($(BIO_CL),false)
DOCSRC  += cl/bio_cl.F90
SUBDIRS += cl
DEFINES += -DBIO_CL
endif

ifneq ($(BIO_IOW),false)
DOCSRC  += iow/bio_iow.F90  
SUBDIRS += iow
DEFINES += -DBIO_IOW
endif

ifneq ($(BIO_SED),false)
DOCSRC  += sed/bio_sed.F90 
SUBDIRS += sed
DEFINES += -DBIO_SED
endif

ifneq ($(BIO_FASHAM),false)
DOCSRC  += fasham/bio_fasham.F90 
SUBDIRS += fasham
DEFINES += -DBIO_FASHAM
endif

ifneq ($(BIO_MAB),false)
DOCSRC  += 
SUBDIRS += mab
DEFINES += -DBIO_MAB
endif

ifneq ($(BIO_ROLM),false)
DOCSRC  += 
SUBDIRS += rolm
DEFINES += -DBIO_ROLM
endif

ifneq ($(BIO_NPZD_FE),false)
DOCSRC  += 
SUBDIRS += npzd_fe
DEFINES += -DBIO_NPZD_FE
endif

ifneq ($(BIO_PHOTO),false)
DOCSRC  += photo/bio_photo.F90  
SUBDIRS += photo
DEFINES += -DBIO_PHOTO
endif



all: $(OBJ1) OBJS $(OBJ2)
	$(MOVE_MODULES_COMMAND)

OBJS: $(SUBDIRS)
	set -e; for i in $(SUBDIRS); do $(MAKE) -C $$i ; done

$(OBJ2): OBJS

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/bio.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
