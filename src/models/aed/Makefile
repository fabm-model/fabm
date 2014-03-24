#
# Makefile to build the aed water quality module
#

include ../../Rules.make

LIB	= $(LIBDIR)/libfabm$(buildtype).a

DOCSRC	= \
aed_util.F90 \
aed_sedflux.F90 \
aed_chlorophylla.F90 \
aed_oxygen.F90 \
aed_silica.F90 \
aed_carbon.F90 \
aed_nitrogen.F90 \
aed_phosphorus.F90 \
aed_organic_matter.F90 \
aed_phyto_utils.F90 \
aed_phytoplankton.F90 \
aed_zoop_utils.F90 \
aed_zooplankton.F90 \
aed_pathogens.F90 \
aed_iron.F90 \
aed_sulfur.F90 \
aed_tracer.F90 \
aed_totals.F90 \
aed_models.F90

OBJS   = \
${LIB}(aed_core.o) \
${LIB}(aed_util.o) \
${LIB}(aed_sedflux.o) \
${LIB}(aed_chlorophylla.o) \
${LIB}(aed_oxygen.o) \
${LIB}(aed_silica.o) \
${LIB}(aed_carbon.o) \
${LIB}(aed_nitrogen.o) \
${LIB}(aed_phosphorus.o) \
${LIB}(aed_organic_matter.o) \
${LIB}(aed_phyto_utils.o) \
${LIB}(aed_phytoplankton.o) \
${LIB}(aed_zoop_utils.o) \
${LIB}(aed_zooplankton.o) \
${LIB}(aed_pathogens.o) \
${LIB}(aed_iron.o) \
${LIB}(aed_sulfur.o) \
${LIB}(aed_tracer.o) \
${LIB}(aed_totals.o) \
${LIB}(aed_models.o)


all: objs
	$(MOVE_MODULES_COMMAND)

${OBJS}: $(FABMBASE)

objs: ${OBJS}

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/fabm.tex
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2008 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
