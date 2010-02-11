#$Id: Makefile,v 1.4 2007-09-25 10:06:10 kbk Exp $
#
# Makefile to build the the airsea interaction library - libairsea.a
#

include ../Rules.make

LIB	= ../../libairsea$(buildtype).a
LIB	:= $(LIBDIR)/libairsea$(buildtype).a

DOCSRC	=\
airsea.F90			\
airsea_variables.F90		\
humidity.F90			\
back_radiation.F90		\
airsea_fluxes.F90		\
kondo.F90			\
fairall.F90			\
short_wave_radiation.F90	

MOD	= \
${LIB}(airsea_variables.o)	\
${LIB}(airsea.o)

OBJ	= \
${LIB}(short_wave_radiation.o)	\
${LIB}(humidity.o) 		\
${LIB}(back_radiation.o) 	\
${LIB}(airsea_fluxes.o) 	\
${LIB}(kondo.o) 		\
${LIB}(fairall.o)

all: $(MOD) $(OBJ)
	$(MOVE_MODULES_COMMAND)

$(OBJ): $(MOD)

tst_airsea: $(MOD)
	$(FC) -o tst_airsea tst_airsea.F90 $(MOD)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/airsea.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/airsea.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 1999 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
