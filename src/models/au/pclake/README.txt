!*                                                                             *
!* Originally developed by :                                                   *
!*     Aquatic Ecology and Water Quality Management                            *
!*     Wageningen University,Netherland                                        *
!* Implemented to FABM and further deveoped by:                                *
!*     Bioscience-Silkeborg                                                    *
!*     Aarhus University,Denmark                                               *
!*     Center of Lake Restoration(CLEAR)                                       *
!*     Southern Denmark University,Denmark                                     *
!*     Copyright and maintaned by the PCLake_FABM group @AU-silkeborg          *
!*     under the GNU Public License - www.gnu.org                              *
!*                                                                             *
!*      Created December 2014                                                  *
!*-----------------------------------------------------------------------------*
This directory is based on $FABMDIR$/src/models/pclake(PCLake reference). The
PCLake reference version will be the benchmark version and starting point of
au-pclake$FABMDIR$/src/models/au/pclake/).
This directory will be a dynamic one that will be under constant update with re-
-gular maintenance and with further development. This directory includes:
1. pclake_model_library.F90, model library for PCLake, make all the modules acce-
  -ssed by FABM.
2. abiotic_water.F90, abiotic processes and state variables in the water column.
3. abiotic_sediment.F90, abiotic processes and state variables in the sediment.
4. phytoplankton_water.F90, phytoplankton (green algae,diatom and blue-greens
   and the related processes in water column.
5. phytoplankton_sediment.F90, settled phytoplankton(green algae,diatom and blue
   -greens)and the related processes in sediment.
6. macrophytes.F90,submerged macrophytes and its related processes.
7. zooplankton.F90,zooplankton in the water column and its related processes.
8. fish.F90, juvenile white fish, adult white fish and piscivorous fish and their
   related processes in the water column.
8. zoobenthos.F90,zoobenthos and the related processes in the sediment
9. auxilary.f90, module describes resuspension, sedimentation and burial proces-
   -ses. au-pcake updated the resuspension method which is related to real time
   shear stress.
10.utility.F90, temperature functions used by different modules.
!*-----------------------------------------------------------------------------*
!* If you have questions and suggestions regarding this model, please write to *
!*      Fenjuan Hu: fenjuan@bios.au.dk                                         *
!*      Dennis Trolle:dtr@bios.au.dk                                           *
!*      Karsten Bolding:bold@bios.au.dk                                        *
!*-----------------------------------------------------------------------------*
!*-----------------------------------------------------------------------------*
!* Development log                                                             *
!*      The history of code changes and its discription                        *
!*-----------------------------------------------------------------------------*
November 17th, 2015
 removed dilution process in auxilary.F90, updated according yaml file
November 18th, 2015
implemented atomospheric deposition process, updated according file
January 7th, 2016
separated fish and zooplankton in pclake_foodweb_water module:
* File foodweb_water.F90 is replaced by zooplankton.F90,fish.F90.
* All zooplankton state variables(sDZoo,sNZoo,sPZoo),external dependencies and rel-
  ated zooplankton processes are kept in zooplankton.F90.
* All fish state variables(sDFiAd,sNFiAd,sPFiAd,sDFiJv,sNFiJv,sPFiJv,sDPisc),exter-
  nal dependencies and related fish processes are kept in fish.F90.
* changed related variable/parameter conventions:
  ** replaced *web* by *zoo* or *fish*
  ** in zooplankton module, wDWebZoo, wNWebZoo and wPWebZoo are replaced by wDZoo
     wNZoo and wPZoo
  ** in fish moudle,
     wDWebFiJv,tDWebFiJv,wNWebFiJv,tNWebFiJv,wPWebFiJv,tPWebFiJv
     wDWebFiAd,tDWebFiAd,wNWebFiAd,tNWebFiAd,wPWebFiAd,tPWebFiAd,wDWebPisc,
     tDWebPisc,wNWebPisc,tNWebPisc,wPWebPisc,tPWebPisc
     are replaced by
     wDFiJv,tDFiJv,wNFiJv,tNFiJv,wPFiJv,tPFiJv
     wDFiAd,tDFiAd,wNFiAd,tNFiAd,wPFiAd,tPFiAd,wDPisc,
     tDPisc,wNPisc,tNPisc,wPPisc,tPPisc
March 2nd, 2016
changed foodweb_sediment to zoobenthos module:
** replaced *web* by *ben*
March 20th, 2016
Implemented fish manipulation process: register fish maninpulation rate as a 
state variable and then can be input as time series file. Fish manipulation can 
be turned on and off, the according yaml file is updated.
April 20th, 2016, implemented modular fluexes diagnostics
modular fluxes: the total change of a state variable from the whole module. 
                for example, sNH4W can be changed in module abiotic_water, 
                abiotic_sediment, phytoplankton_water, zooplankton, fish 
                and auxilary. Each modules total comtribution is called
                modular fluxes for sNH4W.
