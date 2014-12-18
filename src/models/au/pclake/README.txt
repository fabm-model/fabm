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
This directory is similar to pclake directory which is the reference pclake model that kept
the original features of PCLake. The PCLake reference version will be the benchmark version
 and starting point of au-pclake.
This directory will be a dynamic one that will be under constant update with regular maintanence
and with further development. This directory currently includes:
1. au_pclake_model_library.f90, model libray for PCLake, make all the modules accessed by FABM
2. au_abiotic_water.f90, abiotic processes and state varialbles in the water column
3. au_abiotic_sediment.f90, abiotic processes and state varialbles in the sediment
4. au_phytoplankton_water.f90, phytoplankton groups(green algae,diatom and cyanobateria and 
   the related processes in water column
5. au_phytoplankton_sediment.f90, settled phytoplankton groups(green algae,diatom and
   cyanobateria) and the related processes in sediment
6. au_macrophytes.f90,submerged macrophytes and its related processes
7. au_foodweb_water.f90,zooplankton, juvenile white fish, adult white fish and piscivorous
   fish and the related processes in the water column
8. au_foodweb_sediment.f90,zoobenthos and the related processes in the sediment
9. au_auxilary.f90, module describes resuspension, sedimentation and burial processes. au-pcake
   updated the resuspension meothods which is related to real time shear stress
10.au_utility.f90, temperature functions used by different modules
11.au_dummy_state.f90, for the forcing state variable values, testing module for FABM functions.
!*-----------------------------------------------------------------------------*
!* If you have questions and suggestions regarding this model, please write to *
!*      Fenjuan Hu: fenjuan@bios.au.dk                                             *
!*      Dennis Trolle:dtr@bios.au.dk                                           *
!*      Karsten Bolding:bold@bios.au.dk                                        *
!*-----------------------------------------------------------------------------*
