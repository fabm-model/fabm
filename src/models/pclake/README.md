```sh
Developed by :
    Aquatic Ecology and Water Quality Management
    Wageningen University,Netherland
Implemented in FABM by:
     Bioscience-Silkeborg
     Aarhus University,Denmark
     Center of Lake Restoration(CLEAR)
     Southern Denmark University,Denmark
     Copyright by the PCLake_FABM group @AU-silkeborg
     under the GNU Public License - www.gnu.org
     Created December 2014
```
# Directory content:
1. pclake_model_library.f90, model libray for PCLake, make all the modules accessed
   by FABM
2. abiotic_water.f90, abiotic processes and state varialbles in the water column
3. abiotic_sediment.f90, abiotic processes and state varialbles in the sediment
4. phytoplankton_water.f90, phytoplankton groups(green algae,diatom and cyanobateria
   and the related processes in water column)
5. phytoplankton_sediment.f90, settled phytoplankton groups(green algae,diatom and 
   cyanobateria and the related processes in sediment
6. macrophytes.f90,submerged macrophytes and its related processes
7. foodweb_water.f90,zooplankton, juvenile white fish, adult white fish and piscivorous
   fish and the related processes in the water column
8. foodweb_sediment.f90,zoobenthos and the related processes in the water column
9. auxilary.f90, module describes resuspension, sedimentation and burial processes
10.utility.f90, temperature functions used by different modules
11.dummy_state.f90, for the forcing state variable values, testing module for FABM functions.
12. PCLake_original_VS_PCLake_fabm.png, this is the comparision of the fabm-pclake
   model output and pclake-original output(Benchmark test). The fabm-pclake with 0d
   set up reproduced the same results as original pclake model(0d based and only worked
   for 0d). This graph only showed phyto-dw(algae biomass), sdveg(macrophytes biomass),
   sdzoo(zooplankton biomass) and spo4w(phosphate in the water) as representative. But
   all the other model variables indicated the same model behavior for the two versions.
```sh
If you have questions and suggestions regarding this model, please write to 
   fenjuan Hu: fenjuan@bios.au.dk
   trolle@bios.au.dk
   bolding@bios.au.dk
```
# Important notes:
This version is the benchmark version that produce the same model results as the original excel verison.
# Few things has to be specified: 
- the Marsh module in the original PCLake model turned off
- resuspension function is not shear-stress related but kepth the original imperical function.
- water column burial correction is removed since this will be handled by hydrodynamic models
- dredging and fish harvesting is not included
- This version is static, served as reference implementation of PCLake,thus is not update after March, 2015
- If you want to refer and run latest FABM-PCLake model, please use the code under au/pclake
