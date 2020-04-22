!-----------------------------------------------------------------------
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Elizaveta Protsenko,
!                     Jorn Bruggeman
!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:
!
! !INTERFACE:

   module fabm_niva_brom_redox
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_redox
!     Variable identifiers
    type (type_state_variable_id)        :: id_Mn2,id_Mn3,id_Mn4,id_MnCO3
    type (type_state_variable_id)        :: id_H2S
    type (type_state_variable_id)        :: id_Fe2,id_Fe3,id_FeS,id_MnS,id_NO2
    type (type_state_variable_id)        :: id_S0,id_S2O3,id_SO4
    type (type_state_variable_id)        :: id_Baae,id_Bhae,id_Baan,id_Bhan
    type (type_state_variable_id)        :: id_O2,id_NO3,id_NH4
    type (type_state_variable_id)        :: id_PON,id_DON,id_PO4,id_Si,id_Sipart
    type (type_state_variable_id)        :: id_DIC,id_Alk,id_CaCO3,id_FeS2,id_FeCO3,id_CH4

    type (type_dependency_id)            :: id_temp,id_salt
    type (type_dependency_id)            :: id_Kp1,id_Kp2,id_Kp3,id_Knh4,id_Kh2s,id_Hplus,id_KSi,id_Kc0
    type (type_dependency_id)            :: id_Ca,id_CO3,id_Om_Ca,id_Om_Ar,id_pco2

    type (type_diagnostic_variable_id)  :: id_DcPM_O2,id_DcDM_O2,id_DcPM_NOX, id_DcDM_NOX,id_DcPM_SO4, id_DcDM_SO4
    type (type_diagnostic_variable_id)  :: id_DcDM_Fe,id_DcPM_Fe,id_DcDM_Mn,id_DcPM_Mn
    type (type_diagnostic_variable_id)  :: id_ChemBaae, id_ChemBaan,  id_HetBhan,id_HetBhae,id_MortBaan,id_MortBhan,id_MortBaae,id_MortBhae
    type (type_diagnostic_variable_id)  :: id_fe_ox1,id_autolysis,id_anammox,id_fe_ox2,id_feco3_diss,id_feco3_form !id_mns_form
    type (type_diagnostic_variable_id)  :: id_s2o3_rd_DM,id_s2o3_rd_PM,id_s2o3_no3,id_s0_no3,id_so4_rd_PM,id_so4_rd_DM,id_s0_ox,id_s2o3_ox
    type (type_diagnostic_variable_id)  :: id_Nitrif1, id_Nitrif2, id_fe_rd,id_Denitr1_PM,id_s0_disp,id_hs_ox,id_hs_no3,id_so4_rd,id_s2o3_rd
    type (type_diagnostic_variable_id)  :: id_fe_p_compl,id_mn_p_compl,id_fe_si_compl
    type (type_diagnostic_variable_id)  :: id_Denitr1_DM,id_Denitr2_PM,id_Denitr2_DM,id_Denitr1, id_Denitr2
    type (type_diagnostic_variable_id)  :: id_mn_ox2,id_mn_ox1,id_mn_rd1,id_mn_rd2,id_mns_diss,id_mns_form,id_mnco3_diss,id_mnco3_form
    type (type_diagnostic_variable_id)  :: id_DcPM_CH4,id_DcDM_CH4,id_ch4_o2,id_ch4_so4

!     Model parameters
       real(rk) :: Wsed, Wbact, Wm
      !===========================================================================!
      !-------------------------------------------------------------------------
      ! specific rates of biogeochemical processes
       !---- Mn---------!
      real(rk) :: K_mn_ox1, K_mn_ox2, K_mn_rd1, K_mn_rd2, K_mns, K_mns_diss, K_mns_form
      real(rk) :: K_mnco3, K_mnco3_diss, K_mnco3_form, K_mnco3_ox, K_DON_mn, K_PON_mn
      real(rk) :: s_mnox_mn2, s_mnox_mn3, s_mnrd_mn4, s_mnrd_mn3
      !---- Fe---------!
      real(rk) :: K_fe_ox1, K_fe_ox2, K_fe_rd, K_fes, K_fes_form, K_fes_diss, K_fes_ox, K_DON_fe,K_PON_fe
      real(rk) :: K_fes2_form, K_fes2_ox, s_feox_fe2, s_ferd_fe3, K_feco3, K_feco3_diss, K_feco3_form, K_feco3_ox
       !---- S ---------!
      real(rk) :: K_hs_ox, K_s0_ox, K_s2o3_ox, K_so4_rd, K_s2o3_rd, K_s0_disp, K_s0_no3, K_s2o3_no3, K_mnrd_hs, K_ferd_hs
       !---- N---------!
      real(rk) :: K_DON_ox, K_PON_ox, Tda, beta_da, K_omox_o2, K_PON_DON, K_nitrif1, K_nitrif2, K_denitr1, K_denitr2
      real(rk) :: K_omno_no3, K_omno_no2, K_hs_no3, K_annamox
       !---- O2--------!
      real(rk) ::  O2s_nf, O2s_dn, s_omox_o2, s_omno_o2, s_omso_o2, s_omso_no3, K_mnox_o2
       !---- C--------!
      real(rk) ::  K_caco3_diss, K_caco3_form, K_DON_ch4, K_PON_ch4, K_ch4_o2, K_ch4_so4, s_omch_so4
       !---- Si-------!
      real(rk) ::  K_sipart_diss
       !---- Bacteria-!
      real(rk) ::  K_Baae_gro, K_Baae_mrt, K_Baae_mrt_h2s, limBaae
      real(rk) ::  K_Bhae_gro, K_Bhae_mrt, K_Bhae_mrt_h2s, limBhae
      real(rk) ::  K_Baan_gro, K_Baan_mrt, limBaan
      real(rk) ::  K_Bhan_gro, K_Bhan_mrt, K_Bhan_mrt_o2, limBhan
      !===========================================================================!
      !---- Stoichiometric coefficients ----!
      real(rk) ::  r_n_p,r_o_n,r_c_n,r_fe_n,r_mn_n,f
      !---- Partitioning coefficients ----!
      real(rk) ::  r_fe3_p, r_mn3_p, r_fe3_si



   contains
      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_redox), intent(inout), target :: self
   integer,                      intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s):
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   !call self%get_parameter(self%w_bio,'w_bio','m/s','vertical velocity of phytoplankton (<0 for sinking)',default=-1.5_rk)
   !-----Model parameters------
   call self%get_parameter(self%Wsed,   'Wsed', '[1/day]',  'Rate of sinking of detritus (POP, PON)',           default=5.00_rk)
   call self%get_parameter(self%Wbact,  'Wbact','[1/day]',  'Rate of sinking of bacteria (Bhae,Baae,Bhan,Baan)',default=0.4_rk)
   call self%get_parameter(self%Wm,     'Wm',   '[1/day]',  'Rate of accelerated sinking of particles with settled Mn hydroxides', default=7.0_rk)
   ! specific rates of biogeochemical processes
   !---- Mn---------!
   call self%get_parameter(self%K_mn_ox1, 'K_mn_ox1', '[1/day]',  'Specific rate of oxidation of Mn2 to Mn3  with O2',    default=0.01_rk)
   call self%get_parameter(self%K_mn_ox2,'K_mn_ox2','[1/day]',  'Specific rate of oxidation of Mn3 to Mn4  with O2',    default=0.1_rk)
   call self%get_parameter(self%K_mn_rd1, 'K_mn_rd1', '[1/day]',  'Specific rate of reduction of Mn4 to Mn3  with H2S',   default=0.5_rk)
   call self%get_parameter(self%K_mn_rd2,'K_mn_rd2','[1/day]',  'Specific rate of reduction of Mn3 to Mn2  with H2S',   default=0.5_rk)
   call self%get_parameter(self%K_mns,   'K_mns',   '[M]',      'Conditional equilibrium constant for MnS from Mn2 with H2S',   default=0.02_rk)
   call self%get_parameter(self%K_mns_diss, 'K_mns_diss',   '[1/day]',  'Specific rate of dissolution of MnS to Mn2 and H2S',   default=0.0001_rk)
   call self%get_parameter(self%K_mns_form, 'K_mns_form',   '[1/day]',  'Specific rate of formation of MnS from Mn2 with H2S',   default=0.00001_rk)
   call self%get_parameter(self%K_mnco3,    'K_mnco3',      '[M]',      'Conditional equilibrium constant %  1.8e-11 ',     default=1000.0_rk)
   call self%get_parameter(self%K_mnco3_diss, 'K_mnco3_diss', '[1/day]', 'Specific rate of dissolution of MnCO3',   default=2.7e-7_rk)
   call self%get_parameter(self%K_mnco3_form, 'K_mnco3_form', '[1/day]', 'Specific rate of formation of MnCO3',    default=2.7e-7_rk)
   call self%get_parameter(self%K_mnco3_ox, 'K_mnco3_ox',   '[1/day]',  'Specific rate of oxidation of MnCO3 with O2',      default=0.0027_rk)
   call self%get_parameter(self%K_DON_mn,   'K_DON_mn',     '[1/day]',  'Specific rate of oxidation of DON with Mn4',       default=0.001_rk)
   call self%get_parameter(self%K_PON_mn,   'K_PON_mn',     '[1/day]',  'Specific rate of oxidation of PON with Mn4',       default=0.001_rk)
   call self%get_parameter(self%s_mnox_mn2, 's_mnox_mn2',   '[uM Mn]',  'threshold of Mn2 oxidation',       default=0.01_rk)
   call self%get_parameter(self%s_mnox_mn3, 's_mnox_mn3',   '[uM Mn]',  'threshold of Mn3 oxidation',       default=0.01_rk)
   call self%get_parameter(self%s_mnrd_mn4, 's_mnrd_mn4',   '[uM Mn]',  'threshold of Mn4 reduciton',       default=0.01_rk)
   call self%get_parameter(self%s_mnrd_mn3, 's_mnrd_mn3',   '[uM Mn]',  'threshold of Mn3 reduciton',       default=0.01_rk)
   !---- Fe---------!
   call self%get_parameter(self%K_fe_ox1,   'K_fe_ox1',  '[1/day]',  'Specific rate of oxidation of Fe2 to Fe3  with O2',    default=0.5_rk)
   call self%get_parameter(self%K_fe_ox2,   'K_fe_ox2', '[1/day]',  'Specific rate of oxidation of Fe2 to Fe3  with MnO2',  default=0.001_rk)
   call self%get_parameter(self%K_fe_rd,    'K_fe_rd',  '[1/day]',  'Specific rate of reduction of Fe3 to Fe2  with H2S',   default=0.5_rk)
   call self%get_parameter(self%K_fes,      'K_fes',    '[uM]',     'FeS equilibrium constant (Solubility Product Constant)',   default=2510.0_rk)
   call self%get_parameter(self%K_fes_form, 'K_fes_form', '[1/day]','Specific rate of precipitation of FeS from Fe2 with H2S',   default=5.e-5_rk)
   call self%get_parameter(self%K_fes_diss, 'K_fes_diss', '[1/day]','Specific rate of dissollution of FeS to Fe2 and H2S',   default=1.e-6_rk)
   call self%get_parameter(self%K_fes_ox,   'K_fes_ox', '[1/day]',  'Specific rate of oxidation of FeS with O2',   default=0.001_rk)
   call self%get_parameter(self%K_DON_fe,   'K_DON_fe', '[1/day]',  'Specific rate of oxidation of DON with Fe3',   default=0.00005_rk)
   call self%get_parameter(self%K_PON_fe,   'K_PON_fe', '[1/day]',  'Specific rate of oxidation of PON with Fe3',   default=0.00001_rk)
   call self%get_parameter(self%K_fes2_form,'K_fes2_form','[1/day]','Specific rate of FeS2 formation by FeS oxidation by H2S',   default=0.000001_rk)
   call self%get_parameter(self%K_fes2_ox,  'K_fes2_ox','[1/uM/d]', 'Specific rate of pyrite oxidation by O2',   default=0.00044_rk)
   call self%get_parameter(self%s_feox_fe2, 's_feox_fe2', '[uM Fe]','threshold of Fe2 reduciton',   default=0.001_rk)
   call self%get_parameter(self%s_ferd_fe3, 's_ferd_fe3', '[uM Fe]','threshold of Fe3 reduciton',   default=0.01_rk)
   call self%get_parameter(self%K_feco3,    'K_feco3',      '[M]',      'Conditional equilibrium constant %  1.8e-11 ',     default=1000.0_rk)
   call self%get_parameter(self%K_feco3_diss, 'K_feco3_diss', '[1/day]', 'Specific rate of dissolution of FeCO3',   default=2.7e-7_rk)
   call self%get_parameter(self%K_feco3_form, 'K_feco3_form', '[1/day]', 'Specific rate of formation of FeCO3',    default=2.7e-7_rk)
   call self%get_parameter(self%K_feco3_ox, 'K_feco3_ox',   '[1/day]',  'Specific rate of oxidation of FeCO3 with O2',      default=0.0027_rk)
   !---- S---------!
   call self%get_parameter(self%K_hs_ox, 'K_hs_ox', '[1/day]', 'Specific rate of oxidation of H2S to S0  with O2',  default=0.5_rk)
   call self%get_parameter(self%K_s0_ox, 'K_s0_ox', '[1/day]', 'Specific rate of oxidation of S0 with O2',          default=0.02_rk)
   call self%get_parameter(self%K_s2o3_ox,'K_s2o3_ox','[1/day]', 'Specific rate of oxidation of S2O3 with O2',        default=0.01_rk)
   call self%get_parameter(self%K_so4_rd, 'K_so4_rd', '[1/day]', 'Specific rate of OM sulfate reduction with sulfate',default=0.000005_rk)
   call self%get_parameter(self%K_s2o3_rd,'K_s2o3_rd','[1/day]', 'Specific rate of OM sulfate reduction with thiosulfate',default=0.001_rk)
   call self%get_parameter(self%K_s0_disp,'K_s0_disp','[1/day]', 'Specific rate of S0 dispropotionation',default=0.001_rk)
   call self%get_parameter(self%K_s0_no3,'K_s0_no3','[1/day]', 'Specific rate of oxidation of S0 with NO3',default=0.9_rk)
   call self%get_parameter(self%K_s2o3_no3,'K_s2o3_no3','[1/day]', 'Specific rate of oxidation of S2O3 with NO3',default=0.01_rk)
   call self%get_parameter(self%K_mnrd_hs,'K_mnrd_hs','[uM S]', 'half sat. of Mn reduction',default=1.0_rk)
   call self%get_parameter(self%K_ferd_hs,'K_ferd_hs','[uM S]', 'half sat. of Fe reduction',default=1.0_rk)
   !---- N---------!
   call self%get_parameter(self%K_DON_ox,'K_DON_ox','[1/day]', 'Specific rate of oxidation of DON with O2',   default=0.01_rk)
   call self%get_parameter(self%K_PON_ox,'K_PON_ox','[1/day]', 'Specific rate of oxidation of PON with O2',   default=0.002_rk)
   call self%get_parameter(self%Tda,     'Tda',     '[1/day]', 'Temperature control coefficient for OM decay',default=13.0_rk)
   call self%get_parameter(self%beta_da, 'beta_da', '[1/day]', 'Temperature control coefficient for OM decay',default=20.0_rk)
   call self%get_parameter(self%K_omox_o2, 'K_omox_o2', '[uM]','half sat. of o2 for OM mineralization',default=1.0_rk)
   call self%get_parameter(self%K_PON_DON, 'K_PON_DON', '[1/day]','Specific rate of Autolysis of PON to DON',default=0.1_rk)
   call self%get_parameter(self%K_nitrif1, 'K_nitrif1', '[1/day]','Spec.rate of 1st st. of nitrification',default=0.01_rk)
   call self%get_parameter(self%K_nitrif2, 'K_nitrif2', '[1/day]','Spec.rate of 2d st. of nitrification',default=0.1_rk)
   call self%get_parameter(self%K_denitr1, 'K_denitr1', '[1/day]','Spec.rate of 1 stage of denitrif',default=0.20_rk)
   call self%get_parameter(self%K_denitr2, 'K_denitr2', '[1/day]','Spec.rate of 2 stage of denitrif',default=0.25_rk)
   call self%get_parameter(self%K_omno_no3, 'K_omno_no3', '[uM N]','half sat. of no3 for OM denitr.',default=0.001_rk)
   call self%get_parameter(self%K_omno_no2, 'K_omno_no2', '[uM N]','half sat. of no2 for OM denitr.',default=0.001_rk)
   call self%get_parameter(self%K_hs_no3, 'K_hs_no3', '[1/day]','Spec.rate of thiodenitrification.',default=0.8_rk)
   call self%get_parameter(self%K_annamox, 'K_annamox', '[1/day]','Spec.rate of Anammox',default=0.8_rk)
  !---- O2--------!
   call self%get_parameter(self%O2s_nf, 'O2s_nf', '[uM O]','half saturation for nitrification',default=4.488_rk)
   call self%get_parameter(self%O2s_dn, 'O2s_dn', '[uM O]','half saturation for denitrification',default=10.0_rk)
   call self%get_parameter(self%s_omox_o2, 's_omox_o2', '[uM O]','threshold of o2 for OM mineralization',default=0.01_rk)
   call self%get_parameter(self%s_omno_o2, 's_omno_o2', '[uM O]','threshold of o2 for OM denitrification',default=25.0_rk)
   call self%get_parameter(self%s_omso_o2, 's_omso_o2', '[uM O]','threshold of o2 for OM sulfate reduction',default=25.0_rk)
   call self%get_parameter(self%s_omso_no3, 's_omso_no3', '[uM N]','threshold of noX for OM sulfate reduction',default=5.0_rk)
   call self%get_parameter(self%K_mnox_o2, 'K_mnox_o2', '[uM O]','half sat. of Mn oxidation',default=2.0_rk)
  !---- C--------!
   call self%get_parameter(self%K_caco3_diss, 'K_caco3_diss', '[1/day]','CaCO3 dissollution rate constant',default=3.0_rk)
   call self%get_parameter(self%K_caco3_form, 'K_caco3_form', '[1/day]','CaCO3 precipitation rate constant',default=0.0001_rk)
   call self%get_parameter(self%K_DON_ch4, 'K_DON_ch4', '[1/day]','Specific rate of methane production from DON',default=0.00014_rk)
   call self%get_parameter(self%K_PON_ch4, 'K_PON_ch4', '[1/day]','Specific rate of methane production from PON',default=0.00014_rk)
   call self%get_parameter(self%K_ch4_o2, 'K_ch4_o2', '[1/day]','Specific rate of oxidation of CH4 with O2',default=0.14_rk)
   call self%get_parameter(self%K_ch4_so4, 'K_ch4_so4', '[1/day]','Specific rate of anoxic oxidation of CH4 with SO4',default=0.0000274_rk)
   call self%get_parameter(self%s_omch_so4, 's_omch_so4', '[uM S]','threshold of SO4 for methane production from OM',default=15000.0_rk)
  !---- Si-------!
   call self%get_parameter(self%K_sipart_diss, 'K_sipart_diss', '[1/day]','Si dissollution rate constant',default=0.10_rk)
  !---- Bacteria-!
   call self%get_parameter(self%K_Baae_gro, 'K_Baae_gro', '[1/day]','Baae maximum specific growth rate',default=0.019_rk)
   call self%get_parameter(self%K_Baae_mrt, 'K_Baae_mrt', '[1/day]','Baae specific rate of mortality',default=0.005_rk)
   call self%get_parameter(self%K_Baae_mrt_h2s, 'K_Baae_mrt_h2s', '[1/day]','Baae increased specific rate of mortality due to H2S',default=0.899_rk)
   call self%get_parameter(self%limBaae, 'limBaae', '[1/day]','Limiting parameter for nutrient consumprion by Baae',default=2.0_rk)

   call self%get_parameter(self%K_Bhae_gro, 'K_Bhae_gro', '[1/day]','Bhae maximum specific growth rate',default=0.5_rk)
   call self%get_parameter(self%K_Bhae_mrt, 'K_Bhae_mrt', '[1/day]','Bhae specific rate of mortality',default=0.01_rk)
   call self%get_parameter(self%K_Bhae_mrt_h2s, 'K_Bhae_mrt_h2s', '[1/day]','Bhae increased specific rate of mortality due to H2S',default=0.799_rk)
   call self%get_parameter(self%limBhae, 'limBhae', '[1/day]','Limiting parameter for OM consumprion by Bhae',default=5.0_rk)

   call self%get_parameter(self%K_Baan_gro, 'K_Baan_gro', '[1/day]','Baan maximum specific growth rate',default=0.012_rk)
   call self%get_parameter(self%K_Baan_mrt, 'K_Baan_mrt', '[1/day]','Baan specific rate of mortality',default=0.012_rk)
   call self%get_parameter(self%limBaan, 'limBaan', '[1/day]','Limiting parameter for nutrient consumprion by Baan',default=2.0_rk)

   call self%get_parameter(self%K_Bhan_gro, 'K_Bhan_gro', '[1/day]','Bhan maximum specific growth rate',default=0.2_rk)
   call self%get_parameter(self%K_Bhan_mrt, 'K_Bhan_mrt', '[1/day]','Bhan specific rate of mortality',default=0.007_rk)
   call self%get_parameter(self%K_Bhan_mrt_o2, 'K_Bhan_mrt_o2', '[1/day]','Bhan increased specific rate of mortality due to O2',default=0.899_rk)
   call self%get_parameter(self%limBhan, 'limBhan', '[1/day]','Limiting parameter for OM consumprion by Bhan',default=2.0_rk)
  !---- Stoichiometric coefficients ----!
   call self%get_parameter(self%r_n_p,   'r_n_p',  '[-]',      'N[uM]/P[uM]',                  default=16.0_rk)
   call self%get_parameter(self%r_o_n,   'r_o_n',  '[-]',      'O[uM]/N[uM]',                  default=6.625_rk)
   call self%get_parameter(self%r_c_n,   'r_c_n',  '[-]',      'C[uM]/N[uM]',                  default=8.0_rk)
   call self%get_parameter(self%r_fe_n,   'r_fe_n',  '[-]',      'Fe[uM]/N[uM]',      default=26.5_rk)
   call self%get_parameter(self%r_mn_n,   'r_mn_n',  '[-]',      'Mn[uM]/N[uM]',      default=13.25_rk)
   call self%get_parameter(self%f,   'f',  '[-]',      'conversion factor relating solid and dissolved species concentrations', default=0.66_rk)
   call self%get_parameter(self%r_fe3_p,  'r_fe3_p',  '[-]',   'Fe[uM]/P[uM] partitioning coeff. for Fe oxide',   default=2.7_rk)
   call self%get_parameter(self%r_mn3_p,  'r_mn3_p',  '[-]',   'Mn[uM]/P[uM] partitioning coeff. for Mn(III)',     default=0.67_rk)
   call self%get_parameter(self%r_fe3_si,   'r_fe3_si',  '[-]',      'Fe[uM]/Si[uM] partitioning coeff. for Fe oxide',    default=2.7_rk)

   ! Register state variables
   call self%register_state_variable(self%id_Mn2, 'Mn2', 'mmol/m**3','Mn(II)', minimum=0.0_rk)
   call self%register_state_variable(self%id_Mn3, 'Mn3', 'mmol/m**3','Mn(III)',minimum=0.0_rk)
   call self%register_state_variable(self%id_Mn4, 'Mn4', 'mmol/m**3','Mn(IV)', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   call self%register_state_variable(self%id_H2S, 'H2S', 'mmol/m**3','H2S', minimum=0.0_rk)
   call self%register_state_variable(self%id_MnS, 'MnS', 'mmol/m**3','MnS', minimum=0.0_rk)
   call self%register_state_variable(self%id_MnCO3, 'MnCO3', 'mmol/m**3','MnCO3', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)

   call self%register_state_variable(self%id_Fe2, 'Fe2', 'mmol/m**3','Fe(II)', minimum=0.0_rk)
   call self%register_state_variable(self%id_Fe3, 'Fe3', 'mmol/m**3','Fe(III)', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   call self%register_state_variable(self%id_FeS, 'FeS', 'mmol/m**3','FeS', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   !!   call self%register_state_variable(self%id_FeS, 'FeS', 'mmol/m**3','iron sulfide', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   call self%register_state_variable(self%id_FeCO3, 'FeCO3', 'mmol/m**3','FeCO3', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)

   call self%register_state_variable(self%id_NO2, 'NO2', 'mmol/m**3','NO2', minimum=0.0_rk)
   call self%register_state_variable(self%id_S0 , 'S0',  'mmol/m**3','S0', minimum=0.0_rk)
   call self%register_state_variable(self%id_S2O3,'S2O3','mmol/m**3','S2O3', minimum=0.0_rk)
   call self%register_state_variable(self%id_SO4, 'SO4', 'mmol/m**3','SO4', minimum=0.0_rk)
   call self%register_state_variable(self%id_Si, 'Si', 'mmol/m**3','Si', minimum=0.0_rk)
   call self%register_state_variable(self%id_Sipart, 'Sipart', 'mmol/m**3','Si Particulate', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)

   call self%register_state_variable(self%id_Baae, 'Baae', 'mmol/m**3','Aerobic Autotrophic Bacteria', minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
   call self%register_state_variable(self%id_Bhae, 'Bhae', 'mmol/m**3','Aerobic Heterotrophic Bacteria', minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
   call self%register_state_variable(self%id_Baan, 'Baan', 'mmol/m**3','Anaerobic Autotrophic Bacteria', minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
   call self%register_state_variable(self%id_Bhan, 'Bhan', 'mmol/m**3','Anaerobic Heterotrophic Bacteria', minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)

   call self%register_state_variable(self%id_CaCO3, 'CaCO3', 'mmol/m**3','CaCO3', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_FeS2, 'FeS2', 'mmol/m**3','FeS2', minimum=0.0_rk, vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_CH4, 'CH4', 'mmol/m**3','CH4', minimum=0.0_rk)

   ! Register the contribution of all state variables to total nitrogen
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_bio)

   call self%register_dependency(self%id_Hplus, 'Hplus', 'mmol/m**3','H+')
   call self%register_dependency(self%id_Om_Ca,'Om_Ca','-','Omega CaCO3-Calcite')
   call self%register_dependency(self%id_Om_Ar,'Om_Ar','-','Omega CaCO3-Aragonite')
   call self%register_dependency(self%id_CO3,'CO3','mmol/m**3','CO3--')
   call self%register_dependency(self%id_Ca,'Ca','mmol/m**3','Ca++')
   call self%register_dependency(self%id_pco2,'pCO2','ppm','pCO2')
  !!   call self%register_dependency(self%id_Kh2s1,'Kh2s1','-', 'H2S <--> H+ + HS-')

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   call self%register_state_dependency(self%id_DIC,'DIC','mmol/m**3','total dissolved inorganic carbon',required=.false.)
   call self%register_state_dependency(self%id_Alk,'Alk','mmol/m**3','total alkalinity',required=.false.)
   call self%register_state_dependency(self%id_po4,'PO4','mmol/m**3','phosphate',required=.false.)

   call self%register_state_dependency(self%id_O2, 'O2', 'mmol/m**3','dissolved oxygen')
   call self%register_state_dependency(self%id_NH4,'NH4','mmol/m**3','ammonium')
   call self%register_state_dependency(self%id_NO3,'NO3','mmol/m**3','nitrate')
 !   call self%register_state_dependency(self%id_NO2,'NO2','mmol/m**3','nitrite')

   call self%register_state_dependency(self%id_PON,'PON','mmol/m**3','particulate organic nitrogen')
   call self%register_state_dependency(self%id_DON,'DON','mmol/m**3','dissolved organic nitrogen')

       !!!! Register diagnostic variables

   !call self%register_diagnostic_variable(self%id_GPP,'GPP','mmol/m**3',  'gross primary production',           &
   !                  output=output_time_step_integrated)

    call self%register_diagnostic_variable(self%id_DcPM_O2,'DcPM_O2','mmol/m**3',  'POM with O2 mineralization',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDM_O2,'DcDM_O2','mmol/m**3',  'DOM with O2 mineralization',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcPM_NOX,'DcPM_NOX','mmol/m**3',  'POM denitrification (1+2 stage)',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDM_NOX,'DcDM_NOX','mmol/m**3',  'DOM denitrification (1+2 stage)',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcPM_SO4,'DcPM_SO4','mmol/m**3',  'POM sulfatereduction (1+2 stage)',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDM_SO4,'DcDM_SO4','mmol/m**3',  'DOM sulfatereduction (1+2 stage)',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcPM_Mn,'DcPM_Mn','mmol/m**3',  'POM with Mn(IV) mineralization ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDM_Mn,'DcDM_Mn','mmol/m**3',  'DOM with Mn(IV) mineralization',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcPM_Fe,'DcPM_Fe','mmol/m**3',  'POM with Fe(III) mineralization',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDM_Fe,'DcDM_Fe','mmol/m**3',  'DOM with Fe(III) mineralization  ',           &
                output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_DcPM_CH4,'DcPM_CH4','mmol/m**3',  'CH4 production from PON',           &
    !            output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_DcDM_CH4,'DcDM_CH4','mmol/m**3',  'CH4 production from DON',           &
    !            output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_ChemBaae,'ChemBaae','mmol/m**3',  'Growth of  Aerobic Autotrophic Bacteria',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_HetBhae,'HetBhae','mmol/m**3',  'Growth of  Aerobic Heterotrophic Bacteria',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_autolysis,'Autolysis','mmol/m**3',  'Autolysis',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mn_ox1,'mn_ox1','mmol/m**3',  'Mn(II) with O2 oxidation ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mn_ox2,'mn_ox2','mmol/m**3',  'Mn(III) with O2 oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mns_diss,'mns_diss','mmol/m**3',  'MnS dissolution',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mns_form,'mns_form','mmol/m**3',  'MnS formation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mnco3_diss,'mnco3_diss','mmol/m**3',  'MnCO3 dissolusion',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mnco3_form,'mnco3_form','mmol/m**3',  'MnCO3 formation ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mn_rd1,'mn_rd1','mmol/m**3',  'Mn(IV) with H2S reduction',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mn_rd2,'mn_rd2','mmol/m**3',  'Mn(III) with H2S reduction',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_HetBhan,'HetBhan','mmol/m**3',  'Growth of Anaerobic Heterotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_ChemBaan,'ChemBaan','mmol/m**3',  'Growth of Anaerobic Autotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_MortBhan,'MortBhan','mmol/m**3',  'Mortality of Anaerobic Heterotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_MortBaan,'MortBaan','mmol/m**3',  'Mortality of Anaerobic Autotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_MortBaae,'MortBaae','mmol/m**3',  'Mortality of Aerobic Autotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_MortBhae,'MortBhae','mmol/m**3',  'Mortality of Aerobic Heterotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Nitrif1,'Nitrif1','mmol/m**3',  'Nitrification 1 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Nitrif2,'Nitrif2','mmol/m**3',  'Nitrification 2 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_anammox ,'Anammox ','mmol/m**3',  'Anammox',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr1_PM ,'Denitr1_PM ','mmol/m**3',  'POM denitrification 1 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr1_DM ,'Denitr1_DM ','mmol/m**3',  'DOM denitrification 1 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr2_PM ,'Denitr2_PM ','mmol/m**3',  'POM denitrification 2 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr2_DM ,'Denitr2_DM ','mmol/m**3',  'DOM denitrification 2 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr1 ,'Denitr1 ','mmol/m**3',  '(POM+DOM) denitrification 1 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr2 ,'Denitr2 ','mmol/m**3',  '(POM+DOM) denitrification 2 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self% id_fe_ox1,'fe_ox1','mmol/m**3',  'Fe(II) with O2 oxidation ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self% id_fe_ox2,'fe_ox2','mmol/m**3',  'Fe(II)  with Mn(IV) oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_fe_rd,'fe_rd','mmol/m**3',  'Fe (III) with H2S reduction',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_feco3_diss,'feco3_diss','mmol/m**3',  'FeCO3 dissolusion',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_feco3_form,'feco3_form','mmol/m**3',  'FeCO3 formation ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_s2o3_rd_PM,'s2o3_rd_PM','mmol/m**3',  'POM sulfatereduction 2d stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_s2o3_rd_DM,'s2o3_rd_DM','mmol/m**3',  'DOM sulfatereduction 2d stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_s2o3_no3,'s2o3_no3','mmol/m**3',  ' S2O3 with NO3 oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_s0_no3,'s0_no3','mmol/m**3',  'S0 with NO3 oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_s2o3_ox,'s2o3_ox','mmol/m**3',  'S2O3  with O2oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_s0_ox,'s0_ox','mmol/m**3',  'S0 with O2 oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_so4_rd,'so4_rd','mmol/m**3',  '(POM+DOM) sulfatereduction 1st stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_s2o3_rd,'s2o3_rd','mmol/m**3',  '(POM+DOM) sulfatereduction 2d stage ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_s0_disp,'s0_disp','mmol/m**3',  'S0 disproportionation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_hs_ox,'hs_ox','mmol/m**3',  'H2S with O2 oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_hs_no3,'hs_no3','mmol/m**3',  'H2S with NO3 oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_so4_rd_PM,'so4_rd_PM','mmol/m**3',  'POM sulfatereduction 1st stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_so4_rd_DM,'so4_rd_DM','mmol/m**3',  'DOM sulfatereduction 1st stage',           &
                output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_ch4_o2,'ch4_o2','mmol/m**3',  'CH4 with O2 oxidation',           &
    !            output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_ch4_so4,'ch4_so4','mmol/m**3',  'CH4 with SO4 oxidation',           &
    !            output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_fe_p_compl,'fe_p_compl','mmol/m**3',  'complexation of P with Fe(III)',  &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_mn_p_compl,'mn_p_compl','mmol/m**3',  'complexation of P with Mn(III)',  &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_fe_si_compl,'fe_si_compl','mmol/m**3',  'complexation of Si with Fe(III)', &
                output=output_time_step_integrated)
   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)

   ! Register dependencies on equilibrium constants for phosphoric acid, ammonia, hydrogen sulfide.
   call self%register_dependency(self%id_Kp1,  'Kp1',  '-', '[H+][H2PO4-]/[H3PO4]')
   call self%register_dependency(self%id_Kp2,  'Kp2',  '-', '[H][HPO4]/[H2PO4]')
   call self%register_dependency(self%id_Kp3,  'Kp3',  '-', '[H][PO4]/[HPO4]')
   call self%register_dependency(self%id_Knh4, 'Knh4', '-', '[H+][NH3]/[NH4]')
   call self%register_dependency(self%id_Kh2s,'Kh2s','-', '[H+][HS-]/[H2S]')
   call self%register_dependency(self%id_KSi,  'KSi','-','[H+][H3SiO4-]/[Si(OH)4]')
   call self%register_dependency(self%id_Kc0,  'Kc0','-','Henry''s constant')

   ! Specify that are rates computed in this module are per day (default: per second)
   self%dt = 86400.

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_redox),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk)                   :: O2,PON,DON,Si,Sipart
   real(rk)                   :: Mn2,Mn3,Mn4,H2S,MnS,MnCO3
   real(rk)                   :: Fe2,Fe3,FeS,FeS2,FeCO3
   real(rk)                   :: NO2,NO3,NH4, PO4
   real(rk)                   :: S0,S2O3,SO4
   real(rk)                   :: Baae,Bhae,Baan,Bhan
   real(rk)                   :: DIC,Alk,Hplus,CaCO3,Ca,CO3,Om_Ca,Om_Ar,pCO2,CH4
   real(rk)                   :: temp,salt

   real(rk) :: Autolysis,DcDM_O2,DcPM_O2, DcPM_NOX,DcDM_NOX
   real(rk) :: mn_ox1,mn_ox2,mn_rd1,mn_rd2,Om_MnS,mns_form,mns_diss,DcDM_Mn,DcPM_Mn,Om_MnCO3,mnco3_form,mnco3_diss,mnco3_ox
   real(rk) :: fe_ox1,fe_rd,fe_ox2,Om_FeS,fes_form,fes_diss,fes_ox,DcDM_Fe,DcPM_Fe,fes2_form,fes2_ox,Om_FeCO3,feco3_form,feco3_diss,feco3_ox
   real(rk) :: Nitrif1,Nitrif2,Anammox,Denitr1_PM,Denitr1_DM,Denitr2_PM,Denitr2_DM,Denitr1,Denitr2
   real(rk) :: s0_disp,hs_ox,s0_ox,s0_no3,s2o3_ox,s2o3_no3,hs_no3,so4_rd_PM,so4_rd_DM,s2o3_rd_PM,s2o3_rd_DM,so4_rd,s2o3_rd,DcPM_SO4,DcDM_SO4
   real(rk) :: fe_p_compl,mn_p_compl,fe_si_compl
   real(rk) :: ChemBaae,MortBaae,MortBhae,ChemBaan,MortBaan,MortBhan,HetBhan,HetBhae
   real(rk) :: Knh4,Kp1,Kp2,Kp3,Kh2s,KSi,Kc0
   real(rk) :: dAlk, Dc_OM_total,caco3_form,caco3_diss,DcDM_CH4,DcPM_CH4,ch4_o2,ch4_so4

!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_O2,O2)

   _GET_(self%id_PON,PON)
   _GET_(self%id_DON,DON)

   _GET_(self%id_Mn2,Mn2)
   _GET_(self%id_Mn3,Mn3)
   _GET_(self%id_Mn4,Mn4)
   _GET_(self%id_MnS,MnS)
   _GET_(self%id_MnCO3,MnCO3)
   _GET_(self%id_H2S,H2S)

   _GET_(self%id_Fe2,Fe2)
   _GET_(self%id_Fe3,Fe3)
   _GET_(self%id_FeS,FeS)
   _GET_(self%id_FeS2,FeS2)
   _GET_(self%id_FeCO3,FeCO3)

   _GET_(self%id_PO4,PO4)

   _GET_(self%id_NH4,NH4)
   _GET_(self%id_NO2,NO2)
   _GET_(self%id_NO3,NO3)

   _GET_(self%id_S0,  S0)
   _GET_(self%id_S2O3,S2O3)
   _GET_(self%id_SO4, SO4)

   _GET_(self%id_Baae,Baae)
   _GET_(self%id_Bhae,Bhae)
   _GET_(self%id_Baan,Baan)
   _GET_(self%id_Bhan,Bhan)

   _GET_(self%id_DIC,DIC)
   _GET_(self%id_Alk,Alk)
   _GET_(self%id_Hplus,Hplus)
   _GET_(self%id_CaCO3,CaCO3)
   _GET_(self%id_CH4,CH4)
   _GET_(self%id_CO3,CO3)
   _GET_(self%id_Om_Ca,Om_Ca)
   _GET_(self%id_Om_Ar,Om_Ar)
   _GET_(self%id_Ca,Ca)
   _GET_(self%id_pCO2,pCO2)

   ! Get equilibrium constants
   _GET_(self%id_Kp1,  Kp1)
   _GET_(self%id_Kp2,  Kp2)
   _GET_(self%id_Kp3,  Kp3)
   _GET_(self%id_Kh2s,Kh2s)
   _GET_(self%id_Knh4, Knh4)
   _GET_(self%id_KSi, KSi)
   _GET_(self%id_Kc0, Kc0)
   _GET_(self%id_Si,Si)
   _GET_(self%id_Sipart,Sipart)

   ! Retrieve current environmental conditions.
   !_GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_salt,salt)              ! salinity

!%%
!%!-------------------------------------------------------------------------
!%!========PON & DON========================================================
!%!-------------------------------------------------------------------------
     Autolysis=self%K_PON_DON*PON
!%(CH2O)106(NH3)16H3PO4 + 106 O2 -> 106 CO2 + 106 H2O + 16 NH3 + H3PO4 :
     DcDM_O2=self%K_DON_ox*DON*O2/(O2+self%K_omox_o2) &
     *(1.+self%beta_da*yy(self%tda,temp))
     !+++! *exp(0.15*t) !(Savchuk, Wulff,1996)
!%(CH2O)106(NH3)16H3PO4 + 106 O2 -> 106 CO2 + 106 H2O + 16 NH3 + H3PO4 :
     DcPM_O2=self%K_PON_ox*PON*O2/(O2+self%K_omox_o2) &
     *(1.+self%beta_da*yy(self%tda,temp))
!%!-------------------------------------------------------------------------
!%!========Mn===============================================================
!%!-------------------------------------------------------------------------
!% Mn2 oxidation: 4Mn2+ + O2 + 4H+ -> 4Mn3+ + 2H2O (Canfield, 2005):
     mn_ox1 =max(0._rk, 0.5*(1.+tanh(Mn2-self%s_mnox_mn2))*self%K_mn_ox1*Mn2*o2/(o2+self%K_mnox_o2) )
!% Mn3 oxidation: 2Mn3+ + 0.5O2 + 3H20 -> 2MnO2 + 6H+ (Tebo, 1997)	:
     mn_ox2=max(0._rk, 0.5*(1.+tanh(Mn3-self%s_mnox_mn3))*self%K_mn_ox2*Mn3*o2/(o2+self%K_mnox_o2) )
!% Mn4 reduction: 2MnO2 + 7H+ + HS- -> 2Mn3+ + 4H2O + S0 :
     mn_rd1 =max(0._rk, 0.5*(1.+tanh(Mn4-self%s_mnrd_mn4))*self%K_mn_rd1*Mn4*h2s/(h2s+self%K_mnrd_hs) )
!% Mn3 reduction: 2Mn3+ + HS- -> 2Mn2+ + S0 + H+ :
     mn_rd2=max(0._rk, 0.5*(1.+tanh(Mn3-self%s_mnrd_mn3))*self%K_mn_rd2*Mn3*h2s/(h2s+self%K_mnrd_hs) )

!MnS formation/dissollution (dSED) Mn2+  +  HS-  =  MnS(s)  +  H+
          Om_MnS=H2S*Mn2/(self%K_MnS*Hplus*1000000.)
!% Mn2+  +  HS-  ->  MnS(s)  +  H+
         mns_form=self%K_mns_form*max(0._rk,(Om_MnS-1._rk))
!%   MnS(s)  +  H+ -> Mn2+  +  HS-
         mns_diss=self%K_MnS_diss*MnS*max(0._rk,(1._rk-Om_MnS))
!% MnS + 2O2 ->   Mn2+ + SO42-  ?

!  MnCO3 precipitation/dissolution
 !   Om_MnCO3=Mn2*HCO3/(self%K_MnCO3*CO2)
        Om_MnCO3=Mn2*CO3/(self%K_MnCO3)
!% ! Mn2+ + CO3-- <-> MnCO3 (vanCappelen,96):
        mnco3_form=self%K_mnco3_form*max(0._rk,(Om_MnCO3-1._rk))
!%
        mnco3_diss=self%K_mnco3_diss*MnCO3*max(0._rk,(1._rk-Om_MnCO3))
 !% 2 MnCO3(s)  +  O2  +  2 H2O  =  2 MnO2(s)  +  2 HCO3-  +  2 H+ (Morgan,05)
        mnco3_ox=self%K_mnco3_ox*MnCO3*O2

!% (CH2O)106(NH3)16(H3PO4) +212MnO2 + 318CO2 +106H2O ->  424HCO3- + 212Mn2+ +16NH3 +H3PO4  (Boudreau, 1996)   ! in N units
     DcDM_Mn=max(0._rk,self%K_DON_mn*DON &
                        *Mn4/(Mn4 +0.5) &
                        *(1.-0.5*(1+tanh(o2-self%O2s_dn)))  )
!% (CH2O)106(NH3)16(H3PO4) +212MnO2 + 318CO2 +106H2O ->  424HCO3- + 212Mn2+ +16NH3 +H3PO4  (Boudreau, 1996)
     DcPM_Mn=max(0._rk,self%K_PON_mn*PON &
                        *Mn4/(Mn4 +0.5) &
                        *(1.-0.5*(1+tanh(o2-self%O2s_dn)))   )

!%!-------------------------------------------------------------------------
!%!========Fe===============================================================
!%!-------------------------------------------------------------------------
!% Fe2 oxidation1:4Fe2+ + O2 +10H2O-> 4Fe(OH)3 +8H+    (vanCappelen,96):
    fe_ox1= 0.5*(1.+tanh(Fe2-self%s_feox_fe2))*self%K_fe_ox1*o2*Fe2
!% Fe2 oxidation2: Fe2+ + MnO2 + 4H+ -> Fe3+ + Mn2+ + 2H2O (vanCappelen,96):
    fe_ox2=0.5*(1.+tanh(Fe2-self%s_feox_fe2))*self%K_fe_ox2*Mn4*Fe2
!% Fe3 reduction:  2Fe(OH)3+HS-+5H+ ->  2Fe2++S0+6H2O
    fe_rd=0.5*(1.+tanh(Fe3-self%s_ferd_fe3))*self%K_fe_rd*Fe3 *h2s/(h2s+self%K_ferd_hs) !*Kh2s1/(Kh2s1+Hplus)?

!FeS formation/dissollution (Bektursunova, 11)
          Om_FeS=H2S*Fe2/(self%K_fes*Hplus*1000000.)
!% FeS formation Fe2+ + HS- -> FeS + H+ (Bektursunova, 11)
    fes_form=self%K_fes_form*max(0._rk,(Om_FeS-1._rk))
!% FeS dissollution FeS + H+ -> Fe2+ + HS (Bektursunova, 11)
    fes_diss=self%K_fes_diss*FeS*max(0._rk,(1._rk-Om_FeS))
!% FeS oxidation: FeS + 2.25O2 +H2O -> 0.5Fe2O3 + 2H+ +SO42- (Soetaert,07) or FeS + 2O2 -> Fe2+ + SO42-(Bektursunova, 11) :
    fes_ox=self%K_fes_ox*FeS*O2

! Pyrite formation by FeS oxidation by H2S
! FeS + H2S -> FeS2 + H2     (Rickard,97)
      fes2_form=self%K_fes2_form*H2S*FeS
 !Pyrite oxidation by O2
 !  FeS2 + 3.5 O2 + H2O = Fe2+ + 2SO42- + 2H+  (Wijsman,02).
      fes2_ox=self%K_fes2_ox*FeS2*o2

      !  FeCO3 precipitation/dissolution
 !   Om_FeCO3=Fe2*HCO3/(self%K_FeCO3*CO2)
        Om_FeCO3=Fe2*CO3/(self%K_FeCO3)
!% ! Fe2+ + CO3-- <-> FeCO3 (vanCappelen,96?):
        feco3_form=self%K_feco3_form*max(0._rk,(Om_FeCO3-1._rk))
!%
        feco3_diss=self%K_FeCO3_diss*FeCO3*max(0._rk,(1._rk-Om_FeCO3))
 !%  FeCO3(s)  + O2  +  2 H2O  =   Fe2O3(s)  +   HCO3-  +   H+ (Morgan,05)=
        feco3_ox=self%K_feco3_ox*FeCO3*O2

!% (CH2O)106(NH3)16H3PO4 + 424Fe(OH)3 + 742CO2 -> 848HCO3-+ 424Fe2+ +318 H2O +16NH3 +H3PO4  (Boudreau, 1996)  Fe units
     DcDM_Fe=self%K_DON_fe*DON &
                        *Fe3 /(Fe3 +self%K_omno_no3) &
                        *(1.-0.5*(1+tanh(o2-self%O2s_dn)))
!% (CH2O)106(NH3)16H3PO4 + 424Fe(OH)3 + 742CO2 -> 848HCO3-+ 424Fe2+ +318 H2O +16NH3 +H3PO4  (Boudreau, 1996)
     DcPM_Fe=self%K_PON_fe*PON &
                        *Fe3 /(Fe3 +self%K_omno_no3) &
                        *(1.-0.5*(1+tanh(o2-self%O2s_dn)))
!%!-------------------------------------------------------------------------
!%!========N================================================================
!%!-------------------------------------------------------------------------
!%Nitrification 1st stage: NH4+ + 1.5 O2 -> NO2- + 2H+ + H2O (Canfield,05)
     Nitrif1 = self%K_nitrif1*NH4*o2*0.5*(1.+tanh(o2-self%O2s_nf))
!%Nitrification 2d stage: NO2- + 0.5 O2 -> NO3- (Canfield,05):
     Nitrif2 = self%K_nitrif2*NO2*o2*0.5*(1.+tanh(o2-self%O2s_nf))
!%!---------------------------------------- in suboxic conditions:
!%! Anammox : NO2- + NH4+ -> N2 + 2H2O (Canfield,05):
     Anammox = self%K_annamox*NO2*NH4 &
                      *(1.-0.5*(1+tanh(o2-self%O2s_dn)))
!%! OM denitrification (Richards, 1965): (CH2O)106(NH3)16H3PO4 + 84.8HNO3 = 106CO2 + 42.4N2 + 148.4H2O + 16NH3 + H3PO4
!%! POM denitrification (1st stage) (Anderson,82): 1/2CH2O + NO3- -> NO2- + 1/2H2O + 1/2CO2 :
     Denitr1_PM    = self%K_denitr1*PON &
                      *(1.-0.5*(1+tanh(o2-self%O2s_dn)))  &
                      *NO3/(NO3 +self%K_omno_no3)
!%!----------------------------------------
!%! DOM denitrification (1st stage): 1/2CH2O + NO3- -> NO2- + 1/2H2O + 1/2CO2 (Anderson,82):
     Denitr1_DM    = self%K_denitr1*DON &
                      *(1.-0.5*(1.+tanh(o2-self%O2s_dn)))  &
                      *NO3/(NO3+self%K_omno_no3)
!%!----------------------------------------
!%! POM denitrification (2d stage): 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 (Anderson,82):
     Denitr2_PM    = self%K_denitr2*PON &
                      *(1.-0.5*(1.+tanh(o2-self%O2s_dn)))  &
                      *NO2/(NO2+self%K_omno_no2)
!%!----------------------------------------
!%! DOM denitrification (2d stage): 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 (Anderson,82):
     Denitr2_DM    = self%K_denitr2*DON &
                      *(1.-0.5*(1+tanh(o2-self%O2s_dn)))  &
                      *NO2/(NO2+self%K_omno_no2)
     Denitr1       =Denitr1_PM + Denitr1_DM
     Denitr2       =Denitr2_PM + Denitr2_DM

!% ! from (Anderson, 1982) and (Richards, 1965): Denitrification POM (1+2 stage)
     DcPM_NOX      =16./212.*Denitr1_PM + 16./141.3*Denitr2_PM
   ! from (Anderson, 1982) and (Richards, 1965): Denitrification DOM (1+2 stage)
     DcDM_NOX      =16./212.*Denitr1_DM + 16./141.3*Denitr2_DM

!%!-------------------------------------------------------------------------
!%!========P================================================================
!%!-------------------------------------------------------------------------
!% complexation of P with Mn(III)
     mn_p_compl    =(mn_ox2+mn_rd2-mn_ox1-mn_rd1)/self%r_mn3_p
!% complexation of P with Fe(III)
     fe_p_compl    =(fe_rd-fe_ox1-fe_ox2+4.*DcDM_Fe+4.*DcPM_Fe)/self%r_fe3_p
!%%
!%!-------------------------------------------------------------------------
!%!========Si===============================================================
!%!-------------------------------------------------------------------------
!% complexation of Si with Fe(III)
     fe_si_compl   =(fe_rd-fe_ox1-fe_ox2+4.*DcDM_Fe+4.*DcPM_Fe)/self%r_fe3_si
!%%
!%!-------------------------------------------------------------------------
!%!========S================================================================
!%!-------------------------------------------------------------------------
!%! S0 disportionation: 4S0 + 3H2O -> 2H2S + S2O32= + 2H+ :
     s0_disp =self%K_s0_disp*S0
!%!----------------------------------------
!%! HS oxidation with O2: 2H2S + O2 -> 2S0 + 2H2O :
     hs_ox=self%K_hs_ox*o2*H2S
!%!----------------------------------------
!%! S0 oxidation with O2: 2S0 + O2 + H2O -> S2O32= + 2H+ :
     s0_ox=self%K_s0_ox*o2*S0
!%!----------------------------------------
!%! S0 oxidation with NO3: 4S0 + 3NO3- + 7H2O -> 4SO4= + 3NH4+ + 2H+  :
     s0_no3=self%K_s0_no3*NO3*S0
!%!----------------------------------------
!%! S2O3 oxidation with O2: S2O32= + 2O2 + 2OH- -> 2SO42= + H2O :
     s2o3_ox=self%K_s2o3_ox*o2*S2O3
!%!----------------------------------------
!%! S2O3 oxidation with NO3: S2O3= + NO3- + 2H2O --> 2SO4= + NH4+ :
     s2o3_no3=self%K_s2o3_no3*NO3*S2O3
!%!----------------------------------------
!%! Thiodenitrification: 3H2S + 4NO3- + 6OH- -> 3SO4= + 2N2 + 6H2O (Volkov, 1984):
     hs_no3       =self%K_hs_no3*H2S*NO3
!%!---------------------------------------- in anoxic conditions:
!%! OM sulfatereduction (Boudreau, 1996)
!% (CH2O)106(NH3)16H3PO4 + 53SO42- = 106HCO3- + 16NH3 + H3PO4 + 53H2S

!%! POM sulfatereduction (1st stage):
     so4_rd_PM      = (1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
                     *(1.-0.5*(1.+tanh(NO3-self%s_omso_no3)))  &
                      *self%K_so4_rd*SO4*PON
!%! DOM sulfatereduction (1st stage):
     so4_rd_DM      =(1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
                     *(1.-0.5*(1.+tanh(NO3-self%s_omso_no3)))  &
                      *self%K_so4_rd*SO4*DON

            if (o2.gt.10.) then
                        so4_rd_PM=0.
                        so4_rd_DM=0.
            endif

!%! POM sulfatereduction (2d stage):
     s2o3_rd_PM     =(1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
                     *(1.-0.5*(1.+tanh(NO3-self%s_omso_no3)))  &
                      *self%K_s2o3_rd*S2O3*PON
!%! DOM sulfatereduction (2d stage):
     s2o3_rd_DM     =(1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
                     *(1.-0.5*(1.+tanh(NO3-self%s_omso_no3)))  &
                      *self%K_s2o3_rd*S2O3*DON

           if (o2.gt.10.) then
                        s2o3_rd_PM=0.
                        s2o3_rd_DM=0.
           endif

     so4_rd         =so4_rd_PM  + so4_rd_DM  ! in S units
     s2o3_rd        =s2o3_rd_PM + s2o3_rd_DM

     DcPM_SO4      =16./53.*(so4_rd_PM+s2o3_rd_PM)  ! in N units
     DcDM_SO4      =16./53.*(so4_rd_DM+s2o3_rd_DM)
!%%
!%!-------------------------------------------------------------------------
!%!========C================================================================
!%!-------------------------------------------------------------------------
!%! CaCO3 precipitation/dissolution

     caco3_form =self%K_caco3_form*max(0._rk,(Om_Ar-1._rk))              ! (Luff et al., 2001)
     caco3_diss =CaCO3*self%K_caco3_diss*(max(0._rk,(1._rk-Om_Ar)))**4.5 ! (Luff et al., 2001)
!%!-------------------------------------------------------------------------
!%! CH4  production from PON and DON
!% (CH2O)106(NH3)16H3PO4 -> 53 CO2 + 53 CH4 + 16 NH3 + H3PO4 :
     DcDM_CH4=  (1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
              *(1.-0.5*(1.+tanh(NO3-self%s_omso_no3)))  &
              *(1.-0.5*(1.+tanh(SO4-self%s_omch_so4)))  &
              * self%K_DON_ch4*DON
!% (CH2O)106(NH3)16H3PO4 -> 53 CO2 + 53 CH4 + 16 NH3 + H3PO4 :
     DcPM_CH4= (1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
              *(1.-0.5*(1.+tanh(NO3-self%s_omso_no3)))  &
              *(1.-0.5*(1.+tanh(SO4-self%s_omch_so4)))  &
              * self%K_PON_ch4*PON
!%!-------------------------------------------------------------------------
!%! CH4  oxidation with O2
!% CH4 + 2 O2 = CO2 + 2 H2O
     ch4_o2= self%K_ch4_o2*CH4*O2
!%! CH4 anoxic oxidation with SO4
!% CH4 +  SO42- + 2 H+  =  CO2 + H2S + 2 H2O
     ch4_so4= self%K_ch4_so4*CH4*SO4
!%!---------------------------------------------------------------------------
!%!========B=a=c=t(ci)============================================================
!%!---------------------------------------------------------------------------

!%!---------------------------------------- OXIC CONDITIONS
!%!---------------------------------------- aerobic autotrophs

     ChemBaae =(Nitrif1+Nitrif2+mn_ox1+fe_ox1+s2o3_ox+s0_ox +Anammox) &
              *self%K_Baae_gro*Baae*min(yy(self%limBaae,NH4/(Baae+0.0001)),yy(self%limBaae,PO4/(Baae+0.0001)))

     MortBaae=(self%K_Baae_mrt+ self%K_Baae_mrt_h2s*(0.5*(1.-tanh(1.-H2S))))*Baae*Baae
!%!---------------------------------------- aerobic heterotroph

     HetBhae  =(DcPM_O2+DcDM_O2) &
             *self%K_Bhae_gro*Bhae*yy(self%limBhae,DON/(Bhae+0.0001))

     MortBhae =(self%K_Bhae_mrt+ self%K_Bhae_mrt_h2s*(0.5*(1.-tanh(1.-H2S))))*Bhae


!%!---------------------------------------- ANOXIC CONDITIONS
!%!------------------------------------- anaerobic autotrophs

     ChemBaan =(mn_rd1+mn_rd2+fe_rd+hs_ox +hs_no3) &
              *self%K_Baan_gro*Baan*min(yy(self%limBaan,NH4/(Baan+0.0001)),yy(self%limBaan,PO4/(Baan+0.0001)))

     MortBaan = self%K_Baan_mrt*Baan*Baan

!%!---------------------------------------- anaerobic heterotroph
     HetBhan  =(DcPM_NOX+DcDM_NOX+DcPM_SO4+DcDM_SO4+DcDM_Mn+DcPM_Mn+DcPM_Fe+DcDM_Fe+DcPM_ch4) &
                    *self%K_Bhan_gro*Bhan*yy(self%limBhan,DON/(Bhan+0.0001))

     MortBhan =(self%K_Bhan_mrt+ self%K_Bhan_mrt_o2*(0.5+0.5*(tanh(1.-O2))))*Bhan

!%!---------------------------------------- Summariazed OM mineralization
     Dc_OM_total=DcDM_O2+DcPM_O2+DcPM_NOX+DcDM_NOX+DcDM_Mn+DcPM_Mn+DcDM_Fe+DcPM_Fe &
                 +DcDM_SO4+DcPM_SO4+0.5*(DcDM_ch4+DcPM_ch4)


!%!----------------------------------------

          _ADD_SOURCE_(self%id_Mn2,-mn_ox1+mn_rd2-mns_form+mns_diss-mnco3_form+mnco3_diss+0.5*fe_ox2+(DcDM_Mn+DcPM_Mn)*self%r_mn_n)
          _ADD_SOURCE_(self%id_Mn3, mn_ox1-mn_ox2+mn_rd1-mn_rd2)
          _ADD_SOURCE_(self%id_Mn4, mn_ox2-mn_rd1-0.5*fe_ox2+mnco3_ox-(DcDM_Mn+DcPM_Mn)*self%r_mn_n)
          _ADD_SOURCE_(self%id_MnS, mns_form-mns_diss)
          _ADD_SOURCE_(self%id_MnCO3, mnco3_form-mnco3_diss-mnco3_ox)

          _ADD_SOURCE_(self%id_Fe2,-fe_ox1-fe_ox2+fe_rd-fes_form+fes_diss-feco3_form+feco3_diss+(DcDM_Fe+DcPM_Fe)*4.*self%r_fe_n+feS2_ox)
          _ADD_SOURCE_(self%id_Fe3,fe_ox1+fe_ox2-fe_rd+fes_ox+feco3_ox-(DcDM_Fe+DcPM_Fe)*4.*self%r_fe_n)
          _ADD_SOURCE_(self%id_FeS,fes_form-fes_diss-fes_ox-feS2_form)
          _ADD_SOURCE_(self%id_FeS2,feS2_form-feS2_ox)
          _ADD_SOURCE_(self%id_FeCO3, feco3_form-feco3_diss-feco3_ox)

          _ADD_SOURCE_(self%id_H2S,-0.5*mn_rd1-0.5*mn_rd2-0.5*fe_rd -hs_ox-fes_form-mns_form+mns_diss +0.5*s0_disp -hs_no3+s2o3_rd+fes_diss-feS2_form)
          _ADD_SOURCE_(self%id_S0, hs_ox+0.5*mn_rd1+0.5*mn_rd2+0.5*fe_rd -s0_ox -s0_disp-s0_no3)
          _ADD_SOURCE_(self%id_S2O3,0.5*s0_ox-s2o3_ox+0.25*s0_disp +0.5*so4_rd-0.5*s2o3_rd-s2o3_no3)
          _ADD_SOURCE_(self%id_SO4,hs_no3-so4_rd+0.5*s2o3_ox+s0_no3+2.*s2o3_no3+fes_ox+2.*feS2_ox-ch4_so4)

          _ADD_SOURCE_(self%id_O2,(-DcDM_O2-DcPM_O2)*self%r_o_n -0.25*mn_ox1-0.25*mn_ox2-0.25*fe_ox1 -0.5*hs_ox-0.5*s0_ox-0.5*s2o3_ox -1.5*Nitrif1-0.5*Nitrif2-2.25*fes_ox-3.5*feS2_ox-0.5*mnco3_ox+feco3_ox-2.*ch4_o2)
          _ADD_SOURCE_(self%id_DON,(Autolysis-DcDM_O2-DcDM_NOX-DcDM_SO4 -DcDM_Mn -DcDM_Fe -HetBhae-HetBhan))
          _ADD_SOURCE_(self%id_PON,(-Autolysis-DcPM_O2-DcPM_NOX-DcPM_SO4 -DcPM_Mn -DcPM_Fe+MortBaae+MortBaan+MortBhae+MortBhan))
          _ADD_SOURCE_(self%id_NH4,(Dc_OM_total-Nitrif1-Anammox +0.75*s0_no3+s2o3_no3-ChemBaae-ChemBaan))
          _ADD_SOURCE_(self%id_NO2,(Nitrif1-Nitrif2+Denitr1-Denitr2-Anammox))
          _ADD_SOURCE_(self%id_NO3,(Nitrif2-Denitr1-1.6*hs_no3-0.75*s0_no3-s2o3_no3) )
          _ADD_SOURCE_(self%id_DIC,(Dc_OM_total-ChemBaae-ChemBaan)*self%r_c_n-caco3_form+caco3_diss-mnco3_form+mnco3_diss+mnco3_ox -feco3_form+feco3_diss+feco3_ox)
          _ADD_SOURCE_(self%id_Si,self%K_sipart_diss*Sipart+fe_si_compl)
          _ADD_SOURCE_(self%id_Sipart,-self%K_sipart_diss*Sipart)
          _ADD_SOURCE_(self%id_PO4,((Dc_OM_total-ChemBaae-ChemBaan)/self%r_n_p+fe_p_compl+mn_p_compl))
          _ADD_SOURCE_(self%id_Baae,ChemBaae-MortBaae)
          _ADD_SOURCE_(self%id_Baan,ChemBaan-MortBaan)
          _ADD_SOURCE_(self%id_Bhae,HetBhae-MortBhae)
          _ADD_SOURCE_(self%id_Bhan,HetBhan-MortBhan)
          _ADD_SOURCE_(self%id_caco3,caco3_form-caco3_diss)
          _ADD_SOURCE_(self%id_CH4,0.5*(DcDM_CH4+DcPM_CH4)-ch4_o2-ch4_so4)

  !! Alkalinity changes due to redox reactions:

         dAlk = (  & ! 0.01*dCc(k,i_DON) & ! part of OM charged (EPS etc.)
                    !+( Dc_OM_total-ChemBaae-ChemBaan) *Knh4/(Knh4+Hplus)  &  ! +/- NH3(tot) &
                    -(ChemBaae-ChemBaan) & ! +/- NH3(tot) &
                    !+(-Anammox  & !NO2- + NH4+ -> N2 + 2H2O                       ! +/- NH3(tot)
                    !  +0.75*s0_no3 & ! 4S0 + 3NO3- + 7H2O -> 4SO4-- + 3NH4+ + 2H+ ! +/- NH3(tot)
                    !  +s2o3_no3 &  ! S2O3-- + NO3- + 2H2O -> 2SO4-- + NH4+         ! +/- NH3(tot)
                    ! ) *Knh4/(Knh4+Hplus) &                                       ! +/- NH3(tot)
                    - 1.* Nitrif1 & !*(1.-Knh4/(Knh4+Hplus))  &    ! NH4+ + 1.5 O2 -> NO2- + 2H+ + H2O ! W-Gladrow,07   ! +/- HNO2t
!                    - 2.* Nitrif2 & ! NO2- + 0.5 O2 -> NO3  !                                                          ! +/- HNO2t
!: 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 or! 5 CH2O + 4 H+ + 4 NO3- -> 2 N2 + 5 CO2 + 7H2O ! W-Gladrow,07
!                    + 0.*(Denitr1_PM +Denitr1_DM)  &! &    ! +/- HNO2t
                    + 1.*(Denitr2_PM +Denitr2_DM)  &! &    ! +/- HNO2t
                    !+((Dc_OM_total-ChemBaae-ChemBaan)/self%r_n_p &                                            ! +/- PO4t
                    !    +(fe_rd-fe_ox1-fe_ox2+4.*DcDM_Fe+4.*DcPM_Fe)/2.7-(mn_ox2+mn_rd2-mn_ox1-mn_rd1)/0.67) & ! +/- PO4t
                    !    *((Kp1*Kp2-Hplus*Hplus)*Hplus+2e0*Kp1*Kp2*Kp3)/(((Hplus+Kp1) &                      ! +/- PO4t
                    !    *Hplus+Kp1*Kp2)*Hplus+Kp1*Kp2*Kp3)  &                                               ! +/- PO4t
                    !+(-0.5*mn_rd2-0.5*mn_rd2-0.5*fe_rd -hs_ox-fes_form+fes_diss-mns_form  &           ! +/- H2S (tot)
                    !    +0.5*s0_disp -hs_no3 +s2o3_rd-feS2_form)*Kh2s1/(Kh2s1+Hplus)  &              ! +/- H2S (tot)
                    + 2.*(so4_rd + s2o3_rd) &      ! (CH2O)106(NH3)16H3PO4 + 53SO42- = 106HCO3- + 16NH3 + H3PO4 + 53H2S (Boudreau, 1996)
                    !+(self%K_sipart_diss*Sipart)*KSi/(KSi+Hplus) &      !  ! +/- Si (tot)

                    +     mn_ox1   &  ! 4Mn2+ + O2 + 4H+ -> 4Mn3+ + 2H2O
                    - 3.*mn_ox2  &   ! 2Mn3+ + 3H2O  + 0.5 O2 -> 2MnO2 + 6H+
                    + 3.*mn_rd1   &   ! 2MnO2 + 7H+ + HS- -> 2Mn3+ + 4H2O + S0
                    - 1.*mn_rd2  &   ! 2Mn3+ + HS- -> 2Mn2+ + S0 + H+
                    - 2.*mns_form  & ! Mn2+ + H2S <-> MnS + 2H+                                           !!!!
                    + 2.*mns_diss  &  !
                    - 2.*mnco3_form & ! Mn2+ + CO3-- <-> MnCO3 (vanCappelen,96):
                    + 2.*mnco3_diss &
                    + 26.5*(DcDM_Mn +DcPM_Mn ) & ! & !DcDM_Mn is in N-units, i.e. 424/16
!!% (CH2O)106(NH3)16(H3PO4) +212MnO2 + 318CO2 +106H2O ->  424HCO3- + 212Mn2+ +16NH3 +H3PO4  (Boudreau, 1996)
                    - 2.*fe_ox1  &    ! 4Fe2+ + O2 +10H2O-> 4Fe(OH)3 +8H+    (vanCappelen,96):
                    - 1.*fe_ox2  &   ! 2Fe2+ + MnO2 +4H2O -> 2Fe(OH)3 + Mn2+ +2H+  (vanCappelen,96):                                          !!!!
                    + 2.*fe_rd  &   ! 2Fe(OH)3 + HS- + 5H+ ->  2Fe2++S0+6H2O  (here and below d(AlK_H2S) is excluded, as give before)                                        !!!!
                    - 1.*fes_form  & ! Fe2+ + H2S  <-> FeS + H+ (Bektursunova, 11)
                    + 1.*fes_diss  & !
                    - 2.*fes_ox &  ! FeS + 2.25O2 +H2O -> 0.5Fe2O3 + 2H+ +SO42- (Soetaert,07)
                                   ! FeS + H2S -> FeS2 + H2     (Rickard, 1997,Meysmann,2003,Soetaert,07)
                    - 2.*fes2_ox & !  FeS2 + 3.5 O2 + H2O -> Fe2+ + 2 SO42- + 2 H+ (Wijsman,02)
                    - 2.*feco3_form & ! Fe2+ + CO3-- <-> FeCO3 (vanCappelen,96):
                    + 2.*feco3_diss &
                    + 53.*(DcDM_Fe +DcPM_Fe) & ! !DcDM_Fe is in N-units, i.e. 848/16
!% (CH2O)106(NH3)16H3PO4 + 424Fe(OH)3 + 742CO2 -> 848HCO3-+ 424Fe2+ +318 H2O +16NH3 +H3PO4  (Boudreau, 1996)
                    - 0.5*s0_disp  & ! 4S0 + 3H2O -> 2H2S + S2O3-- + 2H+                                         !!!!
                    - 1.*(-s0_ox)  & ! 2S0 + O2 + H2O -> S2O3-- + 2H+
                    - 0.5*s0_no3  &  ! 4S0 + 3NO3- + 7H2O -> 4SO4-- + 3NH4+ + 2H+
                    - 1.*s2o3_ox  &   ! S2O3-- + 2O2 + 2OH- -> 2SO4-- + H2O
                    - 0.4*hs_no3 &  ! 5H2S + 8NO3- + 2OH+ -> 5SO4-- + 4N2 + 6H2O (Volkov, 1984)
                    - 2.*caco3_form &   ! Ca2+ + CO32- -> CaCO3                                        !!!!
                    + 2.*caco3_diss &   ! CaCO3 -> Ca2+ + CO32-
!                    + 0.15*(Autolysis-DcDM_O2-DcDM_NOX-DcDM_SO4 -DcDM_Mn -DcDM_Fe -HetBhae-HetBhan)*7. & ! part of OM charged (EPS etc.)
                    )
    _ADD_SOURCE_(self%id_Alk,dAlk)

    _SET_DIAGNOSTIC_(self%id_DcPM_O2,DcPM_O2)
    _SET_DIAGNOSTIC_(self%id_DcDM_O2,DcDM_O2)
    _SET_DIAGNOSTIC_(self%id_DcPM_NOX,DcPM_NOX)
    _SET_DIAGNOSTIC_(self%id_DcDM_NOX,DcDM_NOX)
    _SET_DIAGNOSTIC_(self%id_DcDM_SO4,DcPM_SO4)
    _SET_DIAGNOSTIC_(self%id_DcPM_SO4,DcPM_SO4)
    _SET_DIAGNOSTIC_(self%id_DcPM_Mn,DcPM_Mn)
    _SET_DIAGNOSTIC_(self%id_DcDM_Mn,DcDM_Mn)
    _SET_DIAGNOSTIC_(self%id_DcDM_Fe,DcPM_Fe)
    _SET_DIAGNOSTIC_(self%id_DcPM_Fe,DcPM_Fe)
    _SET_DIAGNOSTIC_(self%id_ChemBaae,ChemBaae)
    _SET_DIAGNOSTIC_(self%id_HetBhae,HetBhae)
    _SET_DIAGNOSTIC_(self%id_ChemBaan,ChemBaan)
    _SET_DIAGNOSTIC_(self%id_HetBhan,HetBhan)
    _SET_DIAGNOSTIC_(self%id_MortBaan,MortBaan)
    _SET_DIAGNOSTIC_(self%id_MortBhan,MortBhan)
    _SET_DIAGNOSTIC_(self%id_MortBaae,MortBaae)
    _SET_DIAGNOSTIC_(self%id_MortBhae,MortBhae)
    _SET_DIAGNOSTIC_(self%id_autolysis,Autolysis)
    _SET_DIAGNOSTIC_(self%id_anammox ,Anammox)
    _SET_DIAGNOSTIC_(self%id_Denitr1_PM ,Denitr1_PM)
    _SET_DIAGNOSTIC_(self%id_Denitr1_DM ,Denitr1_DM)
    _SET_DIAGNOSTIC_(self%id_Denitr2_PM ,Denitr2_PM)
    _SET_DIAGNOSTIC_(self%id_Denitr2_DM ,Denitr2_DM)
    _SET_DIAGNOSTIC_(self%id_Denitr1 ,Denitr1)
    _SET_DIAGNOSTIC_(self%id_Denitr2 ,Denitr2)
    _SET_DIAGNOSTIC_(self%id_fe_ox1,fe_ox1)
    _SET_DIAGNOSTIC_(self%id_fe_rd,fe_rd)
    _SET_DIAGNOSTIC_(self%id_s2o3_no3,s2o3_no3)
    _SET_DIAGNOSTIC_(self%id_s0_no3,s0_no3)
    _SET_DIAGNOSTIC_(self%id_s2o3_ox,s2o3_ox)
    _SET_DIAGNOSTIC_(self%id_s0_ox,s0_ox)
    _SET_DIAGNOSTIC_(self%id_s2o3_rd_PM,s2o3_rd_PM)
    _SET_DIAGNOSTIC_(self%id_s2o3_rd_DM,s2o3_rd_DM)
    _SET_DIAGNOSTIC_(self%id_so4_rd_PM,so4_rd_PM)
    _SET_DIAGNOSTIC_(self%id_so4_rd_DM,so4_rd_DM)
    _SET_DIAGNOSTIC_(self%id_so4_rd,so4_rd)
    _SET_DIAGNOSTIC_(self%id_s2o3_rd,s2o3_rd)
    _SET_DIAGNOSTIC_(self%id_s0_disp,s0_disp)
    _SET_DIAGNOSTIC_(self%id_hs_ox,hs_ox)
    _SET_DIAGNOSTIC_(self%id_hs_no3,hs_no3)
    _SET_DIAGNOSTIC_(self%id_fe_ox2,fe_ox2)
    _SET_DIAGNOSTIC_(self%id_mn_ox1,mn_ox1)
    _SET_DIAGNOSTIC_(self%id_mn_ox2,mn_ox2)
    _SET_DIAGNOSTIC_(self%id_mn_rd1,mn_rd1)
    _SET_DIAGNOSTIC_(self%id_mn_rd2,mn_rd2)
    _SET_DIAGNOSTIC_(self%id_Nitrif1,Nitrif1)
    _SET_DIAGNOSTIC_(self%id_Nitrif2,Nitrif2)
    _SET_DIAGNOSTIC_(self%id_mns_diss,mns_diss)
    !_SET_DIAGNOSTIC_(self%id_FeS_ox,FeS_ox)
    !_SET_DIAGNOSTIC_(self%id_fes_form,fes_form)
    !!_SET_DIAGNOSTIC_(self%id_mns_ox,mns_ox)
    _SET_DIAGNOSTIC_(self%id_mns_form,mns_form)
    !_SET_DIAGNOSTIC_(self%id_caco3_diss,caco3_diss)
    _SET_DIAGNOSTIC_(self%id_mnco3_diss,mnco3_diss)
    _SET_DIAGNOSTIC_(self%id_mnco3_form,mnco3_form)
    _SET_DIAGNOSTIC_(self%id_feco3_diss,feco3_diss)
    _SET_DIAGNOSTIC_(self%id_feco3_form,feco3_form)
    _SET_DIAGNOSTIC_(self%id_fe_p_compl,fe_p_compl)
    _SET_DIAGNOSTIC_(self%id_mn_p_compl,mn_p_compl)
    _SET_DIAGNOSTIC_(self%id_fe_si_compl,fe_si_compl)
    _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!----------------------------------------------------------------------
! SUBROUTINE: Saturation function squared
!-----------------------------------------------------------------------
!
! !INTERFACE:
   real(rk) function yy(a,x)
!
! !DESCRIPTION:
! This is a squared Michaelis-Menten type of limiter:
! \begin{equation}\label{Y}
! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
! \end{equation}
!
! !INPUT PARAMETERS:
!  real(rk), intent(in)
!  real :: a
   real(rk)  :: a,x
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   yy=x**2/(a**2+x**2)

   end function yy
!EOC

   end module fabm_niva_brom_redox
