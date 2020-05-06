#include "fabm_driver.h"

module csiro_seagrass

   ! This is a FABM implementation of the CSIRO seagrass model described in
   ! Baird et al. (2016) A biophysical representation of seagrass growth for application in a complex shallow-water biogeochemical model. Ecol. Modell. 325, 13–27.
   ! http://dx.doi.org/10.1016/j.ecolmodel.2015.12.011
   !
   ! Original authors: Matthew Adams, Jorn Bruggeman

   use fabm_types
   use fabm_particle

   implicit none

   private

   type,extends(type_particle_model),public :: type_csiro_seagrass
      ! State variables
      type (type_bottom_state_variable_id) :: id_SG_A
      type (type_bottom_state_variable_id) :: id_SG_B

      ! Diagnostics
      type (type_bottom_diagnostic_variable_id) :: id_E_par_below

      ! External sinks/sources
      type (type_state_variable_id)        :: id_DIC
      type (type_state_variable_id)        :: id_O2
      type (type_state_variable_id)        :: id_D_A_C, id_D_A_N, id_D_A_P
      type (type_bottom_state_variable_id) :: id_D_B_C, id_D_B_N, id_D_B_P
      type (type_bottom_state_variable_id) :: id_N_source, id_P_source

      ! Environmental dependencies
      type (type_bottom_dependency_id) :: id_N, id_P
      type (type_bottom_dependency_id) :: id_E_par
      type (type_dependency_id)        :: id_T

      ! Parameters
      real(rk) :: mu_max_SG
      real(rk) :: Omega_SG
      real(rk) :: A_par
      real(rk) :: f_below
      real(rk) :: tau_tran
      real(rk) :: K_SG_P
      real(rk) :: K_SG_N
      real(rk) :: E_comp
      real(rk) :: zeta_SG_A
      real(rk) :: zeta_SG_B
      real(rk) :: f_seed
      real(rk) :: z_root
      real(rk) :: sin_beta_blade
      real(rk) :: Q_10
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   ! Molar masses
   real(rk), parameter :: M_N = 14.007_rk
   real(rk), parameter :: M_C = 12.011_rk
   real(rk), parameter :: M_P = 30.974_rk

   ! Stoichiometry (Atkinson ratio)
   real(rk), parameter :: N_to_P = 30._rk
   real(rk), parameter :: C_to_P = 550._rk

contains

   subroutine initialize(self, configunit)
      class (type_csiro_seagrass), intent(inout), target :: self
      integer,                     intent(in)            :: configunit

      ! Parameters
      call self%get_parameter(self%mu_max_SG,      'mu_max_SG',     'd-1',                'maximum growth rate of above-ground seagrass',default=0.4_rk)
      call self%get_parameter(self%Omega_SG,       'Omega_SG',      '(g N m-2)-1',        'nitrogen-specific area of seagrass',          default=1.5_rk)
      call self%get_parameter(self%A_par,          'A_par',         '-',                  'leaf absorbance',                             default=0.7_rk)
      call self%get_parameter(self%f_below,        'f_below',       '-',                  'fraction biomass below ground',               default=0.5_rk)
      call self%get_parameter(self%tau_tran,       'tau_tran',      'd-1',                'translocation rate',                          default=0.033_rk)
      call self%get_parameter(self%K_SG_P,         'K_SG_P',        'g P m-3',            'half-saturation P uptake',                    default=0.096_rk, scale_factor=1000._rk/M_N)
      call self%get_parameter(self%K_SG_N,         'K_SG_N',        'g N m-3',            'half-saturation N uptake',                    default=0.420_rk, scale_factor=1000._rk/M_P)
      call self%get_parameter(self%E_comp,         'E_comp',        'mol photons m-2 d-1','compensation scalar PAR irradiance',          default=4.5_rk)
      call self%get_parameter(self%zeta_SG_A,      'zeta_SG_A',     'd-1',                'leaf loss rate',                              default=0.03_rk)
      call self%get_parameter(self%zeta_SG_B,      'zeta_SG_B',     'd-1',                'root loss rate',                              default=0.004_rk)
      call self%get_parameter(self%f_seed,         'f_seed',        '-',                  'seed biomass as a fraction of 63 % cover',    default=0.01_rk)
      call self%get_parameter(self%z_root,         'z_root',        'm',                  'seagrass root depth',                         default=0.15_rk)
      call self%get_parameter(self%sin_beta_blade, 'sin_beta_blade','-',                  'sine of the nadir blade angle',               default=0.5_rk)
      call self%get_parameter(self%Q_10,           'Q_10',          '-',                  'Q_10 factor',                                 default=2.0_rk)

      ! Internal state variables
      call self%register_state_variable(self%id_SG_A,'SG_A','g N m-2','above-ground seagrass biomass',initial_value=self%f_seed/self%Omega_SG*(1.0_rk-self%f_below))
      call self%register_state_variable(self%id_SG_B,'SG_B','g N m-2','below-ground seagrass biomass',initial_value=self%f_seed/self%Omega_SG*self%f_below)

      ! Register contributions to mass budgets (all in mmol rather than g, hence multiply by 1000 and divide by atomic mass of nitrogen)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_SG_A, scale_factor=1000._rk/M_N)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,    self%id_SG_A, scale_factor=1000._rk/M_N/N_to_P*C_to_P)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_SG_A, scale_factor=1000._rk/M_N/N_to_P)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_SG_B, scale_factor=1000._rk/M_N)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,    self%id_SG_B, scale_factor=1000._rk/M_N/N_to_P*C_to_P)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_SG_B, scale_factor=1000._rk/M_N/N_to_P)

      ! Register dependencies on external state variables
      call self%register_state_dependency(self%id_DIC,  'DIC',  'mmol C m-3', 'dissolved inorganic carbon')
      call self%register_state_dependency(self%id_O2,   'O2',   'mmol O2 m-3','oxygen')
      call self%register_state_dependency(self%id_N_source,'N_source', 'mmol N m-2', 'dissolved inorganic nitrogen source')
      call self%register_state_dependency(self%id_P_source,'P_source', 'mmol P m-2', 'dissolved inorganic phosphorus source')
      call self%register_state_dependency(self%id_D_A_C,'D_A_C','mmol C m-3', 'sink for carbon in leaf detritus')
      call self%register_state_dependency(self%id_D_A_N,'D_A_N','mmol N m-3', 'sink for nitrogen in leaf detritus')
      call self%register_state_dependency(self%id_D_A_P,'D_A_P','mmol P m-3', 'sink for phosphorus in leaf detritus')
      call self%register_state_dependency(self%id_D_B_C,'D_B_C','mmol C m-2', 'sink for carbon in root detritus')
      call self%register_state_dependency(self%id_D_B_N,'D_B_N','mmol N m-2', 'sink for nitrogen in root detritus')
      call self%register_state_dependency(self%id_D_B_P,'D_B_P','mmol P m-2', 'sink for phosphorus in root detritus')

      call self%request_coupling_to_model(self%id_D_A_C,'D_A',standard_variables%total_carbon)
      call self%request_coupling_to_model(self%id_D_A_N,'D_A',standard_variables%total_nitrogen)
      call self%request_coupling_to_model(self%id_D_A_P,'D_A',standard_variables%total_phosphorus)
      call self%request_coupling_to_model(self%id_D_B_C,'D_B',standard_variables%total_carbon)
      call self%request_coupling_to_model(self%id_D_B_N,'D_B',standard_variables%total_nitrogen)
      call self%request_coupling_to_model(self%id_D_B_P,'D_B',standard_variables%total_phosphorus)

      ! Register environmental dependencies
      call self%register_dependency(self%id_T,     standard_variables%temperature)
      call self%register_dependency(self%id_E_par, 'E_par', 'W m-2',      'photosynthetically active radiation at top of canopy')
      call self%register_dependency(self%id_N,     'N',     'mmol N m-3', 'dissolved inorganic nitrogen in porewater')
      call self%register_dependency(self%id_P,     'P',     'mmol P m-3', 'dissolved inorganic phosphorus in porewater')

      ! Register diagnostics
      call self%register_diagnostic_variable(self%id_E_par_below,'E_par_below','W m-2','photosynthetically active radiation below canopy')
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_csiro_seagrass),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: SG_A,SG_B,N,P,T,E_par
      real(rk) :: mu_SG_A,k_I,k_resp,Tau,dD_A_N,dD_B_N,T_fac

      real(rk), parameter :: T_ref       = 20._rk
      real(rk), parameter :: photon_to_P = 5500._rk
      real(rk), parameter :: O2_to_P     = 716._rk
      real(rk), parameter :: photon_to_W = 86400._rk/563900._rk ! from W m-2 to mol m-2 d-1

      _BOTTOM_LOOP_BEGIN_

      ! Retrieve current state and environment
      _GET_BOTTOM_(self%id_SG_A,SG_A)    ! leaf biomass in g N m-2
      _GET_BOTTOM_(self%id_SG_B,SG_B)    ! root biomass in g N m-2
      _GET_BOTTOM_(self%id_N,N)          ! porewater DIN in mmol m-3
      _GET_BOTTOM_(self%id_P,P)          ! porewater DIP in mmol m-3
      _GET_BOTTOM_(self%id_E_par,E_par)  ! PAR at top of canopy in W m-2
      _GET_(self%id_T,T)                 ! bottom temperature in degrees Celsius

      T_fac = self%Q_10**((T-T_ref)/10._rk)

      k_resp = 2*(self%E_comp*self%A_par*self%Omega_SG*self%sin_beta_blade - photon_to_P/N_to_P/M_N*T_fac*(self%zeta_SG_A+self%f_below/(1._rk-self%f_below)*self%zeta_SG_B))*SG_A
      k_I = E_par*photon_to_W*(1.0_rk-exp(-self%A_par*self%Omega_SG*self%sin_beta_blade*SG_A))
      mu_SG_A = min(T_fac*self%mu_max_SG*N/(self%K_SG_N+N), T_fac*self%mu_max_SG*P/(self%K_SG_P+P), N_to_P/photon_to_P*M_N*max(0.0_rk,k_I-k_resp)/SG_A)
      Tau = T_fac*(self%f_below - SG_B/(SG_A+SG_B))*(SG_A+SG_B)*self%tau_tran
      _SET_BOTTOM_DIAGNOSTIC_(self%id_E_par_below,E_par*exp(-self%A_par*self%Omega_SG*self%sin_beta_blade*SG_A))

      _ADD_BOTTOM_SOURCE_(self%id_N_source,         -1._rk/M_N*mu_SG_A*SG_A/86400._rk*1000._rk)  ! in mmol m-2 s-1
      _ADD_BOTTOM_SOURCE_(self%id_P_source,  -1._rk/N_to_P/M_N*mu_SG_A*SG_A/86400._rk*1000._rk)  ! in mmol m-2 s-1
      _ADD_BOTTOM_FLUX_(self%id_DIC,        -C_to_P/N_to_P/M_N*mu_SG_A*SG_A/86400._rk*1000._rk)  ! in mmol m-2 s-1
      _ADD_BOTTOM_FLUX_(self%id_O2,         O2_to_P/N_to_P/M_N*mu_SG_A*SG_A/86400._rk*1000._rk)  ! in mmol O2 m-2 s-1
      _ADD_BOTTOM_SOURCE_(self%id_SG_A,(mu_SG_A*SG_A - T_fac*self%zeta_SG_A*(SG_A-self%f_seed/self%Omega_SG*(1.0_rk-self%f_below)) - Tau)/86400._rk)
      _ADD_BOTTOM_SOURCE_(self%id_SG_B,(             - T_fac*self%zeta_SG_B*(SG_B-self%f_seed/self%Omega_SG*self%f_below         ) + Tau)/86400._rk)

      ! Detritus production (separate leaf and root)
      ! Also convert from d-1 to s-1, and from g N to mmol N
      dD_A_N = T_fac*(self%zeta_SG_A*(SG_A-self%f_seed/self%Omega_SG*(1.0_rk-self%f_below)))/86400._rk*1000._rk/M_N
      dD_B_N = T_fac*(self%zeta_SG_B*(SG_B-self%f_seed/self%Omega_SG*self%f_below         ))/86400._rk*1000._rk/M_N
      _ADD_BOTTOM_FLUX_(self%id_D_A_N,                dD_A_N)
      _ADD_BOTTOM_FLUX_(self%id_D_A_P,   1._rk/N_to_P*dD_A_N)
      _ADD_BOTTOM_FLUX_(self%id_D_A_C,  C_to_P/N_to_P*dD_A_N)
      _ADD_BOTTOM_SOURCE_(self%id_D_B_N,              dD_B_N)
      _ADD_BOTTOM_SOURCE_(self%id_D_B_P, 1._rk/N_to_P*dD_B_N)
      _ADD_BOTTOM_SOURCE_(self%id_D_B_C,C_to_P/N_to_P*dD_B_N)

      _BOTTOM_LOOP_END_

   end subroutine do_bottom

end module
