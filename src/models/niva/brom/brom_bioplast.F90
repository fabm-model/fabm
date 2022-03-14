#include "fabm_driver.h"

module fabm_niva_brom_bioplast

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_niva_brom_bioplast
      
       ! Add variable identifiers and parameters here.
       type (type_state_variable_id)       :: id_mp_free, id_mp_biof, id_mp_het, id_mp_det
       type (type_state_variable_id)       :: id_Phy, id_Het, id_POM, id_DOM
       ! diagnostic dependences defined here
       type (type_diagnostic_variable_id)  :: id_MP_TOT, id_free2biof, id_free2het, id_biof2het
       type (type_diagnostic_variable_id)  :: id_het2det, id_het2biof, id_det2free, id_det2het        
       type (type_diagnostic_variable_id)  :: id_dMP_free,id_dMP_biof,id_dMP_het,id_dMP_det,id_dMP_tot
       ! diagnostic dependences defined in other modules
       type (type_dependency_id)           :: id_MortHet, id_GrazPhy, id_GrazPOM, id_Autolys
       type (type_dependency_id)           :: id_GrowthPhy, id_RespPhy, id_ExcrPhy, id_MortPhy
       type (type_dependency_id)           :: id_POM_decay_ox, id_POM_decay_denitr, id_RespHet 
      ! global dependencies
       type (type_dependency_id)           :: id_temp,id_par,id_depth
       ! Parameters
       real(rk) :: MP_free_decay, MP_biof_decay, MP_det_decay
       real(rk) :: MP_free_biof        ! Flow from MP_free to MP_biof
       real(rk) :: thr_free, thr_biof  ! threshold for ingestion of MP_free and MP_biof by Het
       real(rk) :: max_free_het, max_biof_het, spec_het_det, hz, uz
       real(rk) :: rate_het_filtr, Kmp_phy
       
       ! Sinking
       real(rk) :: Wfree, Wbiof, Whet, Wpom
       
   contains
      procedure :: initialize
      
      ! Reference model procedures here.
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_niva_brom_bioplast), intent(inout), target :: self
      integer,                         intent(in)            :: configunit
    
      ! Local variables
      real(rk),parameter :: d_per_s = 1.0_rk/86400.0_rk
      
      ! Register model parameters and variables here.

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_MP_TOT,'MP_TOT','g/m^3',  'MP_TOT: Total Microplastic')
      call self%register_diagnostic_variable(self%id_free2biof,'free2biof','g/m^3',  'free2biof: MP_free to MP_biof')
      call self%register_diagnostic_variable(self%id_free2het,'free2het','g/m^3',  'free2het: MP_free to MP_het')
      call self%register_diagnostic_variable(self%id_biof2het,'biof2het','g/m^3',  'biof2het: MP_biof to MP_het')
      !call self%register_diagnostic_variable(self%id_det2het,'det2het','g/m^3',  'det2het: MP_det to MP_het')
      call self%register_diagnostic_variable(self%id_het2det,'het2det','g/m^3',  'het2det: MP_het to MP_det')
      call self%register_diagnostic_variable(self%id_het2biof,'het2biof','g/m^3',  'het2biof: MP_het to MP_biof')
      call self%register_diagnostic_variable(self%id_det2free,'det2free','g/m^3',  'det2free: MP_det to MP_free')
      call self%register_diagnostic_variable(self%id_dMP_free,'dMP_free','g/m^3',  'dMP_free: increament of MP_free')
      call self%register_diagnostic_variable(self%id_dMP_biof,'dMP_biof','g/m^3',  'dMP_biof: increament of MP_biof')
      call self%register_diagnostic_variable(self%id_dMP_het,'dMP_het','g/m^3',  'dMP_het: increament of MP_het')
      call self%register_diagnostic_variable(self%id_dMP_det,'dMP_det','g/m^3',  'dMP_det: increament of MP_det')
      call self%register_diagnostic_variable(self%id_dMP_tot,'dMP_tot','g/m^3',  'dMP_tot: increament of MP_tot')

      ! Dependency
      call self%register_state_dependency(self%id_Phy, 'PHY', 'mmol/m**3','Phy')
      call self%register_state_dependency(self%id_Het, 'HET', 'mmol/m**3','Het')
      call self%register_state_dependency(self%id_POM,'POM','mmol/m**3','particulate organic matter')
      call self%register_state_dependency(self%id_DOM,'DOM','mmol/m**3','dissolved organic matter')
      
      ! Diagnostic dependency
      call self%register_dependency(self%id_MortHet, 'MortHet', '1/d','MortHet')
      call self%register_dependency(self%id_POM_decay_ox, 'POM_decay_ox', '1/d','POM_decay_ox')
      call self%register_dependency(self%id_POM_decay_denitr, 'POM_decay_denitr', '1/d','POM_decay_denitr')
      call self%register_dependency(self%id_GrazPhy, 'GrazPhy', '1/d','GrazPhy')
      call self%register_dependency(self%id_GrazPOM, 'GrazPOM', '1/d','GrazPOM')
      call self%register_dependency(self%id_Autolys, 'Autolys', '1/d','Autolys')
      call self%register_dependency(self%id_RespHet, 'RespHet', '1/d','RespHet')
      call self%register_dependency(self%id_MortPhy, 'MortPhy', '1/d','MortPhy')
      call self%register_dependency(self%id_GrowthPhy, 'GrowthPhy', '1/d','GrowthPhy')
      call self%register_dependency(self%id_RespPhy, 'RespPhy', '1/d','RespPhy')
      call self%register_dependency(self%id_ExcrPhy, 'ExcrPhy', '1/d','ExcrPhy')

      !Register parameters
      call self%get_parameter(self%MP_free_decay, 'MP_free_decay', '1/d', &
                             'MP_free_decay: Rate of free microplastic degradation',      0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%MP_biof_decay, 'MP_biof_decay', '1/d', &
                             'MP_biof_decay: Rate of biofouled microplastic degradation', 0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%MP_det_decay, 'MP_det_decay', '1/d', &
                             'MP_det_decay: Rate of microplastic in detrit degradation',  0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%MP_free_biof, 'MP_free_biof', '1/d', &
                             'MP_free_biof: Rate of biofouling of MP_free',  0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%rate_het_filtr, 'rate_het_filtr', 'm**3/d/mg', &
                             'rate_het_filtr: Rate of filtration by het', 0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%max_free_het, 'max_free_het', '1/d', 'maximum rate of ingestion of free by het', 0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%max_biof_het, 'max_biof_het', '1/d', 'maximum rate of ingestion of biof by het', 0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%thr_free, 'thr_free', 'mg/L', 'thr constant for ingestion of MP_free by het', 0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%thr_biof, 'thr_biof', 'mg/L', 'thr constant for ingestion of MP_biof by het', 0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%spec_het_det, 'spec_het_det', '1/d', 'specific rate of het death', 0.000001_rk, scale_factor=d_per_s)
      call self%get_parameter(self%Kmp_phy, 'Kmp_phy', 'nd', 's"partitioning" coefficient for MP at Phy', default=10000._rk)
      call self%get_parameter(self%uz,        'uz',        'nd',  'Food absorbency for Het',                             default=0.5_rk)
      call self%get_parameter(self%hz,        'hz',        'nd',  'Ratio betw. diss. and part. excretes of Het ',        default=0.5_rk)
   
      ! sinking
      call self%get_parameter(self%Wfree, 'Wfree', 'm/s', 'vertical velocity of MP_free (<0 for sinking)', default=-0.1_rk,scale_factor=d_per_s)
      call self%get_parameter(self%Wbiof, 'Wbiof', 'm/s', 'vertical velocity of MP_biof (<0 for sinking)', default=-0.1_rk,scale_factor=d_per_s)
      call self%get_parameter(self%Whet, 'Whet', 'm/s', 'vertical velocity of het (<0 for sinking)', default=-0.1_rk,scale_factor=d_per_s)
      call self%get_parameter(self%Wpom, 'Wpom', 'm/s', 'vertical velocity of POM (<0 for sinking)', default=-1.0_rk,scale_factor=d_per_s)
      
      ! Register state variables
      call self%register_state_variable(self%id_mp_free, 'MP_free', 'mg/m^3', 'MP_free: Free Microplastic',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Wfree)
      call self%register_state_variable(self%id_mp_biof, 'MP_biof', 'mg/m^3', 'MP_biof: Microplastic with biofouling',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Wbiof)
      call self%register_state_variable(self%id_mp_het,  'MP_het', 'mg/m^3',  'MP_het: Microplastic, injested by HET',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Whet)
      call self%register_state_variable(self%id_mp_det,  'MP_det', 'mg/m^3', 'MP_det: Microplastic in detrit',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Wpom) ! get verical velocity from POM
      
      ! Register state dependencies
      
      call self%register_dependency(self%id_temp, standard_variables%temperature) 
      call self%register_dependency(self%id_depth,standard_variables%pressure)  
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      
   end subroutine initialize

! Add model subroutines here.
   
   subroutine do(self,_ARGUMENTS_DO_)
       class (type_niva_brom_bioplast),intent(in) :: self
       _DECLARE_ARGUMENTS_DO_
       
! Local 
       real(rk) :: MP_free, MP_biof, MP_het, MP_det
       real(rk) :: MP_biof_conc, MP_het_conc, MP_det_conc
       real(rk) :: dMP_biof_conc, dMP_het_conc, dMP_det_conc, dMP_free_conc
       real(rk) :: MP_TOT, MP_free_old, MP_biof_old, MP_free_new, MP_biof_new !, dMP_biof,dMP_free
       real(rk) :: sha_free, sha_bio
       real(rk) :: Phy, Het, POM, DOM
       real(rk) :: temp, depth, Iz
       
       ! Diagnostic variable dependencies
       real(rk) :: MortHet
       real(rk) :: POM_decay_ox, POM_decay_denitr
       real(rk) :: GrazPhy, GrazPOM, Autolys, RespHet, MortPhy 
       real(rk) :: GrowthPhy, RespPhy, ExcrPhy
       
       real(rk) :: mu_free, mu_biof, mu_het, mu_det  ! rate of microplastic change
       real(rk) :: thr_free, thr_biof  ! min threshold for ingestion of MP_free and MP_biof by Het
       real(rk) :: threshold_f  ! threshold for ingestion of MP_free and MP_biof by Het
       real(rk) :: rate_het_filtr
       real(rk) :: free2het, det2het, het2det, het2biof, det2free, biof2het, free2biof
       real(rk) :: dMP_free, dMP_biof, dMP_het, dMP_det, dMP_tot
       real(rk) :: Decay_MP_free, Decay_MP_biof, Decay_MP_het, Decay_MP_det

       ! Parameters
       real(rk) :: max_free_het ! maximum rate of ingestion by het
       real(rk) :: thr_free_het ! thr constant for ingestion by het
       real(rk) :: spec_het_det ! specific rate of het death ! NOT NEEDED
       real(rk) :: Kmp_phy ! "partitioning" coefficient for MP at Phy
       real(rk) :: hz, uz
       real(rk) :: uMn2lip=0.0084 !coeff.to transfer POM (umol/l N)->(g DryWeight/l) 
       real(rk) :: rho_B_OXYDEP_Phy= 1.4E7 ! Density of (living) phytoplankton [mmolN/m3] (default = 1.4E6 mmolN/m3 from PON default)
       
       _LOOP_BEGIN_

! Obtain concentration of microplastic.
       _GET_(self%id_mp_free, MP_free)
       _GET_(self%id_mp_biof, MP_biof)
       _GET_(self%id_mp_het, MP_het)
       _GET_(self%id_mp_det, MP_det)
       
! Get dependencies
       _GET_(self%id_Phy, Phy) 
       _GET_(self%id_Het, Het)
       _GET_(self%id_POM, POM)
       _GET_(self%id_DOM, DOM)
       
! Environment
       _GET_(self%id_temp,temp)              ! temperature
       _GET_(self%id_depth,depth)            ! depth
       _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
        
! Get diagnostic dependencies
       _GET_(self%id_MortHet, MortHet)
       _GET_(self%id_GrazPhy, GrazPhy)
       _GET_(self%id_Autolys, Autolys)
       _GET_(self%id_RespHet, RespHet)
       _GET_(self%id_GrazPOM, GrazPOM)
       _GET_(self%id_MortPhy, MortPhy)
       _GET_(self%id_GrowthPhy, GrowthPhy)
       _GET_(self%id_RespPhy, RespPhy)
       _GET_(self%id_ExcrPhy, ExcrPhy)
       _GET_(self%id_POM_decay_ox, POM_decay_ox)
       _GET_(self%id_POM_decay_denitr, POM_decay_denitr)

! MP_free biofouling with Phy=============================
       if (Phy.lt.0.00001.or.MP_free.le.0.00001) then 
           free2biof = 0.0_rk
       else
!           threshold_f = MP_free/(Het+0.000001)+self%thr_free
!           free2biof = (GrowthPhy-RespPhy-ExcrPhy-MortPhy-GrazPhy)*0.001*(MP_free/MP_free+1.0_rk)
           free2biof = (GrowthPhy-RespPhy-ExcrPhy-MortPhy-GrazPhy)*self%MP_free_biof*0.5*(1-tanh(0.002-MP_free))
!           free2biof = (GrowthPhy-RespPhy-ExcrPhy-MortPhy-GrazPhy)/Phy*max(0.00000000001, MP_biof)
!           Phy*MP_free*Phy*0.5*(1-tanh(self%thr_free-(MP_free/2)**0.75-Phy))
!var3!           free2het =   0.5*(1-tanh(self%thr_free-(MP_free/2)**0.75-Het))*Het
!var2!           free2het =    Het*1/(1 + exp(threshold_f/(MP_free+0.0001)))
!var1!           free2het = self%rate_het_filtr*MP_free/(MP_free+threshold_f)*Het
!           free2het = self%rate_het_filtr*MP_free*MP_free/(MP_free*MP_free+threshold_f)*Het
       endif
       ! here we used the "partitioning approach" adduming that "concentration" of Phy "in" microplastic is Kmp_phy higher that in the water
       !MP_free_old = MP_free
       !MP_biof_old = MP_biof
       !if(Phy<=0.00000001) then 
       ! sha_bio=0. 
       !else
       ! sha_bio=  1000.0_rk*Phy/rho_B_OXYDEP_Phy! Volume occupied by Phy
       !endif 
       !sha_free = 1000.0_rk-sha_bio
       !MP_free_new = max(0.,(MP_free_old+MP_biof_old)*sha_free /(sha_free +sha_bio/self%Kmp_phy)) !
       !MP_biof_new = (MP_free_old+MP_biof_old)-MP_free_new
       !dMP_free = -MP_free_old + MP_free_new
       !dMP_biof = -MP_biof_old + MP_biof_new
! MP_free ingessioon by Het===============================
       if (MP_free.le.0.00001) then 
           free2het = 0.0_rk
       else
!var4!          free2het =  self%max_free_het*Het*MP_free/(MP_free+self%thr_free) 
            free2het =  self%max_free_het*self%Uz*GrazPhy*MP_free/(MP_free+Het*self%thr_free) 
!           free2het =  self%Uz*(GrazPhy+GrazPOM)*self%max_free_het*MP_free/(MP_free+Het*self%thr_free) 
!           free2het =   self%max_free_het/(1.0_rk +exp(-self%rate_het_filtr*(Het-(self%thr_free-MP_free/(MP_free+1._rk)))))
!var3!           free2het =   0.5*(1-tanh(self%thr_free-(MP_free/2)**0.75-Het))*Het
!var2!           free2het =    Het*1/(1 + exp(threshold_f/(MP_free+0.0001)))
!var1!           free2het = self%rate_het_filtr*MP_free/(MP_free+threshold_f)*Het
       endif

! MP_biof ingestion by Het===============================
       if (Het.lt.0.00001.or.MP_biof.le.0.0000001) then 
           biof2het = 0.0_rk
       else
!           biof2het = 0.0_rk
!           threshold_f = MP_biof/(Het+0.000000001)+self%thr_free
!          biof2het =   self%max_biof_het*MP_biof*Het*0.5*(1-tanh(self%thr_biof-(MP_biof/2)**0.75-Het))
!           biof2het =   self%max_biof_het*Het*0.5*(1-tanh(MP_biof-self%thr_biof))
            biof2het =  self%max_free_het*self%Uz*GrazPhy*MP_biof/(MP_biof+Het*self%thr_biof) 
       endif
! MP_det ingessioon by Het===============================
       if (MP_det.le.0.000001) then 
           det2het = 0.0_rk
       else
           det2het =  self%max_free_het*self%Uz*GrazPOM*MP_det/(MP_det+Het*self%thr_biof)           
 !var4!           det2het =  self%max_free_het*self%Uz*(GrazPhy+GrazPOM)*MP_det !/(MP_free+Het*self%thr_free) 
       endif

! MP_det release from MP_het=============================
       if (Het.gt.0.000001) then 
  !var4!          het2det =  MP_het*(self%max_free_het*(GrazPhy+GrazPOM)*(1.-self%uz)*(1.-self%hz)+MortHet)/Het !+increase due to respiraion ("concentrtation")
  !var4!          het2biof = MP_het*(self%max_free_het*GrazPhy+(GrazPOM)*(1.-self%uz)*self%hz+RespHet)/Het 
           het2det =  self%max_free_het*(free2het+det2het)*(1.-self%uz)*(1.-self%hz)+MP_het*MortHet/Het !+increase due to respiraion ("concentrtation")
           het2biof = self%max_free_het*(free2het+det2het)*(1.-self%uz)*self%hz     +MP_het*RespHet/Het 
       else
           het2det  = 0._rk
           het2biof = 0._rk
       endif

! MP_free release from MP_detritus========================
       if (POM.gt.0.000001) then
           det2free = self%max_free_het*(det2het+MP_det*(POM_decay_ox+POM_decay_denitr+Autolys)/POM)
       else
           det2free = 0._rk
       endif
   
! MP vatiables degradation================================
       Decay_MP_free = -self%MP_free_decay*Iz/25.*exp(1._rk-Iz/25._rk)*MP_free !depends on light
       Decay_MP_biof = -(self%MP_free_decay*Iz/25.*exp(1._rk-Iz/25._rk)+self%MP_biof_decay)*MP_biof
       Decay_MP_het  = -self%MP_det_decay *MP_het 
       Decay_MP_det =  -self%MP_det_decay *MP_det

! ========================================================
       MP_TOT = MP_free + MP_het + MP_det + MP_biof !from previous timestep
       
        dMP_free = Decay_MP_free - free2biof + det2free - free2het + het2biof
        dMP_biof = Decay_MP_biof + free2biof - biof2het ! + het2biof
        dMP_het  = Decay_MP_het + free2het + biof2het - het2det - het2biof
        dMP_det  = Decay_MP_det + het2det  - det2free
        dMP_tot  = dMP_free + dMP_biof + dMP_het + dMP_det
        
! Send rates of change to FABM.
       _ADD_SOURCE_(self%id_mp_free, Decay_MP_free - free2biof + det2free - free2het + het2biof  )
       _ADD_SOURCE_(self%id_mp_biof, Decay_MP_biof + free2biof - biof2het) ! + het2biof)
       _ADD_SOURCE_(self%id_mp_het,  Decay_MP_het + free2het + biof2het - het2det - het2biof + det2het )
       _ADD_SOURCE_(self%id_mp_det,  Decay_MP_det + het2det  - det2free - det2het)
       
! Export diagnostic variables
       _SET_DIAGNOSTIC_(self%id_MP_TOT,MP_TOT)
       _SET_DIAGNOSTIC_(self%id_free2biof,free2biof)
       _SET_DIAGNOSTIC_(self%id_biof2het,biof2het)
       _SET_DIAGNOSTIC_(self%id_free2het,free2het)
       _SET_DIAGNOSTIC_(self%id_het2det,het2det)
       _SET_DIAGNOSTIC_(self%id_het2biof,het2biof)
       _SET_DIAGNOSTIC_(self%id_det2free,det2free)   
       _SET_DIAGNOSTIC_(self%id_dMP_free,dMP_free)  
       _SET_DIAGNOSTIC_(self%id_dMP_biof,dMP_biof)  
       _SET_DIAGNOSTIC_(self%id_dMP_het,dMP_het)  
       _SET_DIAGNOSTIC_(self%id_dMP_det,dMP_det)  
       _SET_DIAGNOSTIC_(self%id_dMP_tot,dMP_tot)  
       
       _LOOP_END_
    end subroutine do

end module
