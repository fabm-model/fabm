#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_klimacampus_phy_feedback --- is the implementation of a
! PND-model including 3 different biological-physical feedback mechanisms
! as described in Sonntag and Hense (2011). Here, only the N2-fixing
! cyanobacteria are included.
!
! !INTERFACE:
   module fabm_klimacampus_phy_feedback
!
! !DESCRIPTION:
!  This PND (phytoplankton (cyanobacteria)-nutrient-detritus) model
!  described here consists of 3 state variables.
!  Cyanobacteria growth is limited by light and temperature.
!  Cyanobacteria (PHY) die and the biomass goes into
!  the detritus (DET) which can be remineralized to nutrients (NUT)
!  We assume that the pool of N2 is infinite and therefore do not
!  introduce an additional state variable for N2. Thus, the mass
!  balance is not closed!
!  For all details of the PND-feedback model implemented here,
!  see Sonntag and Hense (2011).
!
! !USES:
   use fabm_types
!
   implicit none
!
!  default: all is private.
   private
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day  = 86400._rk, secs_pr_hour = 3600._rk
!
! !REVISION HISTORY:!
!  Original author(s): Inga Hense
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_klimacampus_phy_feedback
!     Variable identifiers
      type (type_state_variable_id)      :: id_nut,id_phy,id_det
      type (type_dependency_id)          :: id_par,id_temp
      type (type_diagnostic_variable_id) :: id_NFIX, id_dPAR

!     Model parameters
      real(rk) :: muemax_phy,alpha,mortphy,rem, &
                  topt,tl1,tl2,depo,nbot,albedo_bio,drag_bio
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: get_albedo
      procedure :: get_drag
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the phy-feedback (PND) model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the phy_feedback model namelist is read and the variables
!  exported by the model are registered in FABM.
!
! !INPUT PARAMETERS:
   class (type_klimacampus_phy_feedback), intent(inout),target :: self
   integer,                               intent(in)           :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Inga Hense
!
! !LOCAL VARIABLES:
   real(rk)                  :: nut_initial=4.5   ! nutrients
   real(rk)                  :: phy_initial=0.    ! phytolankton (cyanobacteria)
   real(rk)                  :: det_initial=4.5   ! detritus
   real(rk)                  :: rkc=0.03          ! attenuation coefficient (organic matter)
   real(rk)                  :: alpha=0.3         ! initial slope of the PI-curve
   real(rk)                  :: muemax_phy=0.25   ! maximum growth rate
   real(rk)                  :: mortphy=0.04      ! mortality rate
   real(rk)                  :: rem=0.048         ! remineralization rate
   real(rk)                  :: topt=25.          ! optimum temperature for growth
   real(rk)                  :: tl1=2.            ! slope (temp. function)
   real(rk)                  :: tl2=3.            ! slope (temp. function)
   real(rk)                  :: w_phy=1.          ! positive buoyancy
   real(rk)                  :: w_det=-18.        ! sinking of detritus
   real(rk)                  :: depo=0.05         ! detritus burial rate
   real(rk)                  :: nbot=35.          ! bottom nutrient concentration
   real(rk)                  :: albedo_bio=0.002  ! factor for albedo changes through surface phy
   real(rk)                  :: drag_bio=0.05     ! factor for drag coef. changes through surface phy

   namelist /klimacampus_phy_feedback/ &
              nut_initial,phy_initial,det_initial,rkc,alpha,   &
              muemax_phy,mortphy,rem,topt,tl1,tl2,w_phy,w_det, &
              depo,nbot
!EOP
!-----------------------------------------------------------------------
!BOC
!  Read the namelist
   if (configunit>0) read(configunit,nml=klimacampus_phy_feedback,err=99,end=100)

!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   self%alpha       = alpha/secs_pr_day
   self%muemax_phy  = muemax_phy/secs_pr_day
   self%mortphy     = 1.0_rk/mortphy/secs_pr_day
   self%rem         = rem/secs_pr_day
   self%topt        = topt
   self%tl1         = tl1
   self%tl2         = tl2
   self%depo        = depo/secs_pr_day
   self%nbot        = nbot
   self%albedo_bio  = albedo_bio
   self%drag_bio    = drag_bio

!  Register state variables
   call self%register_state_variable(self%id_nut,'nut','mmol/m**3','nutrients',     &
                                    nut_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_phy,'phy','mmol/m**3','phytoplankton', &
                                    phy_initial,minimum=0.0_rk,vertical_movement=    &
                                    w_phy/secs_pr_day,specific_light_extinction=rkc)
   call self%register_state_variable(self%id_det,'det','mmol/m**3','detritus',      &
                                    det_initial,minimum=0.0_rk,vertical_movement=    &
                                    w_det/secs_pr_day,specific_light_extinction=rkc)

!  Register diagnostic variables
   call self%register_diagnostic_variable(self%id_NFIX,'NFIX','mmol/m**3',        &
                     'nitrogen fixation',                                            &
                     time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2',            &
                     'photosynthetically active radiation',                          &
                     time_treatment=time_treatment_averaged)

!  Register contribution to conserved quantities
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_nut)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_phy)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_det)

!  Register environmental dependencies
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   return

99 call self%fatal_error('klimacampus_phy_feedback_init','Error reading namelist klimacampus_phy_feedback.')

100 call self%fatal_error('klimacampus_phy_feedback_init','Namelist klimacampus_phy_feedback was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of  model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
! The sources and sinks of the PND model are calculated. Please note,
! that phy are cyanobacteria which are assumed to fix nitrogen and bring
! nitrogen into the system.
!
! !INPUT PARAMETERS:
   class (type_klimacampus_phy_feedback),       intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Inga Hense
!
! !LOCAL VARIABLES:
   real(rk)                   :: nut,phy,det,par,temp
   real(rk)                   :: ta_phy,mue_par_phy,mue_phy
!EOP
!-----------------------------------------------------------------------
!BOC
!  Enter spatial loops (if any)
   _LOOP_BEGIN_

!  Retrieve current (local) state variable values.
   _GET_(self%id_nut,nut) ! nutrient
   _GET_(self%id_phy,phy) ! phytoplankton
   _GET_(self%id_det,det) ! detritus

!  Retrieve current environmental conditions.
   _GET_ (self%id_par,par)    ! local photosynthetically active radiation
   _GET_ (self%id_temp,temp)  ! local temperature

!  growth rates and limitation functions for temperature and light.
   ta_phy      = exp(-(temp-self%topt)**4/(self%tl1-sign(self%tl2,temp-self%topt))**4) ! temperature dependence

   mue_par_phy = (self%alpha * par) /          &                                       ! light limitation
                    (self%muemax_phy**2 + (self%alpha*par)**2)**0.5
   mue_phy     = self%muemax_phy * ta_phy * mue_par_phy                                ! actual growth rate

!  Set temporal derivatives
  _SET_ODE_(self%id_phy, + mue_phy*phy - self%mortphy*phy                  )
  _SET_ODE_(self%id_det,               + self%mortphy*phy - self%rem*det   )
  _SET_ODE_(self%id_nut,                                  + self%rem*det   )

!  Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_NFIX,mue_phy*phy)

!  Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the albedo coefficient
!
! !DESCRIPTION:
! This is the second feedback considered here: feedback through changes in the albedo,
! leading to changes in light reflection and thus temperature.
!
! !INTERFACE:
   subroutine get_albedo(self,_ARGUMENTS_GET_ALBEDO_)
!
! !INPUT PARAMETERS:
   class (type_klimacampus_phy_feedback), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_ALBEDO_
!
! !REVISION HISTORY:
!  Original author(s): Inga Hense
!
! !LOCAL VARIABLES:
   real(rk)                     :: phys
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

!  Retrieve current (local) state variable values.
   _GET_(self%id_phy,phys)                 ! surface phytoplankton

!  Changes in Albedo due to phytoplankton surface scums.
!  Phy must be the surface concentration!
   _SET_ALBEDO_(self%albedo_bio*phys)

!  Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine get_albedo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the drag coefficient
!
! !DESCRIPTION:
!  This is the third feedback considered here: feedback through changes in the drag
!  coefficient, leading to changes in the momentum flux.
!
! !INTERFACE:
   subroutine get_drag(self,_ARGUMENTS_GET_DRAG_)
!
! !INPUT PARAMETERS:
   class (type_klimacampus_phy_feedback), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_DRAG_
!
! !REVISION HISTORY:
!  Original author(s): Inga Hense
!
! !LOCAL VARIABLES:
   real(rk)                     :: phys
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

!  Retrieve current (local) state variable values.
   _GET_(self%id_phy,phys) ! surface phytoplankton

!  Changes in drag coefficient due to surface scums.
!  Phy must be the surface concentration!
   _SCALE_DRAG_(max(1.0_rk-self%drag_bio*phys,0.0_rk))

!  Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine get_drag
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
! Detritus burial and nutrient restoring are considered to achieve a quasi-steady state.
!
! !INPUT PARAMETERS:
   class (type_klimacampus_phy_feedback), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s): Inga Hense
!
! !LOCAL VARIABLES:
   real(rk)                   :: deb,nub
!EOP
!-----------------------------------------------------------------------
!BOC
!  Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

!  Retrieve current (local) state variable values.
   _GET_(self%id_nut,nub)
   _GET_(self%id_det,deb)

   _SET_BOTTOM_EXCHANGE_(self%id_det,-self%depo*deb*deb)           !detritus burial
   _SET_BOTTOM_EXCHANGE_(self%id_nut,(self%nbot-nub)/secs_pr_hour) !nutrient restoring

!  Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------

   end module fabm_klimacampus_phy_feedback

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
