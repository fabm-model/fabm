#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
!                  MSI-ERGOM (Baltic Sea Ecosystem Model)                       
!                                                                            
!  Coded by Gennadi Lessin (PML) based upon code provided by Thomas Neumann
!  (IOW) and updates the older version of ERGOM included in GOTM-BIO.
!  This older version is still provided in FABM under the name gotm/ergom for           
!  historical/educational purposes.                                           
!                                                                            
!  A detailed description of the original model is given                     
!  in T. Neumann, W. Fennel and C. Kremp, 2002. Experimental simulations     
!  with an ecosystem model of the Baltic Sea: a nutrient load reduction      
!  experiment. Global Biogeochemical Cycles, 16 (3).                         
!  http://dx.doi.org/10.1029/2001GB001450.                                   
!                                                                            
!  The present version adds oxygen-dependent phosphorus dynamics     
!  between sediment and water and the effect of bio-resuspension, as         
!  described in T. Neumann and G. Schernewski, 2008. Eutrophication in the   
!  Baltic Sea and shifts in nitrogen fixation analyzed with a 3D ecosystem   
!  model, J. Mar. Sys., 74 (1.2), pp. 592-602. 
!
!  Revision history:
!  September 2015, by Dennis Trolle (AU):
!  Implemented a switch for choosing between fresh and marine (default) environments.
!  If "fresh" is selected, oxygen debt (negative O2) due to H2S production is disabled.
!  Added a sediment burial process, and a range of additional diagnostic variables to output, 
!  incl. chlorophyll a, oxygen and nutrients in mass concentration units. 
!  Updated yaml input file with new entries (e.g., sediment burial rate, and phytoplankton
!  carbon:chla ratios for individual phyto groups)
!  May 2016, by Dennis Trolle (AU):
!  Added option to switch on or off n-fixation by cyanobacteria
!  Added settling of diatoms to bottom sediments, where diatoms are converted to fluff once settled
!
! !MODULE: selma
!
! !INTERFACE:
   MODULE selma_zooplankton
!
! !DESCRIPTION:
!
! !USE:
   use fabm_types

   implicit none

   private
!
! !PUBLIC_DERIVED_TYPES:
  type,extends(type_base_model),public :: type_selma_zooplankton
      ! Variable identifiers
      type (type_state_variable_id)             :: id_c
      type (type_state_variable_id),allocatable :: id_prey(:)
      type (type_state_variable_id)             :: id_aa, id_po, id_dd, id_dic, id_o2
      type (type_dependency_id)                 :: id_temp

      ! Model parameters
      integer :: nprey
      real(rk), allocatable :: pref(:)
      real(rk) :: nue,sigma_b
      real(rk) :: iv,graz,toptz,zcl1
      real(rk) :: rfr,rfc
   contains
      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the selma/zooplankton model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!   Here, parameter values are read and the variables exported by the model are registered with FABM
!
! !INPUT PARAMETERS:
   class(type_selma_zooplankton),intent(inout),target :: self
   integer,               intent(in)           :: configunit
   real(rk),parameter :: secs_per_day = 86400._rk
   real(rk) :: c0
   integer :: iprey
   character(len=8) :: strprey

   call self%get_parameter(c0,       'c0',  'mmol N/m3',   'background concentration',    default=0._rk)
   call self%get_parameter(self%rfr, 'rfr', 'mol P/mol N', 'phosphorus : nitrogen ratio', default=0.0625_rk)
   call self%get_parameter(self%rfc, 'rfc', 'mol C/mol N', 'carbon : nitrogen ratio',     default=6.625_rk)
   call self%register_state_variable(self%id_c,'c','mmol N/m3', 'concentration', minimum=0.0_rk, background_value=c0)

   call self%get_parameter(self%nprey, 'nprey', '', 'number of prey', default=1) 
   allocate(self%id_prey(self%nprey))
   allocate(self%pref(self%nprey))
   do iprey=1,self%nprey
      write (strprey,'(i0)') iprey
      call self%register_state_dependency(self%id_prey(iprey),'prey'//trim(strprey), 'mmol N/m3', 'prey '//trim(strprey))
      call self%get_parameter(self%pref(iprey), 'pref'//trim(strprey), '-', 'preference for prey '//trim(strprey), default=1.0_rk) 
   end do

   call self%get_parameter(self%nue,     'nue',     'm3/d/mmol N',    'respiration rate',                default=0.01_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%sigma_b, 'sigma_b', 'm3/d/mmol N',    'mortality rate',                  default=0.03_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%iv,      'iv',      '1/(mmol N/m3)3', 'Ivlev constant, quadratic',       default=1.2_rk)
   call self%get_parameter(self%graz,    'graz',    '1/d',            'grazing rate',                    default=0.5_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%toptz,   'toptz',   'deg C',          'optimal temperature for grazing', default=20._rk)
   call self%get_parameter(self%zcl1,    'zcl1',    '-',              'closure parameter',               default= 50._rk)

   ! Register state variables
   call self%register_state_dependency(self%id_aa, 'aa', 'mmol N/m3', 'ammonium')
   call self%register_state_dependency(self%id_po, 'po', 'mmol P/m3', 'phosphate')
   call self%register_state_dependency(self%id_dd, 'dd', 'mmol N/m3', 'detritus')
   call self%register_state_dependency(self%id_o2, 'o2', 'mmol O2/m3','oxygen')

   ! Contribute to total nitrogen, phosphorus, carbon
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, self%rfr)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c, self%rfc)

   ! Register link to external DIC pool
   call self%register_state_dependency(self%id_dic,standard_variables%mole_concentration_of_dissolved_inorganic_carbon, required=.false.)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)

   END subroutine initialize
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of selma model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class(type_selma_zooplankton), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
    real(rk) :: c, cg
    real(rk) :: prey
    real(rk) :: temp, o2
    real(rk) :: food, zz0, food_eps, gg, lzd, lzn, graz_z, tlim
    integer  :: iprey
    real(rk),parameter :: epsilon = 0.00000000001_rk
!EOP
!-----------------------------------------------------------------------
!BOC
    ! Enter spatial_loops (if any)
    _LOOP_BEGIN_

      ! Retrieve current (local) state variable values
      _GET_(self%id_c,c)                  ! own concentration
      _GET_WITH_BACKGROUND_(self%id_c,cg) ! own concentration including background

      food = 0
      do iprey=1,self%nprey
         _GET_(self%id_prey(iprey),prey)
         food = food + self%pref(iprey)*prey
      end do

      ! Local environment
      _GET_(self%id_o2,o2)     ! oxygen (mmol O2/m3)
      _GET_(self%id_temp,temp) ! temperature (degrees Celsius)

      zz0 = self%zcl1*c*c

      food_eps = max(food, epsilon) ! Be sure food is positive
      gg = self%graz * (1.0_rk - exp(self%iv * food * food * (-1.0_rk))) !Grazing rate

      lzd = self%sigma_b ! Zooplankton mortality rate
      lzn = self%nue     ! Zooplankton respiration rate

      if (o2 <= 0.0_rk) then
         ! Anoxic conditions
         gg = 0.0_rk                ! No grazing
         lzd = 10._rk *self%sigma_b ! Higher mortality
         lzn = 0.0_rk               ! No respiration
      end if

      ! Zooplankton grazing depends on food availability and temperature
      tlim = (1.0_rk + 2.7183_rk/self%toptz/self%toptz * max(temp,0.0_rk)**2 * exp(1.0_rk - temp * 2.0_rk / self%toptz))
      graz_z = gg * cg/food_eps * tlim

      _ADD_SOURCE_(self%id_o2, - 6.625_rk * lzn * zz0)
      _ADD_SOURCE_(self%id_aa, + lzn * zz0)
      _ADD_SOURCE_(self%id_po, self%rfr * lzn * zz0)
      do iprey=1,self%nprey
         _GET_(self%id_prey(iprey),prey)
         _ADD_SOURCE_(self%id_prey(iprey), - graz_z * self%pref(iprey) * prey)
      end do
      _ADD_SOURCE_(self%id_c, graz_z * food - (lzn + lzd) * zz0)
      _ADD_SOURCE_(self%id_dd, + lzd * zz0)

      if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic, self%rfc * lzn * zz0)

   ! Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do
!EOC

!-----------------------------------------------------------------------

  END MODULE selma_zooplankton

!-----------------------------------------------------------------------
