!-----------------------------------------------------------------------
! OXYDEP_COD is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Elizaveta Protsenko
!-----------------------------------------------------------------------
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:
!
! !INTERFACE:
   module fabm_niva_oxydep_cod
!
! !DESCRIPTION:
!
! OXYDEP_COD parameterize  the Chemical Oxygen Demand, COD (https://en.wikipedia.org/wiki/Chemical_oxygen_demand) 
! OM mineralization, nitrification, and oxidation of reduced specied of S, Mn, Fe, present in suboxic conditions.
! OXYDEP_COD consists of 1 state variables ( in oxygen-units):
! - CHON - is a theoretical organic compound CnHaObNc that can be can be fully oxidized to inorganic nutrients 
!   with a strong oxidizing agent under acidic condition (oxygen).
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev, Elizaveta Protsenko
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_oxydep_cod
!     Variable identifiers
      type (type_state_variable_id)        :: id_oxy, id_nut, id_dom, id_CHON
      type (type_dependency_id)            :: id_par,id_temp, id_salt
      type (type_horizontal_dependency_id) :: id_windspeed
      type (type_diagnostic_variable_id)   :: id_CHON_decay_ox,id_CHON_decay_denitr
!     Model parameters
      !---Organic matter mineralization---- !
       real(rk) :: r_CHON_nut_oxy,r_CHON_nut_nut,r_CHON_dom_oxy, r_CHON_dom_nut, Tda, beta_da
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
! !IROUTINE: Initialise the OXYDEP-COD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_oxydep_cod), intent(inout), target :: self
   integer,                      intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
       real(rk),parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
  ! CHON
   call self%get_parameter(self%r_CHON_nut_oxy, 'r_CHON_nut_oxy',     '1/d', 'Specific rate of CHON mineralization',          default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%beta_da,        'beta_da',       'nd',  'Coefficient for dependence of mineralization on t ', default=20._rk)
   call self%get_parameter(self%Tda,            'Tda',           'nd',  'Coefficient for dependence of mineralization on t ', default=13._rk)
   ! Register state variables
   call self%register_state_variable(self%id_CHON,'CHON','mmol/m**3','CHON  oxidizable compound ', 0.0_rk, minimum=0.0_rk)
   ! Register link to external variables
   call self%register_state_dependency(self%id_oxy,'Oxy','mmol/m**3','OXY')
!   call self%register_state_dependency(self%id_oxy,'NUT','mmol/m**3','NUT')
!   call self%register_state_dependency(self%id_oxy,'DOM','mmol/m**3','DOM')
   ! Register diagnostic variables
call self%register_diagnostic_variable(self%id_CHON_decay_ox,'CHON_decay_ox','mmol/m**3/d',  'CHON_decay_ox,  Mineralization of CHON with oxygen',           &
                    output=output_time_step_integrated)
!call self%register_diagnostic_variable(self%id_CHON_decay_denitr,'CHON_decay_denitr','mmol/m**3/d',  'CHON_decay_denitr,  Mineralization of CHON with nitrate',           &
!                    output=output_time_step_integrated)
   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
!   call self%register_dependency(self%id_windspeed,standard_variables%wind_speed)
!   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

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
   class (type_niva_oxydep_cod),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk)                   :: CHON,oxy,nut,dom,t,iopt

   real(rk) :: doxy,dCHON
 ! Rates of biogeochemical processes
   real(rk) :: CHON_decay_ox            ! oxic mineralization of CHON and ammonification (1/d)
!   real(rk) :: CHON_decay_denitr        ! suboxic mineralization of CHON (denitrification) (1/d)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_CHON,CHON)
!   _GET_(self%id_nut,nut)
!   _GET_(self%id_dom,dom)

   ! Retrieve current environmental conditions.
!   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
   _GET_(self%id_temp,t)              ! temperature

!--------------------------------------------------------------
! Oxic mineralization of CHON depends on T
   CHON_decay_ox   = self%r_CHON_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*CHON
! Suboxic mineralization depends on T,O2,NO3/NO2
   !CHON_decay_denitr = self%r_pom_nut_nut*(1.+self%beta_da*yy(self%tda,t)) &
   !                        * (0.5-0.5*tanh(self%O2LimC-oxy)) &
   !                        * (1-tanh(1.-nut))*pom
! Mineralization of OM, ammonification and growth of NUT


! Now we can summarize processes and write state variables sink/sources:
!--------------------------------------------------------------
! OXY
!--------------------------------------------------------------
! Changes of OXY due to OM production and decay!
   doxy  = -CHON_decay_ox
   dCHON = -CHON_decay_ox
!   dnut = 
!   ddom = 
!--------------------------------------------------------------



!derivatives for FABM
   _SET_ODE_(self%id_oxy, doxy)
   _SET_ODE_(self%id_CHON,dCHON)
!   _SET_ODE_(self%id_nut,dnut)
!   _SET_ODE_(self%id_dom,ddom)

   
   ! Export diagnostic variables
_SET_DIAGNOSTIC_(self%id_CHON_decay_ox,CHON_decay_ox)
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Saturation function squared
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
! !IN2PUT PARAMETERS:
   real(rk), intent(in)                :: a,x
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

   end module fabm_niva_oxydep_cod

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
