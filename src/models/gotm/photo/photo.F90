#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gotm_photo --- GOTM photo inhibition model, ported to FABM by Jorn Bruggeman
!
! !INTERFACE:
   module gotm_photo
!
! !DESCRIPTION:
! This module calculates the inhibition parameter,
! $Y$, and the photosynthetic production rates $P$ according to a
! theory of \cite{DenmanMarra1986a} as adapted by
! \cite{Nagaietal2003a}.  The photosynthetically available radiation
! $PAR$ is given as \begin{equation} PAR(z)=I_0 (1-a)
! \exp\left(\frac{z}{\gamma_2} \right) \; , \end{equation} where $z$
! denotes the vertical coordinate, $(1-a)$ the short-wave fraction of
! visible light identified with $PAR$, $\gamma_2$ the penetration
! depth of the short-wave fraction, and $I_0$ the incoming solar
! radiation. The parameters $a$ and $\gamma_2$ have to be provided
! separately here via the namelist of this model, but they should be
! made consistent with the values chosen in {\tt obs.inp} via the
! Jerlov water class.
!
! Note that radiation units used by \cite{Nagaietal2003a}
! have been converted to units W m$^{-2}$ s$^{-1}$ using the
! conversion relation $1$ W m$^{-2}$ s$^{-1}$ $ \approx 0.2174 \, \mu$E
! m$^{-2}$ s$^{-1}$, which is strictly valid only for a wave length of
! $550$ nm.
!
! Instantaneous photosynthetic production rates for fully uninhibited and
! fully inhibited cells are
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! P_d = P_{dm} \left(1-\exp\left(-\frac{PAR}{E_d}  \right)  \right) \; , \\ \\
! \displaystyle
! P_l = P_{lm} \left(1-\exp\left(-\frac{PAR}{E_l}  \right)  \right) \; ,
! \end{array}
! \end{equation}
! respectively, with the maximum production rates, $P_{dm}$ and
! $P_{lm}$, and the saturation values, $E_d$ and $E_l$. In contrast to
! \cite{Nagaietal2003a}, here production is in carbon units (converted
! from the oxygen units by using the appropriate molar weights).
!
! The local production rate is then calculated by
! \begin{equation}
!   P=P_d + Y (P_l-P_d) \; ,
! \end{equation}
! where $0 \le Y \le 1$ denotes the instantaneous inhibition
! factor determining whether cells suffer from photoinhibition.
! The inhibition parameter $Y$ is computed from a linear first-order
! differential equation,
! \begin{equation}
! \frac{\mbox{d} Y}{\mbox{d} t} = \frac{1}{\gamma} (X-Y) \; ,
! \end{equation}
! where $\gamma$ is the response time scale and $X$ is the
! inhibition factor after full adaption to the local conditions.
! It is computed from
! \begin{equation}
!   X=1-\exp \left(-\left(\frac{\max\{PAR,E_b\}-E_b}{E_b}\right)^2 \right).
! \end{equation}
!
! From these production parameters, the cell carbon content is then computed
! according to
! \begin{equation}
!  \deriv{C}{t} = P - \mu C \; ,
! \end{equation}
! where $\mu$ denotes a constant respiration rate.
!
! If the parameter {\tt photo\_mig} is set {\tt .true.} particles
! migrate vertically according to the following algorithm. From
! 06:00 am to 06:00 pm particles move upwards unless their individual
! photoinhibition parameter, $Y$, exceeds the threshold $Y_h$. Those
! particles move downward, i.e.\ away from bright regions. From
! 06:00 pm to 06:00 am, all particles move downward. For all types
! of motion the migration speed, $w_c$, is identical and constant.
!
! !USES:
   use fabm_types

! !REVISION HISTORY:!
!  Original author(s): Lars Umlauf, Hans Burchard, Karsten Bolding
!
   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model), public :: type_gotm_photo
!     Variable identifiers
      type (type_diagnostic_variable_id)   :: id_X
      type (type_state_variable_id)        :: id_Y
      type (type_diagnostic_variable_id)   :: id_Prod
      type (type_state_variable_id)        :: id_Carb
      type (type_horizontal_dependency_id) :: id_swr0
      type (type_dependency_id)            :: id_depth
      type (type_global_dependency_id)     :: id_yearday

!     Model parameters
      logical  :: photo_mig
      real(rk) :: gamma
      real(rk) :: pdm
      real(rk) :: plm
      real(rk) :: mu
      real(rk) :: ed
      real(rk) :: el
      real(rk) :: eb
      real(rk) :: yl
      real(rk) :: yh
      real(rk) :: w_c
      real(rk) :: aa
      real(rk) :: g2

   contains

      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement

   end type
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the photoinhibition model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_gotm_photo), intent(inout), target :: self
   integer,                 intent(in)            :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%get_parameter(self%photo_mig, 'photo_mig', '', 'vertical migration in response to photoinhibition', default=.false.)
   call self%get_parameter(self%gamma, 'gamma', 's', 'photoresponse time', default=3600._rk)
   call self%get_parameter(self%pdm, 'pdm', 'pg C/cell/hour', 'uninhibited production', default=4.0_rk)
   call self%get_parameter(self%plm, 'plm', 'pg C/cell/hour', 'inhibited production', default=0.25_rk)
   call self%get_parameter(self%mu, 'mu', '1/hour', 'respiration rate', default=0.1_rk)
   call self%get_parameter(self%ed, 'ed', 'W/m2', 'radiation threshold', default=163.05_rk)
   call self%get_parameter(self%el, 'el', 'W/m2', 'radiation threshold', default=163.05_rk)
   call self%get_parameter(self%eb, 'eb', 'W/m2', 'inhibition threshold', default=43.48_rk)
   call self%get_parameter(self%yl, 'yl', '-', 'inhibition threshold for downward migration', default=0.2_rk)
   call self%get_parameter(self%yh, 'yh', '-', 'inhibition threshold for upward migration', default=0.2_rk)
   call self%get_parameter(self%w_c, 'w_c', 'm/s', 'migration speed', default=0.0003_rk)
   call self%get_parameter(self%aa, 'aa', '-', 'non-visible fraction of shortwave radiation', default=0.78_rk)
   call self%get_parameter(self%g2, 'g2', 'm', 'penetration depth of photosynthetically active radiation', default=7.9_rk)

   call self%register_diagnostic_variable(self%id_X, 'X', '-', 'fully adapted inhibition factor', specific_to=id_water_parcel)
   call self%register_state_variable(self%id_Y, 'Y', '-', 'inhibition factor', 0.0_rk, specific_to=id_water_parcel)
   call self%register_diagnostic_variable(self%id_Prod, 'Prod', 'pg C/hour/vol', 'production')
   call self%register_state_variable(self%id_Carb, 'Carb', 'pg C/vol', 'carbon concentration', 0.0_rk, specific_to=id_water_parcel)

   call self%register_dependency(self%id_depth, standard_variables%depth)
   call self%register_dependency(self%id_swr0, standard_variables%surface_downwelling_shortwave_flux)
   call self%register_dependency(self%id_yearday, standard_variables%number_of_days_since_start_of_the_year)

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of photoinhibition model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_gotm_photo), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk)            :: swr0, depth, par, x, y, p, pd, pl, c
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

      _GET_HORIZONTAL_(self%id_swr0,swr0)
      _GET_(self%id_depth,depth)

      par = swr0*((1.0_rk - self%aa)*exp(-depth/self%g2))

      ! inhibition state
      x     = 1._rk - exp(-((max(par,self%eb)-self%eb)/self%eb)**2) ! adapted

      _GET_(self%id_Y,y)                ! individual

      ! production
      pd    = self%pdm*(1._rk-exp(-par/self%ed))               ! uninhibited production

      pl    = self%plm*(1._rk-exp(-par/self%el))               ! inhibited production

      _SET_ODE_(self%id_Y,1.0_rk/self%gamma*(x-y))                  ! update ihibition

      p     = pd+y*(pl-pd)                        ! instantaneous production

!     update carbon content
      _GET_(self%id_Carb,c)                ! individual
      _SET_ODE_(self%id_Carb,(p - self%mu*c)/3600._rk)

      _SET_DIAGNOSTIC_(self%id_X,x)
      _SET_DIAGNOSTIC_(self%id_Prod,p)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_gotm_photo),intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: y, yearday, dayhour, w

      if (.not. self%photo_mig) return

      _LOOP_BEGIN_

         _GET_(self%id_Y,y)   ! get inhibition state

         _GET_GLOBAL_(self%id_yearday,yearday)
         dayhour = 24*mod(yearday, 1.0_rk)

         if (dayhour > 6.0_rk .and. dayhour < 18.0_rk) then

            if (y > self%yh)  then
               w = -self%w_c          ! it's day but it's too bright: move down
            else
               w = self%w_c           ! it's day and not too bright: move up
            endif

         else

            w = -self%w_c             ! it's night: move down

         end if

         _SET_VERTICAL_MOVEMENT_(self%id_Y,w)
         _SET_VERTICAL_MOVEMENT_(self%id_Carb,w)

      _LOOP_END_

   end subroutine get_vertical_movement

!-----------------------------------------------------------------------

   end module gotm_photo

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

