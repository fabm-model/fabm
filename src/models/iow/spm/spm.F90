#include "fabm_driver.h"
#include "fabm.h"

#define STDERR write(*,*)
#define LEVEL0 STDERR
#define LEVEL1 STDERR '   ',
#define LEVEL2 STDERR '       ',
#define LEVEL3 STDERR '           ',
#define LEVEL4 STDERR '               ',
#define FATAL  STDERR 'FATAL ERROR: ',

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_iow_spm --- 1-class SPM model,
!
! !INTERFACE:
   module fabm_iow_spm
!
! !USES:
   use fabm_types

!  default: all is private.
   private
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Manuel Ruiz Villarreal & Richard Hofmeister & Ulf Gräwe
!
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_iow_spm
!     Variable identifiers
      type (type_state_variable_id)                   :: id_spm      !concentrations
      type (type_bottom_state_variable_id)            :: id_pmpool   !sediment pool
      type (type_horizontal_dependency_id)            :: id_taub     !bottom stress
      type (type_dependency_id)                       :: id_temp     !temperature
      type (type_dependency_id)                       :: id_rhow     !density
      type (type_horizontal_diagnostic_variable_id)   :: id_massflux !exchange layer
      type (type_horizontal_diagnostic_variable_id)   :: id_bedload  !bedload

!     Model parameters
      real(rk) :: diameter
      real(rk) :: c_init
      real(rk) :: mass_sed_init
      real(rk) :: tauc_factor
      real(rk) :: erosion_const
      logical  :: cohesive
      integer  :: sinking_method
      integer  :: bottom_stress_method
      logical  :: pm_pool
      real(rk) :: shading
      real(rk) :: tauc_const
      real(rk) :: ws_const
      real(rk) :: rho
      logical  :: consolidate_bed
      integer  :: consolidate_bed_method
      logical  :: bedload
      integer  :: bedload_method
      logical  :: add_to_density

      contains

      procedure :: initialize
      procedure :: do_bottom
      procedure :: get_vertical_movement
      procedure :: get_light_extinction

   end type type_iow_spm
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise the sediment model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the spm namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:

   integer,                        intent(in)            :: configunit
   class (type_iow_spm),           intent(inout), target :: self

!
! !REVISION HISTORY:
!  Original author(s): Ulf Gräwe & Richard Hofmeister
!
! !LOCAL VARIABLES:

   real(rk) :: diameter=100.0_rk         ! mikrons
   real(rk) :: c_init=5.0_rk             ! mg/l
   real(rk) :: mass_sed_init=1000.0_rk   ! kg/m**2
   real(rk) :: tauc_factor=1
   real(rk) :: erosion_const=0.01_rk
   real(rk) :: tauc_const=0.01_rk        ! N/m**2
   real(rk) :: ws_const=0.001_rk         ! m/s
   logical  :: cohesive=.false.
   integer  :: sinking_method=0
   logical  :: bedload=.true.
   integer  :: bedload_method=1
   integer  :: bottom_stress_method=1
   logical  :: pm_pool=.true.
   real(rk) :: shading=1.0_rk            ! 1/m per mg/l
   real(rk) :: rho=2650.0_rk             ! dry bed density kg/m**3
   logical  :: consolidate_bed=.false.   ! mimic consolidation of the bed
   integer  :: consolidate_bed_method=2  ! mimic consolidation of the bed
   logical  :: add_to_density=.false.    ! include density effects in EOS

   namelist /iow_spm/  diameter, &
         c_init, mass_sed_init, tauc_factor, erosion_const, &
         tauc_const, cohesive, sinking_method, bottom_stress_method, &
         pm_pool, shading, ws_const, rho, bedload, bedload_method, &
         consolidate_bed, consolidate_bed_method, add_to_density
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) read(configunit,nml=iow_spm,err=99,end=100)

   ! Store parameter values in our own derived type
   self%diameter               = diameter/1e6_rk ! convert to m
   self%tauc_factor            = tauc_factor
   self%tauc_const             = tauc_const
   self%cohesive               = cohesive
   self%shading                = shading
   self%ws_const               = ws_const
   self%bottom_stress_method   = bottom_stress_method
   self%sinking_method         = sinking_method
   self%pm_pool                = pm_pool
   self%mass_sed_init          = mass_sed_init
   self%rho                    = rho
   self%erosion_const          = erosion_const
   self%consolidate_bed        = consolidate_bed
   self%consolidate_bed_method = consolidate_bed_method
   self%bedload                = bedload
   self%bedload_method         = bedload_method
   self%add_to_density         = add_to_density

   ! do some printing
   LEVEL2 ' particle diameter        : ',real(self%diameter)
   LEVEL2 ' cohesive                 : ',self%cohesive
   LEVEL2 ' bottom stress method     : ',self%bottom_stress_method
   LEVEL2 ' correct EOS              : ',self%add_to_density
   select case (self%bottom_stress_method )
      case ( 0 )
         LEVEL2 ' constant critical bottom stress : ',real(self%tauc_const),' Pa'
      case ( 1 )
         LEVEL2 ' bottom stress method : Soulsby 1990 '
      case ( 2 )
         LEVEL2 ' bottom stress method : van Rijn 1984 '
   end select

   select case (self%sinking_method )
      case ( 0 )
         LEVEL2 ' sinking method : constant sinking speed : ',real(self%ws_const),' m/s'
      case ( 1 )
         if ( self%cohesive ) then
            LEVEL2 ' sinking method : Krone 1963 '
         else
            LEVEL2 ' sinking method : Soulsby 1997 '
         endif
      case ( 2 )
         if ( self%cohesive ) then
            LEVEL2 ' sinking method : Winterwerp 2001 '
         else
            LEVEL2 ' sinking method : Stokes/Newton '
         endif
      case ( 3 )
         if ( self%cohesive ) then
            LEVEL2 ' sinking method : Mehta 1986 '
         else
            LEVEL2 ' sinking method : currently not defined'
         endif
   end select

   if ( self%consolidate_bed ) then
      select case (self%consolidate_bed_method )
         case ( 1 )
            LEVEL2 ' bed consolidation method : low erodibility '
         case ( 2 )
            LEVEL2 ' bed consolidation method : medium erodibility '
         case ( 3 )
            LEVEL2 ' bed consolidation method : high erodibility '
      end select
      LEVEL3 ' Check the code for details'
   endif

   ! bed load is only computed for non-cohesive spm
   if ( self%bedload .and. self%cohesive ) then
      LEVEL2 'Bed load is only computed for non-cohesive spm!'
      self%bedload = .false.
   endif
   
   if ( self%bedload ) then
      select case ( self%bedload_method )
         case ( 0 )
             LEVEL2 ' bedload switched off '
         case ( 1 )
             LEVEL2 ' account for bedload : van Rijn 1984 '
         case ( 2 )
             LEVEL2 ' account for bedload : Nielsen 1992 '
         case ( 3 )
             LEVEL2 ' account for bedload : Engelund & Hansen 1972 '
      end select
   endif

   ! Register state variables
   call self%register_state_variable(self%id_spm,'spm','mg/l','concentration of SPM',     &
                                    c_init,minimum=0.0_rk, &
                                    vertical_movement=self%ws_const, &
                                    no_river_dilution=.true.)

   if ( self%pm_pool ) &
      call self%register_state_variable(self%id_pmpool,'pmpool','kg/m**2','mass/m**2 of PM in sediment',  &
         self%mass_sed_init)

   call self%register_horizontal_diagnostic_variable(self%id_massflux,'massflux','kg/m**2/s', &
      'massflux in the exchange layer', time_treatment=time_treatment_averaged)

   if ( self%bedload ) then
      call self%register_horizontal_diagnostic_variable(self%id_bedload, 'bedload', &
         'kg/m**2/s','massflux due to bed load', time_treatment=time_treatment_averaged)
   endif

   call self%set_variable_property(self%id_spm,'diameter',self%diameter)
   call self%set_variable_property(self%id_spm,'cohesive',self%cohesive)
   call self%set_variable_property(self%id_spm,'add_to_density',self%add_to_density)
   call self%set_variable_property(self%id_spm,'bedload',self%bedload)
   call self%set_variable_property(self%id_spm,'consolidate_bed',self%consolidate_bed)
   call self%set_variable_property(self%id_spm,'density',self%rho)
   call self%set_variable_property(self%id_spm,'spm',.true.)

   ! Register environmental dependencies
   call self%register_dependency(self%id_taub, standard_variables%bottom_stress)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_rhow, standard_variables%density)

   return

99 call self%fatal_error('iow_spm_create','Error reading namelist iow_spm')

100 call self%fatal_error('iow_spm_create','Namelist iow_spm was not found')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Sedimentation/Erosion
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
! Calculating the benthic fluxes
!
   implicit none

! !INPUT PARAMETERS:
   class (type_iow_spm), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

! !LOCAL VARIABLES:
   real(rk)                     :: taub,spm,pmpool
   real(rk)                     :: porosity
   real(rk), parameter          :: rho_0=1025.0_rk
   real(rk)                     :: Erosion_Flux,Sedimentation_Flux
   real(rk)                     :: tauc_erosion, tauc_sedimentation
   real(rk)                     :: temp
   real(rk)                     :: Ds
   real(rk)                     :: visc
   real(rk)                     :: theta
   real(rk)                     :: rhow
   real(rk)                     :: rhop
   real(rk)                     :: erodedmass
   real(rk)                     :: offset
   real(rk)                     :: massflux
   real(rk)                     :: scale_factor=100.0_rk
   real(rk)                     :: stressexponent=1.0_rk
   real(rk)                     :: qstar
   real(rk)                     :: bedload

   real(rk), parameter          :: g=9.81_rk
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Ulf Gräwe
!
!EOP
!-----------------------------------------------------------------------
!BOC

   if ( self%pm_pool ) then

      porosity=0.333_rk

      if ( self%cohesive ) then
         stressexponent = 1.0_rk
      else
         stressexponent = 1.5_rk
      endif

      _FABM_HORIZONTAL_LOOP_BEGIN_

      _GET_(self%id_spm,spm)
      _GET_(self%id_temp,temp)
      _GET_(self%id_rhow,rhow)
      _GET_HORIZONTAL_(self%id_pmpool,pmpool)
      _GET_HORIZONTAL_(self%id_taub,taub)

      select case ( self%bottom_stress_method )
         case ( 0 )
            tauc_erosion = self%tauc_const
         case ( 1 )
            ! Soulsby 1990
            rhop  = self%rho
            ! compute the dynamic viscosity
            visc  = 1.0_rk/rhow*1.9909e-6_rk*exp(1828.4_rk/(temp+273.15_rk))
            Ds    = (g/visc**2*(rhop/rhow-1.0_rk))**(0.3333_rk)*self%diameter
            theta = 0.3_rk/(1.0_rk + 1.2_rk*Ds) + 0.055_rk*(1.0_rk-exp(-0.02_rk*Ds));
            tauc_erosion = g*self%diameter*(rhop-rhow)*theta
         case ( 2 )
            ! van Rijn 1984
            rhop  = self%rho
            visc  = 1.0_rk/rhow*1.9909e-6_rk*exp(1828.4_rk/(temp+273.15_rk))
            Ds    = (g/visc**2*(rhop/rhow-1.0_rk))**(0.3333_rk)*self%diameter
            if ( Ds .le. 4    ) then
               tauc_erosion = 0.24_rk*Ds**(-0.9_rk)
            elseif ( (Ds .gt. 4  ) .and. (Ds .le. 10 ) ) then
               tauc_erosion = 0.14_rk*Ds**(-0.64_rk)
            elseif ( (Ds .gt. 10 ) .and. (Ds .le. 20 ) ) then
               tauc_erosion = 0.04_rk*Ds**(-0.1_rk)
            elseif ( (Ds .gt. 20 ) .and. (Ds .le. 150) ) then
               tauc_erosion = 0.013_rk*Ds**(0.29_rk)
            else
               tauc_erosion = 0.056
            endif
            tauc_erosion = tauc_erosion*g*self%diameter*(rhop-rhow)
      end select

      if ( self%consolidate_bed ) then
         ! mimic a consolidation of the bed : Dickhudt et. al. 2009
         erodedmass = (max(self%mass_sed_init-pmpool,0.0_rk))
         ! convert to kg/m2
         erodedmass = erodedmass/scale_factor
         select case ( self%consolidate_bed_method )
            case ( 1 )
               ! low erodebility
               offset = 0.835_rk*erodedmass**0.508_rk
            case ( 2 )
               ! transitional erodebility
               offset = 0.531_rk*erodedmass**0.646_rk
            case ( 3 )
               ! high erodebility
               offset = 0.243_rk*erodedmass**0.754_rk
            case default
               offset = 0.0_rk
         end select
         tauc_erosion = tauc_erosion + offset
      endif

      tauc_sedimentation = tauc_erosion*self%tauc_factor
      
      Erosion_Flux       = 0.0_rk
      Sedimentation_Flux = 0.0_rk
      
      ! 1-spm_porosity is the fractional bed concentration
      if( (pmpool .gt. 0.1_rk) .and. (taub .gt. tauc_erosion) ) then
         ! if there are sediments in the pool
         Erosion_Flux = self%erosion_const / rho_0                  &
                       *(1.0_rk-porosity) * (taub-tauc_erosion)**stressexponent
         Erosion_Flux = max(Erosion_Flux,0.0_rk)
      endif

      !sedimentation flux:
      Sedimentation_Flux = min(0.0_rk,ws(self,temp,rhow,spm) * spm *(1.0_rk-taub / tauc_sedimentation))

      if ( self%bedload ) then
         rhop  = self%rho
         visc  = 1.0_rk/rhow*1.9909e-6_rk*exp(1828.4_rk/(temp+273.15_rk))
         Ds    = (g/visc**2*(rhop/rhow-1.0_rk))**(0.3333_rk)*self%diameter
         select case ( self%bedload_method )
            case ( 0 )
               bedload = 0.0_rk
            case ( 1 )
               ! van Rijn 1984
               if( taub .gt. tauc_erosion ) then
                  qstar   = 0.053_rk*Ds**(-0.3_rk)*(taub/tauc_erosion-1.0_rk)**(2.1)
                  bedload = qstar*sqrt(g*self%diameter**3*(rhop/rhow-1.0_rk))
               endif
            case ( 2 )
               ! Nielsen 1992
               if( taub .gt. tauc_erosion ) then
                  qstar   = 12.0_rk*(taub-tauc_erosion)/(g*self%diameter*(rhop-rhow))* &
                           sqrt(taub/g*self%diameter*(rhop-rhow))
                  bedload = qstar*sqrt(g*self%diameter**3*(rhop/rhow-1.0_rk))
               endif
            case ( 3 )
               ! Engelund & Hansen 1972
               qstar   = 0.05_rk*(taub/(g*self%diameter*(rhop-rhow)))**2.5
               bedload = qstar*sqrt(g*self%diameter**3*(rhop/rhow-1.0_rk))
            case ( 4 )
               ! Lesser et al. 2004
               bedload = 0.5_rk*rhop*self%diameter*sqrt(taub/rhow)*Ds**(-0.3_rk)*taub/(g*self%diameter*(rhop-rhow))
               bedload = bedload*sqrt(g*self%diameter**2*(rhop/rhow-1.0_rk))
         end select
      else
         bedload = 0.0_rk
      endif

      ! compute mass flux
      massflux = Sedimentation_Flux+Erosion_Flux
      ! unit is g/m**3 * m/s
      _SET_BOTTOM_EXCHANGE_(self%id_spm,massflux)
      ! unit is kg/m**2/s
      _SET_ODE_BEN_(self%id_pmpool,-(massflux)/1000.0_rk)
      ! Export diagnostic variables mass flux - unit is kg/m**2/s
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_massflux,-(massflux)/1000.0_rk)
      ! Export diagnostic variables bed load flux - unit is kg/m**2/s
      if ( self%bedload ) then
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_bedload,bedload)
      endif

      _FABM_HORIZONTAL_LOOP_END_

   end if

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: iow_spm_get_vertical_movement
!
! !INTERFACE:
   subroutine get_vertical_movement(self,_FABM_ARGS_GET_VERTICAL_MOVEMENT_)

   implicit none

! !INPUT PARAMETERS:
   class (type_iow_spm), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_VERTICAL_MOVEMENT_

! !LOCAL VARIABLES
   real(rk)      :: temp
   real(rk)      :: rhow
   real(rk)      :: spm
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Ulf Gräwe
!
   _FABM_LOOP_BEGIN_

   _GET_DEPENDENCY_(self%id_temp,temp)
   _GET_DEPENDENCY_(self%id_rhow,rhow)
   _GET_(self%id_spm,spm)
   _SET_VERTICAL_MOVEMENT_(self%id_spm,ws(self,temp,rhow,spm))

   _FABM_LOOP_END_

   end subroutine get_vertical_movement

!EOP
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: iow_spm_get_light_extinction
!
! !INTERFACE:

   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)

   implicit none

! !INPUT PARAMETERS:

   class (type_iow_spm), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_

! !LOCAL VARIABLES
   real(rk)      :: spm
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   _GET_(self%id_spm,spm)
   _SET_EXTINCTION_(self%shading*spm)

   _LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------

   function ws(self,temp,rhow,spm) result(svs)

   implicit none

   class (type_iow_spm), intent(in) :: self
   real(rk),             intent(in) :: temp
   real(rk),             intent(in) :: rhow     ! water density
   real(rk),             intent(in) :: spm      ! spm concentration

   real(rk)                        :: svs       ! sinking velocity in m/s

   real(rk)                        :: rhop      ! dry bed density of particle
   real(rk)                        :: visc      ! kinematic viscosity, visc=visc(T)
   real(rk)                        :: Ds        !
   real(rk)                        :: pRe       ! particle Reynolds number
   real(rk)                        :: Cd        ! dynamic drag coefficient
   real(rk)                        :: sF        ! shape factor
   real(rk)                        :: rrho      ! reduced density
   real(rk)                        :: phi       ! volumetric concentration
   real(rk)                        :: phistar   ! phi-limiter
   real(rk)                        :: a,b       ! parameter for cohesive sinking velocity
   integer                         :: i

   real(rk), parameter             :: g=9.81_rk

   if ( self%cohesive ) then
      ! start with the cohesive sinking computation
      ! background sinking is based on Krone 1963
      ! spm must be in kg/m**3 and svs is in mm/s
      a   = 0.36_rk
      b   = 1.33_rk
      svs = (a*(spm/1000.0_rk)**b)/1000.0_rk

      select case ( self%sinking_method )
         case ( 0 )
            svs = abs(self%ws_const)
         case ( 1 )
            svs = svs
         case ( 2 )
            ! Winterwerp 2001
            ! compute volumetric concentration
            phi     = min(0.9999999_rk,spm/g/rhop)
            ! do a limiter
            phistar = min(1.0_rk,phi)
            svs     = svs*(1.0_rk-phistar)*(1.0_rk-phi)/(1.0_rk+2.5_rk*phi)
         case ( 3 )
            ! Mehta 1986
            svs  = svs*(max(1.0_rk-0.008*spm/1000.0_rk,0.0_rk))**5

      end select
      svs = max(svs,0.01_rk)/1000.0_rk
      svs = -svs
   else
      ! now do the non-cohesive sinking
      select case ( self%sinking_method )
         case ( 0 )
            svs = -abs(self%ws_const)
         case ( 1 )
            ! Soulsby R (1997) Dynamics of marine sands - a manual for practical
            !   applications. Thomas Telford, London
            rhop = self%rho
            visc = 1.0_rk/rhow*1.9909e-6_rk*exp(1.8284e3_rk/(temp+273.15_rk))
            Ds   = (g/visc**2*(rhop/rhow-1.0_rk))**(0.3333_rk)*self%diameter
            svs  = visc/self%diameter*(sqrt(10.36_rk**2+1.049_rk*Ds**3)-10.36_rk)
            svs = -svs

         case ( 2 )
            ! Stokes / Newton
            sF   = 0.64_rk ! assume spherical particles
            rhop = self%rho
            rrho = max((rhop-rhow)/rhow,0.0_rk)
            visc = 1.0_rk/rhow*1.9909e-6_rk*exp(1.8284e3_rk/(temp+273.15_rk))
            ! at first, estimate Stokes terminal sinking velocity
            svs  = self%diameter**2*g*rrho/18.0_rk/visc
            ! now, do the iteration since the drag depends on the sinking velocity
            ! UG: I think that 7 iterations are enough. One could also do a while loop
            !     and test for the convergence
            do i=1,7
               pRe = sF*svs*self%diameter/visc
               ! The Cd formula is only valid for 1<pRe<1000, which is mostly the case.
               Cd  = 18.5_rk/pRe**0.6_rk;
               svs = sqrt(4.0_rk*g*self%diameter*rrho/3.0_rk/Cd);
            enddo
            svs = -svs

      end select

   endif

   end function ws


   end module fabm_iow_spm
