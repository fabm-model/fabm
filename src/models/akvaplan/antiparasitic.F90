#include "fabm_driver.h"

module akvaplan_antiparasitic

   ! Model for a dissolved compound ("antiparasitic") that partially absorbs onto particulate organic matter (POM).
   ! Partitioning between dissolved and adsorbed phases is assumed to be instantaneous and
   ! controlled by the organic carbon : water partition coefficient.
   ! The compound experiences temperature-dependent degradation at phase-dependent rates.
   ! The model can contain several classes of POM with different sinking velocities.
   ! At the bed, both the compound and POM are incorporated into the sediment; from where they may be resuspended again.
   ! Implementation by Jorn Bruggeman and Ricardo Torres, Plymouth Marine Laboratory

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_antiparasitic
      ! The model's own variables
      type (type_state_variable_id),        allocatable :: id_pom(:)
      type (type_bottom_state_variable_id), allocatable :: id_pom_bot(:)
      type (type_state_variable_id)                     :: id_antiparasitic
      type (type_bottom_state_variable_id)              :: id_antiparasitic_bot

      ! Environmental dependencies
      type (type_dependency_id)        :: id_T     ! temperature
      type (type_dependency_id)        :: id_h     ! cell thickness
      type (type_bottom_dependency_id) :: id_shear ! bottom shear

      ! Identifiers for diagnostic variables
      type (type_diagnostic_variable_id)        :: id_degradation_pom_flux, id_degradation_c_flux          ! pelagic degradation flux for pom and treatment (antiparasitic!)
      type (type_bottom_diagnostic_variable_id) :: id_resuspension_flux_pom, id_deposition_flux_pom        ! deposition and resuspension fluxes of all pom 
      type (type_bottom_diagnostic_variable_id) :: id_resuspension_flux_c, id_deposition_flux_c            ! deposition and resuspension fluxes of antiparasitic treatment 
      type (type_bottom_diagnostic_variable_id) :: id_degradation_pom_bot_flux, id_degradation_c_bot_flux  ! sediment degradation fluxes for pom and treatment (antiparasitic!)

      ! Parameters
      real(rk), allocatable :: w(:)
      real(rk), allocatable :: tau_bot_crit(:)
      real(rk) :: bed_por
      real(rk) :: erate
      real(rk) :: K_oc
      real(rk) :: E_a
      real(rk) :: T_ref
      real(rk) :: k_w
      real(rk) :: k_ads
      real(rk) :: k_pom
      real(rk) :: max_dt    ! maximum time step (s) to be expected - used to limit benthic-pelagic fluxes for numerical stability
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: get_vertical_movement
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: R = 8.3144598_rk    ! universal gas constant (J mol-1 K-1)
   real(rk), parameter :: Kelvin = 273.15_rk  ! offset between degrees Celsius and Kelvin

contains

   subroutine initialize(self, configunit)
      class (type_antiparasitic), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      integer                         :: i, npom
      character(len=8)                :: index
      logical                         :: conserved
      character(len=attribute_length) :: standard_name

      ! Obtain parameter values (convert all rate constants to s-1 to match FABM's internal time unit)
      call self%get_parameter(npom, 'npom', '', 'number of particulate organic matter classes')
      allocate(self%id_pom(npom), self%id_pom_bot(npom), self%w(npom), self%tau_bot_crit(npom))
      do i=1,npom
         write (index,'(i0)') i
         call self%get_parameter(self%w(i), 'w'//trim(index), 'm d-1', 'sinking velocity of particulate organic matter class '//trim(index), scale_factor=1.0_rk/secs_pr_day)
         call self%register_state_variable(self%id_pom(i), 'pom'//trim(index), 'g C m-3', 'concentration of particulate organic matter class '//trim(index), initial_value=0.0_rk, vertical_movement=-self%w(i), minimum=0.0_rk)
         call self%register_state_variable(self%id_pom_bot(i), 'pom'//trim(index)//'_bot', 'g C m-2', 'particulate organic matter class '//trim(index)//' in sediment', initial_value=0.0_rk, minimum=0.0_rk)
      end do
      call self%get_parameter(self%K_oc, 'K_oc', 'm3 g-1', 'organic carbon : water partition coefficient')

      ! Degradation of antiparasitic and POM
      call self%get_parameter(self%E_a, 'E_a', 'J mol-1', 'activation energy for antiparasitic decay (0: no temperature dependence)', minimum=0.0_rk, default=0.0_rk)
      call self%get_parameter(self%T_ref, 'T_ref', 'degree Celsius', 'reference temperature', default=20.0_rk)
      call self%get_parameter(self%k_ads, 'k_ads', 'd-1', 'degradation of adsorbed antiparasitic at reference temperature', default=0.0_rk, scale_factor=1.0_rk/secs_pr_day, minimum=0._rk)
      call self%get_parameter(self%k_w, 'k_w', 'd-1', 'degradation of dissolved antiparasitic at reference temperature', default=0.0_rk, scale_factor=1.0_rk/secs_pr_day, minimum=0._rk)
      call self%get_parameter(self%k_pom, 'k_pom', 'd-1', 'degradation of particulate organic matter', default=0.0_rk, scale_factor=1.0_rk/secs_pr_day, minimum=0._rk)

      ! Resuspension of POM
      call self%get_parameter(self%bed_por, 'bed_por', '-', 'bed porosity')
      call self%get_parameter(self%erate, 'erate', 'g m-2 d-1', 'bed erodibility', default=0.0_rk, scale_factor=1.0_rk/secs_pr_day, minimum=0._rk)
      do i=1,npom
         write (index,'(i0)') i
         call self%get_parameter(self%tau_bot_crit(i), 'tau_bot_crit'//trim(index), 'Pa', 'critical shear stress for particulate organic matter class '//trim(index), default=0.0_rk, scale_factor=1.0_rk/secs_pr_day, minimum=0.0_rk)
      end do

      ! Numerics and debugging
      call self%get_parameter(self%max_dt, 'max_dt', 's', 'maximum time step used to limit vertical exchanges', default=900.0_rk)
      call self%get_parameter(conserved, 'conserved', '', 'treat antiparasitic and POM as conserved quantities (activates budget tracking in host)', default=.false.)

      ! Register state variables for the antiparasitic compound.
      call self%register_state_variable(self%id_antiparasitic,'c','g m-3','total antiparasitic concentration (dissolved + adbsorbed)', initial_value=0.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_antiparasitic_bot,'c_bot','g m-2','total antiparasitic in sediment', minimum=0.0_rk)

      if (conserved) then
         standard_name = get_safe_name(trim(self%get_path())//'_total_ab')
         call self%add_to_aggregate_variable(type_bulk_standard_variable(name=standard_name(2:), units='g m-3', aggregate_variable=.true.,conserved=.true.), self%id_antiparasitic)
         call self%add_to_aggregate_variable(type_bulk_standard_variable(name=standard_name(2:), units='g m-3', aggregate_variable=.true.,conserved=.true.), self%id_antiparasitic_bot)

         standard_name = get_safe_name(trim(self%get_path())//'_total_pom')
         do i=1,npom
           call self%add_to_aggregate_variable(type_bulk_standard_variable(name=standard_name(2:), units='g m-3', aggregate_variable=.true.,conserved=.true.), self%id_pom(i))
           call self%add_to_aggregate_variable(type_bulk_standard_variable(name=standard_name(2:), units='g m-3', aggregate_variable=.true.,conserved=.true.), self%id_pom_bot(i))
         end do
     end if

      ! Register environmental dependencies.
      call self%register_dependency(self%id_T, standard_variables%temperature)
      call self%register_dependency(self%id_shear, standard_variables%bottom_stress)
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_degradation_pom_flux, 'pom_degradation', 'g C m-3 s-1', 'degradation of total particulate organic matter in pelagic')
      call self%register_diagnostic_variable(self%id_degradation_c_flux, 'c_degradation', 'g m-3 s-1', 'degradation of total antiparasitic in pelagic')

      call self%register_diagnostic_variable(self%id_resuspension_flux_pom, 'resuspension_flux_pom', 'g C m-2 s-1', 'total organic matter resuspension')
      call self%register_diagnostic_variable(self%id_deposition_flux_pom, 'deposition_flux_pom', 'g C m-2 s-1', 'total organic matter deposition')
      call self%register_diagnostic_variable(self%id_resuspension_flux_c, 'resuspension_flux_c', 'g m-2 s-1', 'antiparasitic resuspension')
      call self%register_diagnostic_variable(self%id_deposition_flux_c, 'deposition_flux_c', 'g m-2 s-1', 'antiparasitic deposition')

      call self%register_diagnostic_variable(self%id_degradation_pom_bot_flux, 'pom_bot_degradation','g C m-2 s-1', 'degradation of total organic matter in sediment')
      call self%register_diagnostic_variable(self%id_degradation_c_bot_flux, 'c_bot_degradation','g m-2 s-1', 'degradation of total antiparasitic in sediment')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_antiparasitic), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      integer :: i
      real(rk) :: pom(size(self%id_pom)), pom_tot, frac_ads, frac_wat, T, f, c, k_wat, k_ads

      _LOOP_BEGIN_

         do i=1,size(self%id_pom)
            _GET_(self%id_pom(i), pom(i))
         end do
         pom_tot = sum(pom)
         frac_wat = 1.0_rk / (self%K_oc * pom_tot + 1.0_rk)
         frac_ads = 1.0_rk - frac_wat

         ! Specific degradation (s-1) for adsorbed and dissolved antiparasitic
         _GET_(self%id_T, T)
         f = exp(self%E_a / R * (1.0_rk / (Kelvin + self%T_ref) - 1.0_rk / (Kelvin + T)))
         k_ads = self%k_ads * f
         k_wat = self%k_w * f
         _GET_(self%id_antiparasitic, c)
         _ADD_SOURCE_(self%id_antiparasitic, -(frac_ads * k_ads + frac_wat * k_wat) * c)

         ! Degradation (s-1) for POM
         do i=1,size(self%id_pom)
            _ADD_SOURCE_(self%id_pom(i), -self%k_pom * pom(i))
         end do

         ! Save diagnostics for total pom and antiparasitic degradation  
         _SET_DIAGNOSTIC_(self%id_degradation_pom_flux, self%k_pom * pom_tot) 
         _SET_DIAGNOSTIC_(self%id_degradation_c_flux, (frac_ads * k_ads + frac_wat * k_wat) * c) 

      _LOOP_END_

   end subroutine do

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_antiparasitic), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: i
      real(rk) :: h, pom(size(self%id_pom)), pom_bot(size(self%id_pom)), pom_tot, pom_bot_tot, frac_ads, w, c, tau_bot, erosion(size(self%id_pom)), rel_antiparasitic_erosion, T, f, k_ads, c_bot
      real(rk) :: pom_deposition

      _BOTTOM_LOOP_BEGIN_

         _GET_(self%id_h, h)
         do i = 1, size(self%id_pom)
            _GET_(self%id_pom(i), pom(i))
            _GET_BOTTOM_(self%id_pom_bot(i), pom_bot(i))
         end do
         pom_tot = sum(pom)
         pom_bot_tot = sum(pom_bot)

         ! Erosion rate (g m-2 s-1) per POM class
         _GET_BOTTOM_(self%id_shear, tau_bot)
         erosion = min(pom_bot/self%max_dt, self%erate * (1.0_rk - self%bed_por) * max(tau_bot / self%tau_bot_crit - 1.0_rk, 0.0_rk))

         ! Calculate relative erosion for antiparasitic substance.
         if (pom_bot_tot > 0) then
            rel_antiparasitic_erosion = sum(erosion)/pom_bot_tot
         else
            rel_antiparasitic_erosion = 0
         end if

         ! Sedimentation, resuspension, degradation per POM class
         do i = 1, size(self%id_pom)
           ! From pelagic to bottom. Convention is negative out of pelagic, positive into water
            _ADD_BOTTOM_FLUX_(self%id_pom(i), -min(h/self%max_dt, self%w(i)) * pom(i) + erosion(i))
           ! From bottom to pelagic. Convention is positive into bottom, negative into water. Added degradation of pom_bot
            _ADD_BOTTOM_SOURCE_(self%id_pom_bot(i), min(h/self%max_dt, self%w(i)) * pom(i) - erosion(i) - self%k_pom * pom_bot(i))
         end do

         ! Sinking rate of antiparasitic is weighted average of sinking rate per antiparasitic fraction
         ! (those adsorbed to individual POM classes, plus one non-sinking fraction for the dissolved phase)
         frac_ads = 1.0_rk - 1.0_rk/(self%K_oc * pom_tot + 1.0_rk)
         if (pom_tot > 0) then
            w = min(h/self%max_dt, frac_ads * sum(pom * self%w) / pom_tot)
         else
            w = 0
         end if

         ! Specific degradation (s-1) for adsorbed antiparasitic in sediment
         ! (assume dissolved fraction in pore water is negligible)
         _GET_(self%id_T, T)
         f = exp(self%E_a / R * (1.0_rk / (Kelvin + self%T_ref) - 1.0_rk / (Kelvin + T)))
         k_ads = self%k_ads * f

         ! Sedimentation and resuspension for antiparasitic
         _GET_(self%id_antiparasitic, c)
         _GET_BOTTOM_(self%id_antiparasitic_bot, c_bot)
         _ADD_BOTTOM_FLUX_(self%id_antiparasitic, -w * c + rel_antiparasitic_erosion * c_bot)
         _ADD_BOTTOM_SOURCE_(self%id_antiparasitic_bot, w * c - rel_antiparasitic_erosion * c_bot - k_ads * c_bot)

         ! Save diagnostic total pom bottom degradation  
         _SET_BOTTOM_DIAGNOSTIC_(self%id_degradation_pom_bot_flux, self%k_pom * pom_bot_tot) 
         ! Save diagnostic antiparasitic bottom degradation  
         _SET_BOTTOM_DIAGNOSTIC_(self%id_degradation_c_bot_flux, k_ads * c_bot ) 

         ! Save diagnostic total pom bottom deposition
         pom_deposition = 0.0_rk
         do i = 1, size(self%id_pom)
            pom_deposition = pom_deposition + min(h/self%max_dt, self%w(i)) * pom(i) 
         end do
         _SET_BOTTOM_DIAGNOSTIC_(self%id_deposition_flux_pom, pom_deposition)

         ! Save diagnostic total pom bottom resuspension
         _SET_BOTTOM_DIAGNOSTIC_(self%id_resuspension_flux_pom, sum(erosion))

         ! Save diagnostic antiparasitic bottom deposition  
         _SET_BOTTOM_DIAGNOSTIC_(self%id_deposition_flux_c, w * c )

         ! Save diagnostic antiparasitic bottom resuspension  
         _SET_BOTTOM_DIAGNOSTIC_(self%id_resuspension_flux_c, rel_antiparasitic_erosion * c_bot) 

       _BOTTOM_LOOP_END_

   end subroutine do_bottom

   subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_antiparasitic), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      integer :: i
      real(rk) :: pom(size(self%id_pom)), pom_tot, frac_ads, w

      _LOOP_BEGIN_

         ! Sinking rate of antiparasitic is weighted average of sinking rate per antiparasitic fraction
         ! (those adsorbed to individual POM classes, plus one non-sinking fraction for the dissolved phase)
         ! NB we do not need to set sinking rate per POM class - that is handled by the vertical_movement argument specified during initialization
         do i = 1, size(self%id_pom)
            _GET_(self%id_pom(i), pom(i))
         end do
         pom_tot = sum(pom)
         frac_ads = 1.0_rk - 1.0_rk/(self%K_oc * pom_tot + 1.0_rk)
         if (pom_tot > 0) then
            w = frac_ads * sum(pom * self%w) / pom_tot
         else
            w = 0
         end if
         _ADD_VERTICAL_VELOCITY_(self%id_antiparasitic, -w)

      _LOOP_END_

   end subroutine get_vertical_movement

end module
