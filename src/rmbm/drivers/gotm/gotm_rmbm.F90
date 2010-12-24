!$Id$
#include "cppdefs.h"
#include "rmbm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: gotm_rmbm --- Interface to Repository of Marine Biogeochemical Models (RMBM)
!
! !INTERFACE:
   module gotm_rmbm
!
! !DESCRIPTION:
! TODO
! 
! !USES:
   use rmbm
   use rmbm_types
   
   implicit none

   interface   
      subroutine ode_solver(solver,numc,nlev,dt,cc,right_hand_side_rhs,right_hand_side_ppdd)
         integer,  intent(in)                :: solver,nlev,numc
         REALTYPE, intent(in)                :: dt
         REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)

         interface
            subroutine right_hand_side_ppdd(first,numc,nlev,cc,pp,dd)
               logical,  intent(in)                 :: first
               integer,  intent(in)                 :: numc,nlev
               REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
               
               REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
               REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
            end
         end interface

         interface
            subroutine right_hand_side_rhs(first,numc,nlev,cc,rhs)
               logical,  intent(in)                 :: first
               integer,  intent(in)                 :: numc,nlev
               REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
               REALTYPE, intent(out)                :: rhs(1:numc,0:nlev)
            end
         end interface
      end subroutine
   end interface
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_gotm_rmbm,init_var_gotm_rmbm
   public set_env_gotm_rmbm,do_gotm_rmbm
   public clean_gotm_rmbm,save_gotm_rmbm

!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
   ! Namelist variables
   REALTYPE                  :: cnpar
   integer                   :: w_adv_method,w_adv_discr,ode_method,split_factor
   logical                   :: rmbm_calc,bioshade_feedback,repair_state

   ! Model
   type (type_model), pointer :: model
   
   ! Arrays for state and diagnostic variables
   REALTYPE,allocatable,dimension(LOCATION_DIMENSIONS,:),target :: cc
   REALTYPE,allocatable,dimension(LOCATION_DIMENSIONS,:)        :: cc_diag

   ! Arrays for work, vertical movement, and cross-boundary fluxes
   REALTYPE,allocatable,dimension(LOCATION_DIMENSIONS,:) :: ws
   REALTYPE,allocatable,dimension(:)                     :: sfl,bfl,total,cc_diag_hz
   REALTYPE,allocatable _ATTR_DIMENSIONS_1_                :: local
   
   ! Arrays for environmental variables not supplied externally.
   REALTYPE,allocatable,dimension(LOCATION_DIMENSIONS)   :: par,pres
   
   ! External variables
   REALTYPE :: dt,dt_eff   ! External and internal time steps
   integer  :: w_adv_ctr   ! Scheme for vertical advection (0 if not used)
   REALTYPE,pointer,dimension(LOCATION_DIMENSIONS) :: nuh,h,bioshade,rad,w,z
   REALTYPE,pointer _ATTR_LOCATION_DIMENSIONS_HZ_ :: precip,evap

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_gotm_rmbm(namlst,fname)
!
! !DESCRIPTION: 
! Initializes the GOTM-RMBM driver module by reading settings from rmbm.nml.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)        :: namlst
   character(len=*), intent(in)        :: fname

!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!  local variables
   integer                   :: i
   character(len=64)         :: models(256)
   type (type_model),pointer :: childmodel
   namelist /bio_nml/ rmbm_calc,models,                                 &
                      cnpar,w_adv_discr,ode_method,split_factor,        &
                      bioshade_feedback,repair_state
!
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_gotm_rmbm'
   
   ! Initialize RMBM model identifiers to invalid id.
   rmbm_calc         = .false.
   models            = ''
   cnpar             = _ONE_
   w_adv_discr       = 6
   ode_method        = 1
   split_factor      = 1
   bioshade_feedback = .true.
   repair_state      = .false.

   ! Open the namelist file and read the namelist.
   ! Note that the namelist file is left open until the routine terminates,
   ! so RMBM can read more namelists from it during initialization.
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_nml,err=99)

   if (rmbm_calc) then
   
      ! Create model tree
      model => rmbm_create_model()
      do i=1,ubound(models,1)
         if (trim(models(i)).ne.'') &
            childmodel => rmbm_create_model(trim(models(i)),parent=model)
      end do
      
      ! Initialize model tree (creates metadata and assigns variable identifiers)
      call rmbm_init(model,namlst)

      ! Report prognostic variable descriptions
      LEVEL2 'RMBM pelagic state variables:'
      do i=1,ubound(model%info%state_variables,1)
         LEVEL3 trim(model%info%state_variables(i)%name), '  ', &
                trim(model%info%state_variables(i)%units),'  ',&
                trim(model%info%state_variables(i)%longname)
      end do

      LEVEL2 'RMBM benthic state variables:'
      do i=1,ubound(model%info%state_variables_ben,1)
         LEVEL3 trim(model%info%state_variables_ben(i)%name), '  ', &
                trim(model%info%state_variables_ben(i)%units),'  ',&
                trim(model%info%state_variables_ben(i)%longname)
      end do

      ! Report diagnostic variable descriptions
      LEVEL2 'RMBM diagnostic variables defined on the full model domain:'
      do i=1,ubound(model%info%diagnostic_variables,1)
         LEVEL3 trim(model%info%diagnostic_variables(i)%name), '  ', &
                trim(model%info%diagnostic_variables(i)%units),'  ',&
                trim(model%info%diagnostic_variables(i)%longname)
      end do

      LEVEL2 'RMBM diagnostic variables defined on a horizontal slice of the model domain:'
      do i=1,ubound(model%info%diagnostic_variables_hz,1)
         LEVEL3 trim(model%info%diagnostic_variables_hz(i)%name), '  ', &
                trim(model%info%diagnostic_variables_hz(i)%units),'  ',&
                trim(model%info%diagnostic_variables_hz(i)%longname)
      end do

      ! Report type of solver 
      LEVEL2 "Using Eulerian solver"
      select case (ode_method)
         case (1)
            LEVEL2 'Using euler_forward()'
         case (2)
            LEVEL2 'Using runge_kutta_2()'
         case (3)
            LEVEL2 'Using runge_kutta_4()'
         case (4)
            LEVEL2 'Using patankar()'
         case (5)
            LEVEL2 'Using patankar_runge_kutta_2()'
         case (6)
            LEVEL2 'Using patankar_runge_kutta_4()'
         case (7)
            LEVEL2 'Using modified_patankar()'
         case (8)
            LEVEL2 'Using modified_patankar_2()'
         case (9)
            LEVEL2 'Using modified_patankar_4()'
         case (10)
            LEVEL2 'Using emp_1()'
         case (11)
            LEVEL2 'Using emp_2()'
         case (1003)
            LEVEL2 'Using runge_kutta_4() with pp/dd matrices'
         case default
            stop "init_gotm_rmbm: no valid ode_method specified in rmbm.nml!"
      end select
      
      ! Initialize RMBM output (creates NetCDF variables)
      call init_output_gotm_rmbm()

   end if

   ! Close the namelist file
   close(namlst)

   return

98 LEVEL2 'I could not open '//trim(fname)
   LEVEL2 'If thats not what you want you have to supply '//trim(fname)
   LEVEL2 'See the bio example on www.gotm.net for a working '//trim(fname)
   rmbm_calc = .false.
   return
99 FATAL 'I could not read '//trim(fname)
   stop 'init_gotm_rmbm'
   end subroutine init_gotm_rmbm
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise bio variables
!
! !INTERFACE:
   subroutine init_var_gotm_rmbm(LOCATION)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   IMPLICIT NONE
   
   _LOCATION_TYPE_,intent(in) :: LOCATION
!
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i,rc
!
!-----------------------------------------------------------------------
!BOC
   if (.not. rmbm_calc) return
   
   call rmbm_set_domain(model,LOCATION)

   ! Allocate state variable array for pelagic amnd benthos combined and provide initial values.
   ! In terms of memory use, it is a waste to allocate storage for benthic variables across the entire
   ! column (the bottom layer should suffice). However, it is important that all values at a given point
   ! in time are integrated simultaneously in multi-step algorithms. This currently can only be arranged
   ! By storing benthic values together with the pelagic, in a full;y depth=-explicit array.
   allocate(cc(1:ubound(model%info%state_variables,1)+ubound(model%info%state_variables_ben,1),LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (cc)'
   cc = _ZERO_
   do i=1,ubound(model%info%state_variables,1)
      cc(i,:) = model%info%state_variables(i)%initial_value
   end do
   do i=1,ubound(model%info%state_variables_ben,1)
      cc(ubound(model%info%state_variables,1)+i,1) = model%info%state_variables_ben(i)%initial_value
   end do

   ! Allocate diagnostic variable array and set all values to zero.
   ! (needed because time-integrated/averaged variables will increment rather than set the array)
   allocate(cc_diag(1:ubound(model%info%diagnostic_variables,1),LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag)'
   cc_diag = _ZERO_

   ! Allocate diagnostic variable array and set all values to zero.
   ! (needed because time-integrated/averaged variables will increment rather than set the array)
   allocate(cc_diag_hz(1:ubound(model%info%diagnostic_variables_hz,1)),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag_hz)'
   cc_diag_hz = _ZERO_

   ! Allocate array with vertical movement rates (m/s, positive for upwards),
   ! and set these to the values provided by the model.
   allocate(ws(LOCATION_RANGE,1:ubound(model%info%state_variables,1)),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   do i=1,ubound(model%info%state_variables,1)
      ws(:,i) = model%info%state_variables(i)%vertical_movement
   end do

   ! Allocate array for surface fluxes and initialize these to zero (no flux).
   allocate(sfl(1:ubound(model%info%state_variables,1)),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (sfl)'
   sfl = _ZERO_

   ! Allocate array for bottom fluxes and initialize these to zero (no flux).
   allocate(bfl(1:ubound(model%info%state_variables,1)),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (bfl)'
   bfl = _ZERO_

   ! Allocate array for photosynthetically active radiation (PAR).
   ! This will be calculated internally during each time step.
   allocate(par(LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (par)'
   call rmbm_link_data(model,varname_par,par(1:LOCATION))

   ! Allocate array for local pressure.
   ! This will be calculated [approximated] from layer depths internally during each time step.
   allocate(pres(LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (pres)'
   call rmbm_link_data(model,varname_pres,pres(1:LOCATION))

   ! Allocate arrays for storing local and column-integrated values of diagnostic variables.
   ! These are used during each save.
   allocate(total(1:ubound(model%info%conserved_quantities,1)),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (total)'
#ifdef _RMBM_USE_1D_LOOP_
   allocate(local(1:LOCATION,1:ubound(model%info%conserved_quantities,1)),stat=rc)
#else
   allocate(local(1:ubound(model%info%conserved_quantities,1)),stat=rc)
#endif
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (local)'

   end subroutine init_var_gotm_rmbm
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set bio module environment 
!
! !INTERFACE: 
   subroutine set_env_gotm_rmbm(dt_,w_adv_method_,w_adv_ctr_,temp,salt_,rho,nuh_,h_,w_,rad_,bioshade_,I_0,wnd,precip_,evap_,z_)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: dt_
   integer,  intent(in) :: w_adv_method_,w_adv_ctr_
   REALTYPE, intent(in),target _ATTR_LOCATION_DIMENSIONS_    :: temp,salt_,rho,nuh_,h_,w_,rad_,bioshade_,z_
   REALTYPE, intent(in),target _ATTR_LOCATION_DIMENSIONS_HZ_ :: I_0,wnd,precip_,evap_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------!
!BOC
   if (.not. rmbm_calc) return

   ! Provide pointers to arrays with environmental variables to RMBM.
   call rmbm_link_data   (model,varname_temp,   temp)
   call rmbm_link_data   (model,varname_salt,   salt_)
   call rmbm_link_data   (model,varname_dens,   rho)
   call rmbm_link_data_hz(model,varname_wind_sf,wnd)
   call rmbm_link_data_hz(model,varname_par_sf, I_0)
   
   ! Save pointers to external dynamic variables that we need later (in do_gotm_rmbm)
   nuh => nuh_             ! turbulent heat diffusivity [1d array] used to diffuse biogeochemical state variables
   h   => h_               ! layer heights [1d array] needed for advection, diffusion
   w   => w_               ! vertical medium velocity [1d array] needed for advection of biogeochemical state variables
   rad => rad_             ! short wave radiation [1d array] used to calculate photosynthetically active radiation
   bioshade => bioshade_   ! biogeochemical light attenuation coefficients [1d array], output of biogeochemistry, input for physics
   z => z_                 ! depth [1d array], used to calculate local pressure
   precip => precip_       ! precipitation [scalar] - used to calculate dilution of biogeochemical variables due to increased water volume
   evap   => evap_         ! evaporation [scalar] - used to calculate concentration of biogeochemical variables due to decreased water volume
   
   ! Copy scalars that will not change during simulation, and are needed in do_gotm_rmbm)
   dt = dt_
   w_adv_method = w_adv_method_
   w_adv_ctr = w_adv_ctr_

   ! Calculate and save internal time step.
   dt_eff = dt/float(split_factor)
   
   end subroutine set_env_gotm_rmbm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update the RMBM model
!
! !INTERFACE:
   subroutine do_gotm_rmbm(nlev)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   use util,only: flux,Neumann
   !use observations,only: SRelaxTau,sProf
!
   IMPLICIT NONE
!
   integer, intent(in) :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer, parameter        :: adv_mode_0=0
   integer, parameter        :: adv_mode_1=1
   REALTYPE                  :: Qsour(0:nlev),Lsour(0:nlev)
   REALTYPE                  :: RelaxTau(0:nlev)
   REALTYPE                  :: dilution
   integer                   :: i
   integer                   :: split,posconc
!
!-----------------------------------------------------------------------
!BOC

   if (.not. rmbm_calc) return

   ! Set local source terms for diffusion scheme to zero,
   ! because source terms are integrated independently later on.
   Qsour    = _ZERO_
   Lsour    = _ZERO_
   
   ! Disable relaxation to observed local values.
   RelaxTau = 1.e15
   
   ! Calculate local pressure
   pres(1:nlev) = -z(1:nlev)
   
#ifdef RMBM_SINGLE_STATE_VARIABLE_ARRAY
   call rmbm_link_state_data        (model,cc(1:ubound(model%info%state_variables,1),1:nlev))
   call rmbm_link_benthos_state_data(model,cc(ubound(model%info%state_variables,1)+1:,1))
#else
   do i=1,ubound(model%info%state_variables,1)
      call rmbm_link_state_data(model,i,cc(i,1:nlev))
   end do
   do i=1,ubound(model%info%state_variables_ben,1)
      call rmbm_link_benthos_state_data(model,i,cc(ubound(model%info%state_variables,1)+i,1))
   end do
#endif

   ! Get updated vertical movement (m/s, positive for upwards) for biological state variables.
#ifdef _RMBM_USE_1D_LOOP_
   call rmbm_get_vertical_movement(model,1,nlev,ws(1:nlev,:))
#else
   do i=1,nlev
      call rmbm_get_vertical_movement(model,i,ws(i,:))
   end do
#endif

   ! Get updated air-sea fluxes for biological state variables.
   call rmbm_get_surface_exchange(model,nlev,sfl)
   
   ! Calculate dilution due to surface freshwater flux (m/s)
   ! If surface freshwater flux is not specified, but surface salinity is relaxed to observations,
   ! calculate the effective dilution from the relation term, and use that instead.
   dilution = precip+evap
   !if (any(SRelaxTau(1:nlev)<1.e10)) &
   !   dilution = dilution + sum((salt(1:nlev)-sProf(1:nlev))/SRelaxTau(1:nlev)*h(2:nlev+1)) &
   !                        /sum(salt(1:nlev)*h(2:nlev+1)) * sum(h(2:nlev+1))

   do i=1,ubound(model%info%state_variables,1)
      ! Add surface flux due to evaporation/precipitation, unless the model explicitly says otherwise.
      if (.not. model%info%state_variables(i)%no_precipitation_dilution) &
         sfl(i) = sfl(i)-cc(i,nlev)*dilution
   
      ! Determine whether the variable is positive-definite based on its lower allowed bound.
      posconc = 0
      if (model%info%state_variables(i)%minimum.ge._ZERO_) posconc = 1
         
      ! Do advection step due to settling or rising
      call adv_center(nlev,dt,h,h,ws(:,i),flux,                   &
           flux,_ZERO_,_ZERO_,w_adv_discr,adv_mode_1,cc(i,:))
         
      ! Do advection step due to vertical velocity
      if (w_adv_method .ne. 0) &
         call adv_center(nlev,dt,h,h,w,flux,                   &
              flux,_ZERO_,_ZERO_,w_adv_ctr,adv_mode_0,cc(i,:))
      
      ! Do diffusion step
      call diff_center(nlev,dt,cnpar,posconc,h,Neumann,Neumann,&
           sfl(i),bfl(i),nuh,Lsour,Qsour,RelaxTau,cc(i,:),cc(i,:))

   end do

   ! Repair state before calling RMBM
   call do_repair_state(nlev,'gotm_rmbm::do_gotm_rmbm, after advection/diffusion')

   do split=1,split_factor
      ! Update local light field (self-shading may have changed through changes in biological state variables)
      call light(nlev,bioshade_feedback)

      ! Time-integrate one biological time step
      call ode_solver(ode_method,ubound(cc,1),nlev,dt_eff,cc(:,0:nlev),right_hand_side_rhs,right_hand_side_ppdd)

      ! Provide RMBM with (pointers to) updated state variables.
#ifdef RMBM_SINGLE_STATE_VARIABLE_ARRAY
      call rmbm_link_state_data        (model,cc(1:ubound(model%info%state_variables,1), 1:nlev))
      call rmbm_link_benthos_state_data(model,cc(ubound(model%info%state_variables,1)+1:,1))
#else
      do i=1,ubound(model%info%state_variables,1)
         call rmbm_link_state_data(model,i,cc(i,1:nlev))
      end do
      do i=1,ubound(model%info%state_variables_ben,1)
         call rmbm_link_benthos_state_data(model,i,cc(ubound(model%info%state_variables,1)+i,1))
      end do
#endif

      ! Repair state
      call do_repair_state(nlev,'gotm_rmbm::do_gotm_rmbm, after time integration')

      ! Time-integrate diagnostic variables defined on horizontonal slices, where needed.
      do i=1,ubound(model%info%diagnostic_variables_hz,1)
         if (model%info%diagnostic_variables_hz(i)%time_treatment.eq.time_treatment_last) then
            ! Simply use last value
            cc_diag_hz(i) = rmbm_get_diagnostic_data_hz(model,i)
         else
            ! Integration or averaging in time needed: for now do simple Forward Euler integration.
            ! If averaging is required, this will be done upon output by diving by the elapsed period.
            cc_diag_hz(i) = cc_diag_hz(i) + rmbm_get_diagnostic_data_hz(model,i)*dt_eff
         end if
      end do
      
      ! Time-integrate diagnostic variables defined on the full domain, where needed.
      do i=1,ubound(model%info%diagnostic_variables,1)
         if (model%info%diagnostic_variables(i)%time_treatment.eq.time_treatment_last) then
            ! Simply use last value
            cc_diag(i,1:nlev) = rmbm_get_diagnostic_data(model,i)
         else
            ! Integration or averaging in time needed: for now do simple Forward Euler integration.
            ! If averaging is required, this will be done upon output by diving by the elapsed period.
            cc_diag(i,1:nlev) = cc_diag(i,1:nlev) + rmbm_get_diagnostic_data(model,i)*dt_eff
         end if
      end do
   end do

   end subroutine do_gotm_rmbm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks the current values of all state variables
!
! !INTERFACE:
   subroutine do_repair_state(nlev,location)
!
! !DESCRIPTION:
! Checks the current values of all state variables and repairs these
! if allowed and possible.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,         intent(in) :: nlev
   character(len=*),intent(in) :: location
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   logical :: valid
#ifndef _RMBM_USE_1D_LOOP_
   integer :: ci
#endif
!
!-----------------------------------------------------------------------
!BOC
#ifdef _RMBM_USE_1D_LOOP_
   call rmbm_check_state(model,1,nlev,repair_state,valid)
#else
   do ci=1,nlev
      call rmbm_check_state(model,ci,repair_state,valid)
      if (.not.(valid.or.repair_state)) exit
   end do
#endif   
   if (.not. (valid .or. repair_state)) then
      FATAL 'State variable values are invalid and repair is not allowed.'
      FATAL location
      stop 'gotm_rmbm::do_repair_state'
   end if

   end subroutine do_repair_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculates temporal derivatives as a derivative vector
!
! !INTERFACE:
   subroutine right_hand_side_rhs(first,numc,nlev,cc,rhs)
!
! !DESCRIPTION:
! Checks the current values of all state variables and repairs these
! if allowed and possible.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc,nlev
   REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: rhs(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: i,n
!
!-----------------------------------------------------------------------
!BOC
   ! Shortcut to the number of pelagic state variables.
   n = ubound(model%info%state_variables,1)

   ! Provide RMBM with (pointers to) the current state.
#ifdef RMBM_SINGLE_STATE_VARIABLE_ARRAY
   call rmbm_link_state_data        (model,cc(1:n, 1:nlev))
   call rmbm_link_benthos_state_data(model,cc(n+1:,1))
#else
   do i=1,ubound(model%info%state_variables,1)
      call rmbm_link_state_data(model,i,cc(i,1:nlev))
   end do
   do i=1,ubound(model%info%state_variables_ben,1)
      call rmbm_link_benthos_state_data(model,i,cc(n+i,1))
   end do
#endif

   ! If this is not the first step in the (multi-step) integration scheme,
   ! then first make sure that the intermediate state variable values are valid.
   if (.not. first) call do_repair_state(nlev,'gotm_rmbm::right_hand_side_rhs')

   ! Initialization is needed because the different biogeochemical models increment or decrement
   ! the temporal derivatives, rather than setting them directly. This is needed for the simultaenous
   ! running of different coupled BGC models.
   rhs = _ZERO_

   ! Calculate temporal derivatives due to benthic processes.
   call rmbm_do_benthos(model,1,rhs(1:n,1),rhs(n+1:,1))
   
   ! Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
   rhs(1:n,1) = rhs(1:n,1)/h(2)

   ! Add pelagic sink and source terms for all depth levels.
#ifdef _RMBM_USE_1D_LOOP_
   call rmbm_do(model,1,nlev,rhs(1:n,1:nlev))
#else   
   do i=1,nlev
      call rmbm_do(model,i,rhs(1:n,i))
   end do
#endif

   end subroutine right_hand_side_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculates temporal derivatives as production/destruction matrices
!
! !INTERFACE:
   subroutine right_hand_side_ppdd(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! Checks the current values of all state variables and repairs these
! if allowed and possible.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,  intent(in)                 :: first
   integer,  intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: i,n
!
!-----------------------------------------------------------------------
!BOC
   ! Shortcut to the number of pelagic state variables.
   n = ubound(model%info%state_variables,1)

   ! Provide RMBM with (pointers to) the current state.
#ifdef RMBM_SINGLE_STATE_VARIABLE_ARRAY
   call rmbm_link_state_data        (model,cc(1:n,1:nlev))
   call rmbm_link_benthos_state_data(model,cc(n+1:,1))
#else
   do i=1,ubound(model%info%state_variables,1)
      call rmbm_link_state_data(model,i,cc(i,1:nlev))
   end do
   do i=1,ubound(model%info%state_variables_ben,1)
      call rmbm_link_benthos_state_data(model,i,cc(n+i,1))
   end do
#endif

   ! If this is not the first step in the (multi-step) integration scheme,
   ! then first make sure that the intermediate state variable values are valid.
   if (.not. first) call do_repair_state(nlev,'gotm_rmbm::right_hand_side_ppdd')

   ! Initialiaze production and destruction matrices to zero because RMBM
   ! biogeochemical models increment these, rather than set these.
   pp = _ZERO_
   dd = _ZERO_
   
   ! Calculate temporal derivatives due to benthic processes.
   call rmbm_do_benthos(model,1,pp(:,:,1),dd(:,:,1),n)
   
   ! Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
   pp(1:n,:,1) = pp(1:n,:,1)/h(2)
   dd(1:n,:,1) = dd(1:n,:,1)/h(2)

   ! Add pelagic sink and source terms for all depth levels.
#ifdef _RMBM_USE_1D_LOOP_
   call rmbm_do(model,1,nlev,pp(1:n,1:n,1:nlev),dd(1:n,1:n,1:nlev))
#else
   do i=1,nlev
      call rmbm_do(model,i,pp(1:n,1:n,i),dd(1:n,1:n,i))
   end do
#endif

   end subroutine right_hand_side_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish biogeochemical model
!
! !INTERFACE:
   subroutine clean_gotm_rmbm
!
! !DESCRIPTION:
!  Deallocate memory.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'clean_gotm_rmbm'

   ! Deallocate internal arrays
   if (allocated(cc))             deallocate(cc)
   if (allocated(cc_diag))        deallocate(cc_diag)
   if (allocated(cc_diag_hz))     deallocate(cc_diag_hz)
   if (allocated(ws))             deallocate(ws)
   if (allocated(sfl))            deallocate(sfl)
   if (allocated(bfl))            deallocate(bfl)
   if (allocated(total))          deallocate(total)
   if (allocated(local))          deallocate(local)
   if (allocated(par))            deallocate(par)
   if (allocated(pres))           deallocate(pres)
   LEVEL1 'done.'

   end subroutine clean_gotm_rmbm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate light over entire column
!
! !INTERFACE:
   subroutine light(nlev,bioshade_feedback)
!
! !DESCRIPTION:
! Calculate photosynthetically active radiation over entire column
! based on surface radiation, and background and biotic extinction.
!
! !USES:
   use observations, only:A,g2
   
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   logical, intent(in)                 :: bioshade_feedback
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: i
   REALTYPE :: zz,bioext,localext
#ifdef _RMBM_USE_1D_LOOP_
   REALTYPE :: localexts(1:nlev)
#endif
!
!-----------------------------------------------------------------------
!BOC
   zz = _ZERO_
   bioext = _ZERO_

#ifdef _RMBM_USE_1D_LOOP_
   call rmbm_get_light_extinction(model,1,nlev,localexts)
#endif
   do i=nlev,1,-1
#ifdef _RMBM_USE_1D_LOOP_
      localext = localexts(i)
#else
      call rmbm_get_light_extinction(model,i,localext)
#endif
   
      ! Add the extinction of the first half of the grid box.
      bioext = bioext+localext*0.5*h(i+1)

      zz=zz+0.5*h(i+1)
      par(i)=rad(nlev)*(_ONE_-A)*exp(-zz/g2-bioext)

      ! Add the extinction of the second half of the grid box.
      bioext = bioext+localext*0.5*h(i+1)
      
      zz=zz+0.5*h(i+1)
      if (bioshade_feedback) bioshade(i)=exp(-bioext)
   end do

   end subroutine light
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize output
!
! !INTERFACE:
   subroutine init_output_gotm_rmbm()
   
!
! !DESCRIPTION:
!  Initialize the output by defining biogeochemical variables.
!
! !USES:
   use output,  only: out_fmt
#ifdef NETCDF_FMT
   use ncdfout, only: ncid,lon_dim,lat_dim,z_dim,time_dim,dims
   use ncdfout, only: define_mode,new_nc_variable,set_attributes
#endif
!
   IMPLICIT NONE
!
#ifdef NETCDF_FMT
#include "netcdf.inc"
#endif
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: iret,n
!
!-----------------------------------------------------------------------
!BOC
   if (.not. rmbm_calc) return
   
   select case (out_fmt)
      case (NETCDF)
#ifdef NETCDF_FMT
         ! Put NetCDF library in define mode.
         iret = define_mode(ncid,.true.)

         ! Set up dimension indices for 4D variables (longitude,latitude,depth,time).
         dims(1) = lon_dim
         dims(2) = lat_dim
         dims(3) = z_dim
         dims(4) = time_dim

         ! Add a NetCDF variable for each 4D (longitude,latitude,depth,time) biogeochemical state variable.
         do n=1,ubound(model%info%state_variables,1)
            iret = new_nc_variable(ncid,model%info%state_variables(n)%name,NF_REAL, &
                                   4,dims,model%info%state_variables(n)%externalid)
            iret = set_attributes(ncid,model%info%state_variables(n)%externalid,       &
                                  units=model%info%state_variables(n)%units,    &
                                  long_name=model%info%state_variables(n)%longname)
         end do

         ! Add a NetCDF variable for each 4D (longitude,latitude,depth,time) biogeochemical diagnostic variable.
         do n=1,ubound(model%info%diagnostic_variables,1)
            iret = new_nc_variable(ncid,model%info%diagnostic_variables(n)%name,NF_REAL, &
                                   4,dims,model%info%diagnostic_variables(n)%externalid)
            iret = set_attributes(ncid,model%info%diagnostic_variables(n)%externalid,    &
                                  units=model%info%diagnostic_variables(n)%units,        &
                                  long_name=model%info%diagnostic_variables(n)%longname)
         end do

         ! Set up dimension indices for 3D variables (longitude,latitude,time).
         dims(3) = time_dim

         ! Add a NetCDF variable for each 3D (longitude,latitude,time) biogeochemical state variable.
         do n=1,ubound(model%info%state_variables_ben,1)
            iret = new_nc_variable(ncid,model%info%state_variables_ben(n)%name,NF_REAL, &
                                   3,dims,model%info%state_variables_ben(n)%externalid)
            iret = set_attributes(ncid,model%info%state_variables_ben(n)%externalid,    &
                                  units=model%info%state_variables_ben(n)%units,        &
                                  long_name=model%info%state_variables_ben(n)%longname)
         end do

         ! Add a NetCDF variable for each 3D (longitude,latitude,time) biogeochemical diagnostic variable.
         do n=1,ubound(model%info%diagnostic_variables_hz,1)
            iret = new_nc_variable(ncid,model%info%diagnostic_variables_hz(n)%name,NF_REAL, &
                                   3,dims,model%info%diagnostic_variables_hz(n)%externalid)
            iret = set_attributes(ncid,model%info%diagnostic_variables_hz(n)%externalid,    &
                                  units=model%info%diagnostic_variables_hz(n)%units,        &
                                  long_name=model%info%diagnostic_variables_hz(n)%longname)
         end do

         ! Add a variable for each conserved quantity
         do n=1,ubound(model%info%conserved_quantities,1)
            iret = new_nc_variable(ncid,trim(model%info%conserved_quantities(n)%name)//'_tot',NF_REAL, &
                                   3,dims,model%info%conserved_quantities(n)%externalid)
            iret = set_attributes(ncid,model%info%conserved_quantities(n)%externalid,      &
                                  units='m*'//model%info%conserved_quantities(n)%units,    &
                                  long_name=trim(model%info%conserved_quantities(n)%longname)//', depth-integrated')
         end do

         ! Take NetCDF library out of define mode (ready for storing data).
         iret = define_mode(ncid,.false.)
#endif
   end select

   end subroutine init_output_gotm_rmbm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Save values of biogeochemical variables
!
! !INTERFACE:
   subroutine save_gotm_rmbm(nlev)
   
!
! !DESCRIPTION:
!  Save properties of biogeochemical model, including state variable
!  values, diagnostic variable values, and sums of conserved quantities.
!
! !USES:
   use output,  only: nsave,out_fmt
#ifdef NETCDF_FMT
   use ncdfout, only: ncid
   use ncdfout, only: store_data
#endif

   IMPLICIT NONE

#ifdef NETCDF_FMT
#include "netcdf.inc"
#endif
!
! !INPUT PARAMETERS:
   integer, intent(in)                  :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer :: iret,n
!
!-----------------------------------------------------------------------
!BOC
   if (.not. rmbm_calc) return
   
   select case (out_fmt)
      case (NETCDF)
#ifdef NETCDF_FMT
         ! Store pelagic biogeochemical state variables.
         do n=1,ubound(model%info%state_variables,1)
            iret = store_data(ncid,model%info%state_variables(n)%externalid,XYZT_SHAPE,nlev,array=cc(n,0:nlev))
         end do

         ! Store benthic biogeochemical state variables.
         do n=1,ubound(model%info%state_variables_ben,1)
            iret = store_data(ncid,model%info%state_variables_ben(n)%externalid,XYT_SHAPE,1,scalar=cc(ubound(model%info%state_variables,1)+n,1))
         end do

         ! Process and store diagnostic variables defined on the full domain.
         do n=1,ubound(model%info%diagnostic_variables,1)
            ! Time-average diagnostic variable if needed.
            if (model%info%diagnostic_variables(n)%time_treatment==time_treatment_averaged) &
               cc_diag(n,1:nlev) = cc_diag(n,1:nlev)/(nsave*dt)
               
            ! Store diagnostic variable values.
            iret = store_data(ncid,model%info%diagnostic_variables(n)%externalid,XYZT_SHAPE,nlev,array=cc_diag(n,0:nlev))
            
            ! Reset diagnostic variables to zero if they will be time-integrated (or time-averaged).
            if (model%info%diagnostic_variables(n)%time_treatment==time_treatment_averaged .or. &
                model%info%diagnostic_variables(n)%time_treatment==time_treatment_step_integrated) &
               cc_diag(n,1:nlev) = _ZERO_
         end do

         ! Process and store diagnostic variables defined on horizontal slices of the domain.
         do n=1,ubound(model%info%diagnostic_variables_hz,1)
            ! Time-average diagnostic variable if needed.
            if (model%info%diagnostic_variables_hz(n)%time_treatment==time_treatment_averaged) &
               cc_diag_hz(n) = cc_diag_hz(n)/(nsave*dt)
               
            ! Store diagnostic variable values.
            iret = store_data(ncid,model%info%diagnostic_variables_hz(n)%externalid,XYT_SHAPE,nlev,scalar=cc_diag_hz(n))
            
            ! Reset diagnostic variables to zero if they will be time-integrated (or time-averaged).
            if (model%info%diagnostic_variables_hz(n)%time_treatment==time_treatment_averaged .or. &
                model%info%diagnostic_variables_hz(n)%time_treatment==time_treatment_step_integrated) &
               cc_diag_hz(n) = _ZERO_
         end do

         ! Integrate conserved quantities over depth.
#ifdef _RMBM_USE_1D_LOOP_
         call rmbm_get_conserved_quantities(model,1,nlev,local)
         do n=1,ubound(model%info%conserved_quantities,1)
            ! Note: our pointer to h has a lower bound of 1, while the original pointed-to data starts at 0.
            ! We therefore need to increment the index by 1 in order to address original elements >=1!
            total(n) = sum(h(2:nlev+1)*local(1:nlev,n))
         end do
#else
         total = _ZERO_
         do n=1,nlev
            ! Note: our pointer to h has a lower bound of 1, while the original pointed-to data starts at 0.
            ! We therefore need to increment the index by 1 in order to address original elements >=1!
            call rmbm_get_conserved_quantities(model,n,local)
            total = total + h(n+1)*local
         end do
#endif

         ! Store conserved quantity integrals.
         do n=1,ubound(model%info%conserved_quantities,1)
            iret = store_data(ncid,model%info%conserved_quantities(n)%externalid,XYT_SHAPE,1,scalar=total(n))
         end do
#endif
   end select

   end subroutine save_gotm_rmbm
!EOC

   end module gotm_rmbm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------

