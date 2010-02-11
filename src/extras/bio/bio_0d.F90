!$Id: bio_0d.F90,v 1.10 2010-01-20 16:57:51 jorn Exp $
#include"cppdefs.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_0d --- Interface between 1D GOTM and 0D bio models
!
! !INTERFACE:
   module bio_0d
!
! !DESCRIPTION:
! This module is the interface between GOTM and the library of 0d
! biogeochemical models. Note that this model library is accessed via
! the bio_0d_gen module: new 0d models should be registered there,
! and not in the present bio_0d module.
!
! !USES:
   use bio_var
   use bio_types
   use bio_0d_gen
   use time,    only: timestep

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public model,init_bio_0d, init_var_0d,                 &
          light_0d, light_0d_par, &
          surface_fluxes_0d, update_sinking_rates_0d, &
          do_bio_0d_eul, do_bio_0d_eul_rhs, do_bio_0d_par, &
          save_bio_0d, end_bio_0d
!
! !PRIVATE DATA MEMBERS:
   type (type_model_collection) :: model
   
   ! Lagrangian model
   integer :: np,nt
   REALTYPE, allocatable :: shade(:),extinction(:),cc_loc(:,:)
   
   REALTYPE, allocatable :: cc_diag(:,:)

!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
   
   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the 0d biological framework
!
! !INTERFACE:
   subroutine init_bio_0d(namlst,unit)
!
! !DESCRIPTION:
!  TODO
!
! !USES:
   use observations, only: A,g2
!
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: namlst
   integer,          intent(in)   :: unit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer :: i,ioffset
   character(len=64) :: model_name
   integer,parameter :: bio_model2=1002
   character(len=64) :: bio_model2_nml = 'bio_jellyfish.nml'
   !integer,parameter :: bio_model2=1003
   !character(len=64) :: bio_model2_nml = 'bio_co2_sys.nml'

!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_0d'
   
   ! Get actual model info based on the model selected.
   model_name = get_model_name(bio_model)
   
   ! Use single instance
   !model = init_bio_0d_generic(1,(/bio_model/),namlst,(/'bio_'//trim(model_name)//'.nml'/))
   
   ! Use two instances
   model = init_bio_0d_generic(2,(/bio_model,bio_model2/),namlst,(/'bio_'//trim(model_name)//'.nml',trim(bio_model2_nml)/))
   !model = init_bio_0d_generic(2,(/bio_model2,bio_model/),namlst,(/trim(bio_model2_nml),'bio_'//trim(model_name)//'.nml               '/))
   
   ! Allocate global arrays for info on biogeochemical model
   ! Add a variable for particle densities if using Lagragian model
   numc = model%info%state_variable_count
   if (.not. bio_eulerian) numc = numc+1
   call bio_alloc_info

   ! If using Lagrangian model, the first variable will describe the number
   ! particles per unit of volume.
   if (.not. bio_eulerian) then
      var_names(1) = 'Np'
      var_units(1) = 'counts/volume'
      var_long (1) = 'number of particles per volume'
      ioffset = 1
      
      ntype = 1
      nprop = model%info%state_variable_count
   else
      ioffset = 0
   end if

   ! Register the variables used by the biogeochemical model in the global arrays.
   do i=1,model%info%state_variable_count
      var_names(i+ioffset) = model%info%variables(i)%name
      var_units(i+ioffset) = model%info%variables(i)%unit
      var_long (i+ioffset) = model%info%variables(i)%longname
   end do
   
   write (*,*) 'The biogeochemical model has a total of:'
   write (*,*) model%info%state_variable_count,     ' state variables'
   write (*,*) model%info%diagnostic_variable_count,' diagnostic variables'
   write (*,*) model%info%conserved_quantity_count, ' conserved quantities'

   end subroutine init_bio_0d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fill initial particle distribution
!
! !INTERFACE:
   subroutine init_par_bio_0d()
!
! !DESCRIPTION:
! Particles are distributed homogeneously over the whole water column. 
! Indices of vertical grid cells are assigend to all particles, and the
! particles are marked active.  
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf, Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n

!-----------------------------------------------------------------------
!BOC


!  create homogeneous particle distribution
   do j=1,ntype
      do n=1,npar
         par_z(n,j) = zbot + n/float(npar+1)*(ztop-zbot)
      end do
   end do

!  assign cell indices to particles
   do j=1,ntype
      do n=1,npar
         do i=1,nlev
            if (zlev(i) .gt. par_z(n,j)) exit
         end do
         par_ind(n,j) = i
         par_act(n,j)  = .true.
      end do
   end do


   return
   end subroutine init_par_bio_0d
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the biological variables
!
! !INTERFACE:
   subroutine init_var_0d
!
! !DESCRIPTION:
!  TODO
!
! !USES:
   IMPLICIT NONE

! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman

! !LOCAL VARIABLES:
   integer :: j,ioffset

!EOP
!-----------------------------------------------------------------------
!BOC

   ioffset = 0
   if (.not. bio_eulerian) ioffset = 1
   
   do j=1,model%info%state_variable_count
      ws(j+ioffset,0:nlev) = model%info%variables(j)%vertical_movement
      posconc(j+ioffset) = model%info%variables(j)%positive_definite
#if 0
      mussels_inhale(j+ioffset) = model%info%variables(j)%mussels_inhale
#endif
   end do

   if (.not. bio_eulerian) then
      nt = 1

      do j=1,model%info%state_variable_count
         par_prop(:,j,nt) = (ztop-zbot)*model%info%variables(j)%initial_value/npar
      end do
      
      ! Configure particle variable
      cc = _ZERO_
      ws(1,0:nlev) = _ZERO_
      posconc(1) = 1
#if 0
      mussels_inhale(1) = .false.
#endif

      call init_par_bio_0d

      allocate(cc_loc(1:model%info%state_variable_count,0:1))
   else
      do j=1,numc
         cc(j,1:nlev) = model%info%variables(j)%initial_value
      end do
   end if
   
   ! Array for diagnostic variables
   ! NB it needs to have lower bound 0 for the second dimension (depth),
   ! even though the values at index 0 will never be used.
   ! Reason: store_data expects a lower bound of zero.
   allocate(cc_diag(1:model%info%diagnostic_variable_count,0:nlev))
   
   ! Initialize diagnostic variable values to zero, because for time-integrated
   ! and time-averaged variables, the array will be incremented rather than set directly.
   cc_diag = _ZERO_

   ! Arrays for particle model
   allocate(shade(0:nlev))
   allocate(extinction(1:nlev))

   end subroutine init_var_0d
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes
!
! !INTERFACE:
   subroutine surface_fluxes_0d(nlev)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer                              :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   call update_airsea_exchange_bio_0d_generic(model,nlev,numc,sfl)

   end subroutine surface_fluxes_0d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Sinking rates
!
! !INTERFACE:
   subroutine update_sinking_rates_0d(numc,nlev,cc)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer                              :: numc,nlev
  REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                    :: ci
!EOP
!-----------------------------------------------------------------------
!BOC

   ! If sinking rates are constant (time- and space-independent, just return.
   if (model%info%dynamic_vertical_movement.eq.0) return

   ! Iterate over all depth levels
   do ci=1,nlev
      call get_vertical_movement_bio_0d_generic(model,ci,ws(:,ci))
   end do

   end subroutine update_sinking_rates_0d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for the 0d biological framework
!
! !INTERFACE:
   subroutine light_0d(nlev,bioshade_feedback)
!
! !DESCRIPTION:
! TODO
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
! !LOCAL VARIABLES:
   integer :: i
   REALTYPE :: zz,bioext

!EOP
!-----------------------------------------------------------------------
!BOC
   zz = _ZERO_
   bioext = _ZERO_
   do i=nlev,1,-1
      ! Add the extinction of the first half of the grid box.
      bioext = bioext+get_bio_extinction_bio_0d_generic(model,i)*0.5*h(i)

      zz=zz+0.5*h(i)
      par(i)=rad(nlev)*(_ONE_-A)*exp(-zz/g2-bioext)

      ! Add the extinction of the second half of the grid box.
      bioext = bioext+get_bio_extinction_bio_0d_generic(model,i)*0.5*h(i)
      
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade_(i)=exp(-bioext)
   end do

   end subroutine light_0d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for particle version of the 0d biological framework
!
! !INTERFACE:
   subroutine light_0d_par(model,nlev,bioshade_feedback)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),intent(in)        :: model
   integer, intent(in)                 :: nlev
   logical, intent(in)                 :: bioshade_feedback
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer :: ind,i
   REALTYPE :: bioext

!EOP
!-----------------------------------------------------------------------
!BOC
   nt = 1

   ! Get bio-related extinction in each grid box
   extinction = _ZERO_
   do np=1,npar
      if (par_act(np,nt)) then

         ! Get depth index of the current particle
         ind = par_ind(np,nt)

         ! Get the light extinction coefficient combined over all state variables.
         !extinction(ind) = extinction(ind) + get_bio_extinction_bio_0d_generic(model,par_prop(np,:,nt)/h(ind))
      end if
   end do

   ! Calculate shade factor at each depth
   bioext = _ZERO_
   shade(nlev) = _ZERO_
   do i=nlev,1,-1
      bioext = bioext+extinction(i)*h(i)
      shade(i-1) = exp(-bioext) ! Note: this is the shade factor at the top of the grid box
   end do
   
   ! Set feedback of bioturbidity to physics if needed
   if (bioshade_feedback) bioshade_(1:nlev)=shade(0:nlev-1)
   end subroutine light_0d_par
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of eulerian 0D biogeochemical model
!
! !INTERFACE:
   subroutine do_bio_0d_eul(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! TODO
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
   REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                    :: ci
   REALTYPE                   :: diag(1:model%info%diagnostic_variable_count,nlev)
!EOP
!-----------------------------------------------------------------------
!BOC

!KBK - is it necessary to initialise every time - expensive in a 3D model
   pp = _ZERO_
   dd = _ZERO_
   
   ! Iterate over all depth levels
   do ci=1,nlev
      call do_bio_0d_generic(model,ci,numc,model%info%diagnostic_variable_count, &
               pp(:,:,ci),dd(:,:,ci),diag(:,ci))
   end do
   
   if (first) then
      ! First time during this time step that do_bio is called: store diagnostic values.
      ! NB. higher order integration schemes may call this routine multiple times.
      ! In that case only the value at the first call is used (essential because time
      ! integration is done internally).
      
      do ci=1,model%info%diagnostic_variable_count
         if (model%info%diagnostic_variables(ci)%time_treatment.eq.0) then
            ! Simply use last value
            cc_diag(ci,1:nlev) = diag(ci,1:nlev)
         else
            ! Integration or averaging in time needed: for now do simple Forward Euler integration.
            cc_diag(ci,1:nlev) = cc_diag(ci,1:nlev) + diag(ci,1:nlev)*timestep
         end if
      end do
   end if
   
   end subroutine do_bio_0d_eul
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of eulerian 0D biogeochemical model
!
! !INTERFACE:
   subroutine do_bio_0d_eul_rhs(first,numc,nlev,cc,rhs)
!
! !DESCRIPTION:
! TODO
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
! !LOCAL VARIABLES:
   integer                    :: ci
   REALTYPE                   :: diag(1:model%info%diagnostic_variable_count,nlev)
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Initialization is needed because the different biogeochemical models increment or decrement
   ! the temporal derivatives, rather than setting them directly. Tis is needed for the simultaenous
   ! running of different coupled BGC models.
   rhs = _ZERO_

   ! Iterate over all depth levels
   do ci=1,nlev
      call do_bio_0d_generic(model,ci,numc,model%info%diagnostic_variable_count,rhs(:,ci),diag(:,ci))
   end do
   
   if (first) then
      ! First time during this time step that do_bio is called: store diagnostic values.
      ! NB. higher order integration schemes may call this routine multiple times.
      ! In that case only the value at the first call is used (essential because time
      ! integration is done internally).
      
      do ci=1,model%info%diagnostic_variable_count
         if (model%info%diagnostic_variables(ci)%time_treatment.eq.0) then
            ! Simply use last value
            cc_diag(ci,1:nlev) = diag(ci,1:nlev)
         else
            ! Integration or averaging in time needed: for now do simple Forward Euler integration.
            cc_diag(ci,1:nlev) = cc_diag(ci,1:nlev) + diag(ci,1:nlev)*timestep
         end if
      end do
   end if
   
   end subroutine do_bio_0d_eul_rhs
!EOC

!BOP
!
! !IROUTINE: Update particle model described by 0D biogeochemical model.
!
! !INTERFACE:
   subroutine do_bio_0d_par(ode_method,dt)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in) :: ode_method
   REALTYPE,intent(in) :: dt
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                    :: i,j,ind
!EOP
!-----------------------------------------------------------------------
!BOC

   nt = 1
   
   ! First column of state variable array is not used by ode solvers.
   cc_loc(1:model%info%state_variable_count,0) = _ZERO_

   ! For each particle: integrate over the time step
   do np=1,npar
      ! Get the current state
      cc_loc(1:model%info%state_variable_count,1) = par_prop(np,1:model%info%state_variable_count,nt)
      
      ! Call ode solver to get updated state
      !call ode_solver(ode_method,model%info%state_variable_count,1,dt,cc_loc,get_bio_0d_par_rhs)
      
      ! Store updated state as particle properties
      par_prop(np,1:model%info%state_variable_count,nt) = cc_loc(1:model%info%state_variable_count,1)
   end do
   
   ! Initialize array with [Eulerian] summary statistics
   cc = _ZERO_

   do np=1,npar
      if (par_act(np,nt)) then

         ! Get depth index of the current particle
         ind = par_ind(np,nt) 

         ! Count particles per grid volume
         cc(1,ind) = cc(1,ind) + _ONE_
         do j=1,model%info%state_variable_count
            cc(j+1,ind) = cc(j+1,ind) + par_prop(np,j,nt)
         end do

      end if

   end do

   ! Compute volume averages
   do i=1,nlev
      cc(:,i) = cc(:,i)/h(i)
   end do
   
   end subroutine do_bio_0d_par
!EOC

!BOP
!
! !IROUTINE: Get the right-hand side of the ODE system for the current particle
!
! !INTERFACE:
   subroutine get_bio_0d_par_rhs(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   use observations,only:A,g2

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc,nlev
   REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                   :: rat, localshade, diag(model%info%diagnostic_variable_count)
   integer                    :: i
!EOP
!-----------------------------------------------------------------------
!BOC

   i=par_ind(np,nt)
   rat=(par_z(np,nt)-zlev(i-1))/h(i)

   localshade = rat*shade(i)+(1.-rat)*shade(i-1)
   !env_par%par  = rad(nlev)*(_ONE_-A)*exp(par_z(np,nt)/g2)*localshade

   ! Linearly interpolate environmental conditions
   !env_par%z    = par_z(np,nt)
   !env_par%t    = rat*  t(i)+(1.-rat)*  t(i-1)
   !env_par%s    = rat*  s(i)+(1.-rat)*  s(i-1)
   !env_par%nuh  = rat*nuh(i)+(1.-rat)*nuh(i-1)
   !env_par%rho  = rat*rho(i)+(1.-rat)*rho(i-1)
   !call do_bio_0d_generic(model,first,par_prop(np,1:model%info%state_variable_count,nt),env_par,pp(:,:,1),dd(:,:,1),diag)
   
   end subroutine get_bio_0d_par_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine save_bio_0d(first,out_fmt,out_unit,ncid)
!
! !DESCRIPTION:
!  Save additional properties of 0d biogeochemical model
!
! !USES:
   use output,  only: nsave
   use ncdfout, only: lon_dim,lat_dim,z_dim,time_dim,dims
   use ncdfout, only: define_mode,new_nc_variable,set_attributes,store_data

   IMPLICIT NONE

#ifdef NETCDF_FMT
#include "netcdf.inc"
#endif
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: out_fmt,out_unit,ncid
!
! !LOCAL VARIABLES:
   integer :: n,iret,ilev
   REALTYPE :: total(1:model%info%conserved_quantity_count),local(1:model%info%conserved_quantity_count)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (out_fmt)
      case (NETCDF)
#ifdef NETCDF_FMT
         if(first) then
            iret = define_mode(ncid,.true.)

            dims(1) = lon_dim
            dims(2) = lat_dim
            dims(3) = z_dim
            dims(4) = time_dim

            ! Add a variable for each diagnostic variable
            do n=1,model%info%diagnostic_variable_count
               iret = new_nc_variable(ncid,model%info%diagnostic_variables(n)%name,NF_REAL, &
                                      4,dims,model%info%diagnostic_variables(n)%id)
               iret = set_attributes(ncid,model%info%diagnostic_variables(n)%id,       &
                                     units=model%info%diagnostic_variables(n)%unit,    &
                                     long_name=model%info%diagnostic_variables(n)%longname)
            end do

            dims(3) = time_dim

            ! Add a variable for each conserved quantity
            do n=1,model%info%conserved_quantity_count
               iret = new_nc_variable(ncid,'tot_'//model%info%conserved_quantities(n)%name,NF_REAL, &
                                      3,dims,model%info%conserved_quantities(n)%id)
               iret = set_attributes(ncid,model%info%conserved_quantities(n)%id,       &
                                     units='m*'//model%info%conserved_quantities(n)%unit,    &
                                     long_name='depth-integrated '//model%info%conserved_quantities(n)%longname)
            end do

            iret = define_mode(ncid,.false.)
         end if

         do n=1,model%info%diagnostic_variable_count
            if (model%info%diagnostic_variables(n)%time_treatment==2) &
               cc_diag(n,1:nlev) = cc_diag(n,1:nlev)/(nsave*timestep)
            iret = store_data(ncid,model%info%diagnostic_variables(n)%id,XYZT_SHAPE,nlev,array=cc_diag(n,0:nlev))
            if (model%info%diagnostic_variables(n)%time_treatment==2 .or. &
                model%info%diagnostic_variables(n)%time_treatment==3) &
               cc_diag(n,1:nlev) = _ZERO_
         end do

         ! Integrate conserved quantities over depth
         total = _ZERO_
         do ilev=1,nlev
            call get_conserved_quantities_bio_0d_generic(model,ilev,local)
            total = total + h(ilev)*local
         end do

         do n=1,model%info%conserved_quantity_count
            iret = store_data(ncid,model%info%conserved_quantities(n)%id,XYT_SHAPE,1,scalar=total(n))
         end do
#endif
   end select

   end subroutine save_bio_0d
!EOC
   
   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_0d
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
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

   end subroutine end_bio_0d
!EOC

!-----------------------------------------------------------------------

   end module bio_0d

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
