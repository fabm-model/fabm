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
   use util,only: flux,Neumann

   use rmbm
   use rmbm_types
   
   implicit none

   interface   
      subroutine ode_solver(solver,numc,nlev,dt,cc,right_hand_side_rhs,right_hand_side_ppdd)
         integer, intent(in)                 :: solver,nlev,numc
         REALTYPE, intent(in)                :: dt
         REALTYPE, intent(inout)             :: cc(1:numc,0:nlev)

         interface
            subroutine right_hand_side_ppdd(first,numc,nlev,cc,pp,dd)
               logical, intent(in)                  :: first
               integer, intent(in)                  :: numc,nlev
               REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
               
               REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
               REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
            end
         end interface

         interface
            subroutine right_hand_side_rhs(first,numc,nlev,cc,rhs)
               logical, intent(in)                  :: first
               integer, intent(in)                  :: numc,nlev
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
   integer                   :: w_adv_discr,ode_method,split_factor
   logical                   :: rmbm_calc,bioshade_feedback

   ! Model
   type (type_model),pointer :: model
   
   ! Arrays for state and diagnostic variables
   REALTYPE,allocatable,dimension(LOCATION_DIMENSIONS,:),target :: cc
   REALTYPE,allocatable,dimension(LOCATION_DIMENSIONS,:)        :: cc_diag

   ! Arrays for work, vertical movement, and cross-boundary fluxes
   REALTYPE,allocatable,dimension(LOCATION_DIMENSIONS,:) :: work_cc_diag,ws
   REALTYPE,allocatable,dimension(:)                     :: sfl,bfl,total,local
   
   ! Arrays for environmental variables not supplied externally.
   REALTYPE,allocatable,dimension(LOCATION_DIMENSIONS)   :: par,pres
   
   ! External variables
   REALTYPE :: dt          ! External time step
   integer  :: w_adv_ctr   ! Scheme for vertical advection (0 if not used)
   REALTYPE,pointer,dimension(LOCATION_DIMENSIONS) :: nuh,h,bioshade,rad,w,z

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_gotm_rmbm(namlst,fname,unit,nmax_)
!
! !DESCRIPTION: 
! TODO
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: nmax_

!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!  local variables
   integer                   :: i
   integer                   :: models(256)
   type (type_model),pointer :: childmodel
   namelist /bio_nml/ rmbm_calc,models,                                 &
                      cnpar,w_adv_discr,ode_method,split_factor,        &
                      bioshade_feedback
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_gotm_rmbm'
   
   ! Initialize RMBM model identifiers to invalid id.
   rmbm_calc         = .false.
   models            = -1
   cnpar             = _ONE_
   w_adv_discr       = 6
   ode_method        = 1
   split_factor      = 1
   bioshade_feedback = .true.

   ! Open the namelist file and read the namelist
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_nml,err=99)

   if (rmbm_calc) then
      ! Create model tree
      model => rmbm_create_model()
      do i=1,ubound(models,1)
         if (models(i).ne.-1) &
            childmodel => rmbm_create_model(models(i),parent=model)
      end do
      
      ! Initialize model tree (creates metadata and assigns variable identifiers)
      call rmbm_init(model,namlst)

      ! Report prognostic variable descriptions
      LEVEL2 'RMBM state variables:'
      do i=1,model%info%state_variable_count
         LEVEL3 trim(model%info%variables(i)%name), '  ', &
                trim(model%info%variables(i)%units),'  ',&
                trim(model%info%variables(i)%longname)
      end do

      ! Report diagnostic variable descriptions
      LEVEL2 'RMBM diagnostic variables:'
      do i=1,model%info%diagnostic_variable_count
         LEVEL3 trim(model%info%diagnostic_variables(i)%name), '  ', &
                trim(model%info%diagnostic_variables(i)%units),'  ',&
                trim(model%info%diagnostic_variables(i)%longname)
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

98 LEVEL2 'I could not open rmbm.nml'
   LEVEL2 'If thats not what you want you have to supply rmbm.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working rmbm.nml'
   rmbm_calc = .false.
   return
99 FATAL 'I could not read rmbm.nml'
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
   
   LOCATION_TYPE,intent(in) :: LOCATION
!
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!  local variables
   integer                   :: i,rc,varid

!-----------------------------------------------------------------------
!BOC
   if (.not. rmbm_calc) return

   allocate(cc(1:model%info%state_variable_count,LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (cc)'
   do i=1,model%info%state_variable_count
      cc(i,:) = model%info%variables(i)%initial_value
      model%state(i)%data => cc(i,1:LOCATION)
   end do

   allocate(cc_diag(1:model%info%diagnostic_variable_count,LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag)'
   cc_diag = _ZERO_

   allocate(work_cc_diag(1:model%info%diagnostic_variable_count,LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (work_cc_diag)'
   work_cc_diag = _ZERO_

   allocate(ws(LOCATION_RANGE,1:model%info%state_variable_count),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   do i=1,model%info%state_variable_count
      ws(:,i) = model%info%variables(i)%vertical_movement
   end do

   allocate(sfl(1:model%info%state_variable_count),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (sfl)'
   sfl = _ZERO_

   allocate(bfl(1:model%info%state_variable_count),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (bfl)'
   bfl = _ZERO_

   allocate(par(LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (par)'
   varid = rmbm_get_variable_id(model,varname_par,shape3d)
   if (varid.ne.-1) call rmbm_link_variable_data(model,varid,par(1:LOCATION))

   allocate(pres(LOCATION_RANGE),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (pres)'
   varid = rmbm_get_variable_id(model,varname_pres,shape3d)
   if (varid.ne.-1) call rmbm_link_variable_data(model,varid,pres(1:LOCATION))

   allocate(total(1:model%info%conserved_quantity_count))
   allocate(local(1:model%info%conserved_quantity_count))

   end subroutine init_var_gotm_rmbm
!EOC





!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set bio module environment 
!
! !INTERFACE: 
   subroutine set_env_gotm_rmbm(dt_,w_adv_ctr_,temp,salt,rho,nuh_,h_,w_,rad_,bioshade_,I_0,wnd,z_)
!
! !DESCRIPTION:
! TODO
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in) :: dt_
   integer, intent(in) :: w_adv_ctr_
   REALTYPE, intent(in),target ATTR_LOCATION_DIMENSIONS    :: temp,salt,rho,nuh_,h_,w_,rad_,bioshade_,z_
   REALTYPE, intent(in),target ATTR_LOCATION_DIMENSIONS_HZ :: I_0,wnd
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
   integer :: varid
!
!EOP
!-----------------------------------------------------------------------!


!-----------------------------------------------------------------------!
!BOC
   if (.not. rmbm_calc) return

   call rmbm_link_variable_data(model,varname_temp,temp)
   call rmbm_link_variable_data(model,varname_salt,salt)
   call rmbm_link_variable_data(model,varname_dens,rho)
   call rmbm_link_variable_data(model,varname_wind_sf,wnd)
   call rmbm_link_variable_data(model,varname_par_sf,I_0)
   
   nuh => nuh_
   h   => h_
   w   => w_
   rad => rad_
   bioshade => bioshade_
   z => z_
   
   dt = dt_
   w_adv_ctr = w_adv_ctr_

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
   IMPLICIT NONE
!
   integer, intent(in) :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer, parameter        :: adv_mode_0=0
   integer, parameter        :: adv_mode_1=1
   REALTYPE                  :: Qsour(0:nlev),Lsour(0:nlev)
   REALTYPE                  :: RelaxTau(0:nlev)
   REALTYPE                  :: dt_eff
   integer                   :: j
   integer                   :: split,posconc

!-----------------------------------------------------------------------
!BOC

   if (.not. rmbm_calc) return

   Qsour    = _ZERO_
   Lsour    = _ZERO_
   RelaxTau = 1.e15
   
   pres(1:nlev) = -z(1:nlev)
   
   sfl = _ZERO_
   call rmbm_update_air_sea_exchange(model,nlev,sfl)

   do j=1,nlev
      call rmbm_get_vertical_movement(model,j,ws(j,:))
   end do

   do j=1,model%info%state_variable_count
   
      posconc = 0
      if (model%info%variables(j)%minimum.ge._ZERO_) posconc = 1
         
!     do advection step due to settling or rising
      call adv_center(nlev,dt,h,h,ws(:,j),flux,                   &
           flux,_ZERO_,_ZERO_,w_adv_discr,adv_mode_1,cc(j,:))
         
!     do advection step due to vertical velocity
      if(w_adv_ctr .ne. 0) &
         call adv_center(nlev,dt,h,h,w,flux,                   &
              flux,_ZERO_,_ZERO_,w_adv_ctr,adv_mode_0,cc(j,:))
      
!     do diffusion step
      call diff_center(nlev,dt,cnpar,posconc,h,Neumann,Neumann,&
           sfl(j),bfl(j),nuh,Lsour,Qsour,RelaxTau,cc(j,:),cc(j,:))

   end do

   do split=1,split_factor
      dt_eff=dt/float(split_factor)

      call light_0d(nlev,bioshade_feedback)
      
      call ode_solver(ode_method,model%info%state_variable_count,nlev,dt_eff,cc,right_hand_side_rhs,right_hand_side_ppdd)
      
   end do

   end subroutine do_gotm_rmbm
!EOC


   subroutine right_hand_side_ppdd(first,numc,nlev,cc,pp,dd)
      logical, intent(in)                  :: first
      integer, intent(in)                  :: numc,nlev
      REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
      
      REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
      REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
      
      integer :: ci

      pp = _ZERO_
      dd = _ZERO_
      
      ! Iterate over all depth levels
      do ci=1,nlev
         call rmbm_do(model,ci,pp(:,:,ci),dd(:,:,ci),work_cc_diag(:,ci))
      end do
      
      if (first) then
         ! First time during this time step that do_bio is called: store diagnostic values.
         ! NB. higher order integration schemes may call this routine multiple times.
         ! In that case only the value at the first call is used (essential because time
         ! integration is done internally).
         
         do ci=1,model%info%diagnostic_variable_count
            if (model%info%diagnostic_variables(ci)%time_treatment.eq.0) then
               ! Simply use last value
               cc_diag(ci,1:nlev) = work_cc_diag(ci,1:nlev)
            else
               ! Integration or averaging in time needed: for now do simple Forward Euler integration.
               cc_diag(ci,1:nlev) = cc_diag(ci,1:nlev) + work_cc_diag(ci,1:nlev)*dt
            end if
         end do
      end if

   end subroutine

   subroutine right_hand_side_rhs(first,numc,nlev,cc,rhs)
      logical, intent(in)                  :: first
      integer, intent(in)                  :: numc,nlev
      REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
      
      REALTYPE, intent(out)                :: rhs(1:numc,0:nlev)

      integer :: ci

      ! Initialization is needed because the different biogeochemical models increment or decrement
      ! the temporal derivatives, rather than setting them directly. Tis is needed for the simultaenous
      ! running of different coupled BGC models.
      rhs = _ZERO_

      ! Iterate over all depth levels
      do ci=1,nlev
         call rmbm_do(model,ci,rhs(:,ci),work_cc_diag(:,ci))
      end do
      
      if (first) then
         ! First time during this time step that do_bio is called: store diagnostic values.
         ! NB. higher order integration schemes may call this routine multiple times.
         ! In that case only the value at the first call is used (essential because time
         ! integration is done internally).
         
         do ci=1,model%info%diagnostic_variable_count
            if (model%info%diagnostic_variables(ci)%time_treatment.eq.0) then
               ! Simply use last value
               cc_diag(ci,1:nlev) = work_cc_diag(ci,1:nlev)
            else
               ! Integration or averaging in time needed: for now do simple Forward Euler integration.
               cc_diag(ci,1:nlev) = cc_diag(ci,1:nlev) + work_cc_diag(ci,1:nlev)*dt
            end if
         end do
      end if

   end subroutine



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

!  internal arrays
   if (allocated(par))            deallocate(par)
   if (allocated(cc))             deallocate(cc)
   if (allocated(cc_diag))        deallocate(cc_diag)
   if (allocated(work_cc_diag))   deallocate(work_cc_diag)
   if (allocated(ws))             deallocate(ws)
   if (allocated(sfl))            deallocate(sfl)
   if (allocated(bfl))            deallocate(bfl)
   if (allocated(total))          deallocate(total)
   if (allocated(local))          deallocate(local)
   LEVEL1 'done.'

   end subroutine clean_gotm_rmbm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate light over entire column
!
! !INTERFACE:
   subroutine light_0d(nlev,bioshade_feedback)
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
! !LOCAL VARIABLES:
   integer :: i
   REALTYPE :: zz,bioext,localext

!EOP
!-----------------------------------------------------------------------
!BOC
   zz = _ZERO_
   bioext = _ZERO_
   do i=nlev,1,-1
      localext = rmbm_get_bio_extinction(model,i)
   
      ! Add the extinction of the first half of the grid box.
      bioext = bioext+localext*0.5*h(i)

      zz=zz+0.5*h(i)
      par(i)=rad(nlev)*(_ONE_-A)*exp(-zz/g2-bioext)

      ! Add the extinction of the second half of the grid box.
      bioext = bioext+localext*0.5*h(i)
      
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade(i)=exp(-bioext)
   end do

   end subroutine light_0d
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

   IMPLICIT NONE

#ifdef NETCDF_FMT
#include "netcdf.inc"
#endif
!
! !LOCAL VARIABLES:
   integer :: iret,n
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (.not. rmbm_calc) return
   
   select case (out_fmt)
      case (NETCDF)
#ifdef NETCDF_FMT
         iret = define_mode(ncid,.true.)

         dims(1) = lon_dim
         dims(2) = lat_dim
         dims(3) = z_dim
         dims(4) = time_dim

         ! Add a variable for each prognostic variable
         do n=1,model%info%state_variable_count
            iret = new_nc_variable(ncid,model%info%variables(n)%name,NF_REAL, &
                                   4,dims,model%info%variables(n)%id)
            iret = set_attributes(ncid,model%info%variables(n)%id,       &
                                  units=model%info%variables(n)%units,    &
                                  long_name=model%info%variables(n)%longname)
         end do

         ! Add a variable for each diagnostic variable
         do n=1,model%info%diagnostic_variable_count
            iret = new_nc_variable(ncid,model%info%diagnostic_variables(n)%name,NF_REAL, &
                                   4,dims,model%info%diagnostic_variables(n)%id)
            iret = set_attributes(ncid,model%info%diagnostic_variables(n)%id,       &
                                  units=model%info%diagnostic_variables(n)%units,    &
                                  long_name=model%info%diagnostic_variables(n)%longname)
         end do

         dims(3) = time_dim

         ! Add a variable for each conserved quantity
         do n=1,model%info%conserved_quantity_count
            iret = new_nc_variable(ncid,trim(model%info%conserved_quantities(n)%name)//'_tot',NF_REAL, &
                                   3,dims,model%info%conserved_quantities(n)%id)
            iret = set_attributes(ncid,model%info%conserved_quantities(n)%id,       &
                                  units='m*'//model%info%conserved_quantities(n)%units,    &
                                  long_name=trim(model%info%conserved_quantities(n)%longname)//', depth-integrated')
         end do

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
! !LOCAL VARIABLES:
   integer :: iret,ilev,n
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (.not. rmbm_calc) return
   
   select case (out_fmt)
      case (NETCDF)
#ifdef NETCDF_FMT
         ! Store biogeochemical state (prognostic) variables.
         do n=1,model%info%state_variable_count
            iret = store_data(ncid,model%info%variables(n)%id,XYZT_SHAPE,nlev,array=cc(n,0:nlev))
         end do

         ! Time-average and store diagnostic variables.
         do n=1,model%info%diagnostic_variable_count
            if (model%info%diagnostic_variables(n)%time_treatment==2) &
               cc_diag(n,1:nlev) = cc_diag(n,1:nlev)/(nsave*dt)
            iret = store_data(ncid,model%info%diagnostic_variables(n)%id,XYZT_SHAPE,nlev,array=cc_diag(n,0:nlev))
            if (model%info%diagnostic_variables(n)%time_treatment==2 .or. &
                model%info%diagnostic_variables(n)%time_treatment==3) &
               cc_diag(n,1:nlev) = _ZERO_
         end do

         ! Integrate conserved quantities over depth
         total = _ZERO_
         do ilev=1,nlev
            call rmbm_get_conserved_quantities(model,ilev,local)
            total = total + h(ilev+1)*local
         end do

         ! Store conserved quantity integrals
         do n=1,model%info%conserved_quantity_count
            iret = store_data(ncid,model%info%conserved_quantities(n)%id,XYT_SHAPE,1,scalar=total(n))
         end do
#endif
   end select

   end subroutine save_gotm_rmbm
!EOC

   end module gotm_rmbm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------

