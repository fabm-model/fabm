!$Id: bio.F90,v 1.53 2009-11-20 08:16:56 kb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio --- biological model \label{sec:bio}
!
! !INTERFACE:
   module bio
!
! !DESCRIPTION:
! This is the central module for all biogeochemical models. 
! From here, after reading the namelist file {\tt bio.nml},
! the individual biogeochemical model is initialised, the memory
! is allocated, the advection and diffusion is called, the ODE solvers
! for the right hand sides are called, and simple Lagrangian particle
! calculations are managed.
! 
! !USES:
   use bio_var
#ifndef NO_0D_BIO
   use bio_0d, only: model,init_bio_0d, init_var_0d, &
                     light_0d, light_0d_par, &
                     surface_fluxes_0d, update_sinking_rates_0d, &
                     do_bio_0d_eul, do_bio_0d_eul_rhs, do_bio_0d_par
#endif

#ifdef BIO_TEMPLATE
   use bio_template, only : init_bio_template,init_var_template
   use bio_template, only : light_template
#endif

#ifdef BIO_NPZD
   use bio_npzd, only : init_bio_npzd,init_var_npzd
   use bio_npzd, only : light_npzd, do_bio_npzd
#endif

#ifdef BIO_CL
   use bio_cl, only : init_bio_cl,init_var_cl
   use bio_cl, only : do_bio_cl
#endif

#ifdef BIO_IOW
   use bio_iow, only : init_bio_iow,init_var_iow
   use bio_iow, only : light_iow,surface_fluxes_iow,do_bio_iow
#endif

#ifdef BIO_FASHAM
   use bio_fasham, only : init_bio_fasham,init_var_fasham
   use bio_fasham, only : light_fasham,do_bio_fasham
#endif

#ifdef BIO_SED
   use bio_sed, only : init_bio_sed,init_var_sed
   use bio_sed, only : do_bio_sed_eul,do_bio_sed_par
#endif

#ifdef BIO_MAB
   use bio_mab, only : init_bio_mab,init_var_mab
   use bio_mab, only : light_mab,surface_fluxes_mab,do_bio_mab
#endif

#ifdef BIO_ROLM
   use bio_rolm, only : init_bio_rolm,init_var_rolm
   use bio_rolm, only : light_rolm,surface_fluxes_rolm,do_bio_rolm
#endif

#ifdef BIO_NPZD_FE
   use bio_npzd_fe, only : init_bio_npzd_fe,init_var_npzd_fe
   use bio_npzd_fe, only : light_npzd_fe,surface_fluxes_npzd_fe,do_bio_npzd_fe
#endif

#ifdef BIO_PHOTO
   use bio_photo, only : init_bio_photo,init_var_photo
   use bio_photo, only : do_bio_photo_eul,do_bio_photo_par
#endif

#if 0
   use mussels, only : init_mussels, do_mussels, end_mussels
   use mussels, only : mussels_calc,total_mussel_flux
#endif

   use output, only: out_fmt,write_results,ts

   use util

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
   public init_bio,init_var_bio
   public set_env_bio,do_bio,get_bio_updates
   public clean_bio


   logical, public                     :: bio_calc=.false.
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio.F90,v $
!  Revision 1.53  2009-11-20 08:16:56  kb
!  allow compilation excluding 0D BIO framework - set in Rules.make
!
!  Revision 1.52  2009-11-11 13:08:54  kb
!  added chlorination model - Rennau
!
!  Revision 1.50  2009-10-21 08:02:07  hb
!  Fluff layer resuspension added.
!
!  Revision 1.49  2008-11-11 13:40:32  jorn
!  major revision of 0d biogeochemical framework; added output of depth-integrated conserved BGC quantities; added support for running multiple BGC models side-by-side
!
!  Revision 1.48  2008-11-06 13:42:44  jorn
!  several changes to 0d framework for biogeochemical models; added explicit support for time- and space-varying sinking and light extinction
!
!  Revision 1.47  2008-11-05 12:51:38  jorn
!  restructured 0d biogeochemical framework; added experimental support for self-shading in 0d-based particle models
!
!  Revision 1.46  2008-11-03 12:57:34  jorn
!  added support for 0D biogeochemical models in Lagrangian mode
!
!  Revision 1.45  2008-10-31 11:10:32  jorn
!  first version of bio framework for 0D models and 0D NPZD test case
!
!  Revision 1.44  2008-07-08 10:09:05  lars
!  new structure with general particle support
!
!  Revision 1.43  2008-03-26 08:56:53  kb
!  new directory based bio structure
!
!  Revision 1.42  2008-02-20 11:29:59  kb
!  added NPZD iron model - Weber et. all + Inga Hense
!
!  Revision 1.40  2007-11-07 11:14:24  kb
!  no mussesl yet
!
!  Revision 1.39  2007-10-01 12:44:06  kbk
!  added RedOxLayer Model (ROLM)
!
!  Revision 1.38  2007-04-18 07:36:47  kbk
!  mussels will be developed in 4.1.x
!
!  Revision 1.37  2007-04-18 06:57:36  kbk
!  Lagrangian simulations disabled by default
!
!  Revision 1.36  2007-03-14 12:46:07  kbk
!  proper cleaning after simulation
!
!  Revision 1.35  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.34  2007-01-04 12:54:12  hb
!  ifdef LAGRANGE removed
!
!  Revision 1.33  2006-11-17 07:13:17  kbk
!  rho amd wind-speed available via bio_var
!
!  Revision 1.32  2006-11-12 19:42:44  hb
!  vertical advection due to physical vertical velocities enabled for the bio module
!
!  Revision 1.31  2006-11-06 13:36:46  hb
!  Option for conservative vertical advection added to adv_center
!
!  Revision 1.30  2006-10-26 13:12:46  kbk
!  updated bio models to new ode_solver
!
!  Revision 1.29  2005-12-27 11:23:04  hb
!  Weiss 1970 formula now used for surface oxygen saturation calculation in bio_mab.F90
!
!  Revision 1.28  2005-12-27 06:51:49  hb
!  New biomodel bio_mab (bio_iow with additional sediment equation) added
!
!  Revision 1.27  2005-12-02 20:57:27  hb
!  Documentation updated and some bugs fixed
!
!  Revision 1.26  2005-11-18 10:59:35  kbk
!  removed unused variables - some left in parameter lists
!
!  Revision 1.25  2005/11/17 09:58:18  hb
!  explicit argument for positive definite variables in diff_center()
!
!  Revision 1.24  2005/10/11 08:43:44  lars
!  checked new transport routines
!
!  Revision 1.23  2005/09/19 21:07:00  hb
!  yevol replaced by adv_center and diff_center
!
!  Revision 1.22  2005/09/12 14:48:33  kbk
!  merged generic biological module support
!
!  Revision 1.21.2.1  2005/07/06 09:00:19  hb
!  moved bio_save() from do_bio() to time_loop - temporary no NPZD totn calculation
!
!  Revision 1.21  2004/08/18 11:34:14  hb
!  zlev now allocated from 0 to nlev
!
!  Revision 1.20  2004/08/02 11:44:12  kbk
!  bio module compiles and runs with GETM
!
!  Revision 1.19  2004/08/02 08:35:08  hb
!  no need to pass time information
!
!  Revision 1.18  2004/08/01 15:54:49  hb
!  call to light_fasham commented in again
!
!  Revision 1.17  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.16  2004/07/28 11:34:29  hb
!  Bioshade feedback may now be switched on or off, depending on bioshade_feedback set to .true. or .false. in bio.nml
!
!  Revision 1.15  2004/06/29 08:03:16  hb
!  Fasham et al. 1990 model implemented
!
!  Revision 1.14  2004/05/28 13:24:49  hb
!  Extention of bio_iow to fluff layer and surface nutrient fluxes
!
!  Revision 1.13  2004/04/13 09:18:54  kbk
!  size and temperature dependend filtration rate
!
!  Revision 1.12  2004/03/31 12:58:52  kbk
!  lagrangian solver uses - total_mussel_flux
!
!  Revision 1.11  2004/03/30 11:32:48  kbk
!  select between eulerian or lagrangian solver
!
!  Revision 1.10  2003/12/11 09:58:22  kbk
!  now compiles with FORTRAN_COMPILER=IFORT - removed TABS
!
!  Revision 1.9  2003/10/28 10:22:45  hb
!  added support for sedimentation only 1 compartment bio model
!
!  Revision 1.8  2003/10/16 15:42:16  kbk
!  simple mussesl model implemented - filter only
!
!  Revision 1.7  2003/10/14 08:00:09  hb
!  initialise sfl - no special treatment when cc(,) < 0
!
!  Revision 1.6  2003/09/16 12:11:24  hb
!  added new biological model - bio_iow
!
!  Revision 1.5  2003/07/23 12:27:31  hb
!  more generic support for different bio models
!
!  Revision 1.3  2003/04/05 07:01:41  kbk
!  moved bioshade variable to meanflow - to compile properly
!
!  Revision 1.2  2003/04/04 14:25:52  hb
!  First iteration of four-compartment geobiochemical model implemented
!
!  Revision 1.1  2003/04/01 17:01:00  hb
!  Added infrastructure for geobiochemical model
!
!EOP
!-----------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
!  from a namelist
   REALTYPE                  :: cnpar=0.9
   integer                   :: w_adv_discr=6
   integer                   :: ode_method=1
   integer                   :: split_factor=1
   logical                   :: bioshade_feedback=.true.


   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio(namlst,fname,unit,nmax_)
!
! !DESCRIPTION: 
! Here, all initialization procedures are triggered that
! are independent of the use of the BIO module from inside GOTM, or
! from an external calling program, i.e.\ from a 3D model. The
! following steps are subsequently performed.
!
! \begin{enumerate}
!  \item Memory for {\tt nmax} vertical layers is allocated for all 
!  Eulerian fields. If called from an external 3D model, {\tt nmax} 
!  corresponds to the maximum number of layers in \emph{any} of the
!  water columns.
!  \item The namelist {\tt bio.nml} is read. 
!  \item The initialization routine for the selected ecosystem model 
!  is called. In most cases this amounts to not much more than reading 
!  the corresponding namelist, e.g.\ {\tt bio\_npzd.nml}, and assigning
!  the description tags to the ecosystem variables.
! \end{enumerate}
!
! After the initialization is finished, some information about
! the selected ODE-solver is written to the screen.
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
!  Original author(s): Lars Umlauf, Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!  local variables
   integer                   :: rc,i,j,n
   namelist /bio_nml/ bio_calc,bio_model,bio_eulerian,                  &
                      cnpar,w_adv_discr,ode_method,split_factor,        &
                      bioshade_feedback,npar
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_bio'

!  open and read the namelist
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_nml,err=99)
   close(namlst)

   if (bio_calc) then
   
!     read individual model namelists
      select case (bio_model)
         
      case (-1)
#ifdef BIO_TEMPLATE

         call init_bio_template(namlst,'bio_template.nml',unit)
         
#else
         FATAL "bio_template not compiled!!"
         FATAL "set environment variable BIO_TEMPLATE different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (1) ! The NPZD model
#ifdef BIO_NPZD

         call init_bio_npzd(namlst,'bio_npzd.nml',unit)

#else
         FATAL "bio_npzd not compiled!!"
         FATAL "set environment variable BIO_NPZD different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (2)  ! The IOW model
#ifdef BIO_IOW

         call init_bio_iow(namlst,'bio_iow.nml',unit)

#else
         FATAL "bio_iow not compiled!!"
         FATAL "set environment variable BIO_IOW different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (3)  ! The simple sedimentation model
#ifdef BIO_SED

         call init_bio_sed(namlst,'bio_sed.nml',unit)

#else
         FATAL "bio_sed not compiled!!"
         FATAL "set environment variable BIO_SED different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (4)  ! The FASHAM model
#ifdef BIO_FASHAM

         call init_bio_fasham(namlst,'bio_fasham.nml',unit)

#else
         FATAL "bio_fasham not compiled!!"
         FATAL "set environment variable BIO_FASHAM different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (5)  ! The IOW model, modified for MaBenE
#ifdef BIO_MAB

         call init_bio_mab(namlst,'bio_mab.nml',unit)

#else
         FATAL "bio_mab not compiled!!"
         FATAL "set environment variable BIO_MAB different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (6)  ! RedOxLayer Model (ROLM)
#ifdef BIO_ROLM

         call init_bio_rolm(namlst,'bio_rolm.nml',unit)

#else
         FATAL "bio_rolm not compiled!!"
         FATAL "set environment variable BIO_ROLM different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (7)  ! NPZD_FE (Weber_etal2007)
#ifdef BIO_NPZD_FE

         call init_bio_npzd_fe(namlst,'bio_npzd_fe.nml',unit)

#else
         FATAL "bio_bio_npzd_fe not compiled!!"
         FATAL "set environment variable BIO_NPZD_FE different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (8) ! The CL model
#ifdef BIO_CL

         call init_bio_cl(namlst,'bio_cl.nml',unit)

#else
         FATAL "bio_cl not compiled!!"
         FATAL "set environment variable BIO_CL different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif

      case (20)  ! PHOTO (adaption model according to Nagai et al. (2003)
#ifdef BIO_PHOTO

         call init_bio_photo(namlst,'bio_photo.nml',unit)

#else
         FATAL "bio_bio_photo not compiled!!"
         FATAL "set environment variable BIO_PHOTO different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif

#ifndef NO_0D_BIO
      case (1000:)
         call init_bio_0d(namlst,unit)
#endif
   
      case default
         stop "bio: no valid biomodel specified in bio.nml !"
      end select


!     report variable descriptions
      LEVEL2 'updated fields will be:'
      do n=1,numc
         LEVEL3 trim(var_names(n)),'  ',trim(var_units(n)), &
                '  ',trim(var_long(n))
      end do




!     allocate memory for Eulerian fields
      if (nmax_ .gt. 0) then
         nlev = nmax_
         nmax = nmax_
         call alloc_eul
      else
         FATAL 'nmax>0 for successful memory allocation'
         stop 'init_bio()'
      endif

!     report type of solver 
      if ( bio_eulerian ) then
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
               stop "bio: no valid ode_method specified in bio.nml!"
         end select
      else
         LEVEL2 "Using particle solver"
      end if

#if 0
!     Initialise 'mussels' module
      call init_mussels(namlst,"mussels.nml",unit,nlev,h)
#endif
   end if


   return

98 LEVEL2 'I could not open bio.nml'
   LEVEL2 'If thats not what you want you have to supply bio.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio.nml'
   bio_calc = .false.
   return
99 FATAL 'I could not read bio.nml'
   stop 'init_bio'
   end subroutine init_bio
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise bio variables
!
! !INTERFACE:
   subroutine init_var_bio
!
! !DESCRIPTION:
! A call to this routine  initializes the variables of the bio module
! with meaningful values, depending on the water column properties set
! by a previous call to {\tt set\_env\_bio()}. The steps taken here
! include
! \begin{enumerate}
!  \item Allocating memory for the particle properties (if a particle
!  solver has been chosen),
!  \item Computing the position of the grid interfaces (needed for 
!  averaging particle properties over grid cells)
!  \item Calling the initialization routines for the respective 
!  bio modules
! \end{enumerate}
!  
! If the bio module is used from a 3D code outside GOTM this routine 
! should not be called. In this case, all initialization should be 
! done inside the external calling program. Note in particular that 
! the numbe of particles may change during the run, and variable memory
! needs to be allocated.
!
! !USES:
   IMPLICIT NONE
!
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Lars Umlauf, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!  local variables
   integer                   :: rc,i,j,n

!-----------------------------------------------------------------------
!BOC

!  allocate initial memory for particles
   if (.not. bio_eulerian) then

      call alloc_par()

   end if


!  update grid
   zlev(0)= zbot

   do n=1,nlev
      zlev(n)=zlev(n-1)+h(n)
   end do

   ztop         = zlev(nlev)


!  initialize individual models
   if (bio_calc) then

      select case (bio_model)
         
      case (-1)
#ifdef BIO_TEMPLATE
         
         call init_var_template

#else
         FATAL "bio_template not compiled!!"
         FATAL "set environment variable BIO_TEMPLATE different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (1) ! The NPZD model
#ifdef BIO_NPZD

         call init_var_npzd

#else
         FATAL "bio_npzd not compiled!!"
         FATAL "set environment variable BIO_NPZD different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (2)  ! The IOW model
#ifdef BIO_IOW

         call init_var_iow

#else
         FATAL "bio_iow not compiled!!"
         FATAL "set environment variable BIO_IOW different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (3)  ! The simple sedimentation model
#ifdef BIO_SED

         call init_var_sed

#else
         FATAL "bio_sed not compiled!!"
         FATAL "set environment variable BIO_SED different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (4)  ! The FASHAM model
#ifdef BIO_FASHAM

         call init_var_fasham

#else
         FATAL "bio_fasham not compiled!!"
         FATAL "set environment variable BIO_FASHAM different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (5)  ! The IOW model, modified for MaBenE
#ifdef BIO_MAB

         call init_var_mab

#else
         FATAL "bio_mab not compiled!!"
         FATAL "set environment variable BIO_MAB different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (6)  ! RedOxLayer Model (ROLM)
#ifdef BIO_ROLM

         call init_var_rolm

#else
         FATAL "bio_rolm not compiled!!"
         FATAL "set environment variable BIO_ROLM different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (7)  ! NPZD_FE (Weber_etal2007)
#ifdef BIO_NPZD_FE

         call init_var_npzd_fe

#else
         FATAL "bio_bio_npzd_fe not compiled!!"
         FATAL "set environment variable BIO_NPZD_FE different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (8) ! The CL model
#ifdef BIO_CL

         call init_var_cl

#else
         FATAL "bio_cl not compiled!!"
         FATAL "set environment variable BIO_CL different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif
      case (20)  ! PHOTO (adaption model according to Nagai et al. (2003)
#ifdef BIO_PHOTO

         call init_var_photo

#else
         FATAL "bio_bio_photo not compiled!!"
         FATAL "set environment variable BIO_PHOTO different from false"
         FATAL "and re-compile"
         stop "init_bio()"
#endif

#ifndef NO_0D_BIO
      case (1000:)
         call init_var_0d
#endif
         
      case default
         stop "bio: no valid biomodel specified in bio.nml !"
      end select

   end if


   return

98 LEVEL2 'I could not open bio.nml'
   LEVEL2 'If thats not what you want you have to supply bio.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio.nml'
   bio_calc = .false.
   return
99 FATAL 'I could not read bio.nml'
   stop 'init_bio'
   end subroutine init_var_bio
!EOC





!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set bio module environment 
!
! !INTERFACE: 
   subroutine set_env_bio(nlev_,dt_,zbot_,u_taub_,h_,t_,s_,rho_,nuh_,rad_,   &
                          wind_,I_0_,secondsofday_,w_,w_adv_ctr_,npar_)
!
! !DESCRIPTION:
! This routine prepares the environment for the bio module, and keeps 
! track of the memory. Every time this routine is called 
! the module's variables related to meteorology, mixing, grid parameters, 
! etc. are updated with the values supplied as arguments. These updated
! values are then accessible to all module routines if required.
!
! Note that if the bio-module is called from outside GOTM,
! e.g.\ from a 3-D model,
! this routine has to be called every time the environment parameters,
! the number of vertical grid points, or the number of particles in the 
! local water column has changed.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer , intent(in)                :: nlev_
   REALTYPE, intent(in)                :: dt_
   REALTYPE, intent(in)                :: zbot_
   REALTYPE, intent(in)                :: h_(0:nlev)
   REALTYPE, intent(in)                :: u_taub_
   REALTYPE, intent(in)                :: nuh_(0:nlev)
   REALTYPE, intent(in)                :: t_(0:nlev)
   REALTYPE, intent(in)                :: s_(0:nlev)
   REALTYPE, intent(in)                :: rho_(0:nlev)
   REALTYPE, intent(in)                :: rad_(0:nlev)
   REALTYPE, intent(in)                :: wind_
   REALTYPE, intent(in)                :: I_0_
   integer , intent(in)                :: secondsofday_
   REALTYPE, optional, intent(in)      :: w_(0:nlev)
   integer , optional, intent(in)      :: w_adv_ctr_
   integer , optional, intent(in)      :: npar_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Lars Umlauf, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------!


!-----------------------------------------------------------------------!
!BOC

!  assign mandatory arguments
   nlev         = nlev_
   dt           = dt_
   zbot         = zbot_
   h            = h_
   u_taub       = u_taub_
   t            = t_
   s            = s_
   rho          = rho_
   nuh          = nuh_
   rad          = rad_
   wind         = wind_
   I_0          = I_0_
   secondsofday = secondsofday_


!  assign optional arguments

   if (present(w_)        )  w         = w_
   if (present(w_adv_ctr_))  w_adv_ctr = w_adv_ctr_

!  check for sufficient memory
   if (nlev > nmax) then
      FATAL 'nlev=', nlev, ' but memory is allocated only for nmax=', nmax
      FATAL 'call "init\_bio()" with nmax >= nlev'
      stop 'set_env_bio()'
   end if


!  keep track of particle memory
   if (present(npar_) ) then

      if (npar /= npar_) then
         
         npar      = npar_
         par_allocation = .true.
         
      else
         par_allocation = .false.
      end if

   elseif (.not. bio_eulerian) then

      FATAL 'npar not specified'
      stop 'set_env_bio('

   endif
      

   return
   end subroutine set_env_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update the bio model 
!
! !INTERFACE:
   subroutine do_bio()
!
! !DESCRIPTION:
! This routine is a simple wrapper selecting between calls to the
! Eulerian and Lagrangian (particle) routines used to update the 
! bio model. 
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES
   integer                             :: n

!-----------------------------------------------------------------------
!BOC

!  update grid
   zlev(0)= zbot

   do n=1,nlev
      zlev(n)=zlev(n-1)+h(n)
   end do

   ztop         = zlev(nlev)

!  update model fields
   if (bio_calc) then
      

      ! Eulerian model
      if (bio_eulerian) then
         
         call do_bio_eul
         

      ! particle model
      else 
         
         if (npar /= 0) call do_bio_par
         
      end if
      
#if 0
      if (mussels_calc) then
         call do_mussels(numc,dt,t(1))
      end if
#endif
      
   end if

   return
   end subroutine do_bio
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update the Eulerian bio model \label{sec:do-bio}
!
! !INTERFACE:
   subroutine do_bio_eul()
!
! !DESCRIPTION:
! This is the main loop for the Eulerian biogeochemical model. Basically 
! an operational split method is used, with first calculating the
! transport part, and than the reaction part.
! During the transport part, all sinks and sources are set to zero,
! and the surface fluxes are computed by calling the
! model specific surface flux subroutine. Then the mussel module
! is called.  For the Eulerian calculation, vertical advection
! (due to settling or rising or vertical migration), vertical advection due
! to physical velocity and vertical
! diffusion (due to mixing) and afterwards the light 
! calculation (for the PAR) and the ODE solver for the right
! hand sides are called. The vertical advection due to settling and
! rising must be conservative,
! which is ensured by setting the local variable {\tt adv\_mode\_1=1},
! see section \ref{sec:advectionMean} on page \pageref{sec:advectionMean}.
! In contrast to this, the vertical advection due to physical velocities must be
! non-conservative, such that for that the local variable {\tt adv\_mode\_0}
! is set to 0, see  see section \ref{sec:advectionMean} on page
! \pageref{sec:advectionMean}.
! It should be noted here that the PAR and the selfshading effect
! is calculated in a similar way for all biogeochemical models
! implemented in GOTM so far. In the temperature equation the
! absorption of solar radiation, $I(z)$, is the only source term,
! see equation (\ref{Iz}) section \ref{sec:temperature}.
! In (\ref{Iz}), a term $B(z)$ due to bioturbidity is used, which 
! is calculated as a function of the biogeochemical particulate
! matter in the water column:
! \begin{equation}\label{B}
! B(z)=\exp\left(-k_c\int_z^0\left(\sum C_{turb}(\xi)\right)\,d\xi\right),
! \end{equation}
! where $k_c$ is the attenuation constant for self shading and 
! $\sum C_{turb}$ is the sum of the biogeochemical particulate 
! matter concentrations.
! The photosynthetically
! available radiation, $I_{PAR}$, follows from
! \begin{equation}
!   \label{light}
!   I_{PAR}(z)=I_0
! (1-a)\exp\left(\frac{z}{\tilde\eta_2}\right)
!   B(z).
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
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
   integer                   :: j,n
   integer                   :: split
   integer                   :: i,np,nc,vars_zero_d=0
   integer, save             :: count=0
   logical, save             :: set_C_zero=.true.

!-----------------------------------------------------------------------
!BOC

   Qsour    = _ZERO_
   Lsour    = _ZERO_
   RelaxTau = 1.e15
   
   select case (bio_model)
   case (-1)
   case (1)
   case (2)
#ifdef BIO_IOW
      call surface_fluxes_iow(nlev,t(nlev),s(nlev),wind)
      if (numc.eq.10) vars_zero_d=1
#endif
   case (3)
   case (4)
   case (5)
#ifdef BIO_MAB
      call surface_fluxes_mab(nlev,t(nlev),s(nlev))
#endif
   case (6)
#ifdef BIO_ROLM
      call surface_fluxes_rolm(nlev,t(nlev))
#endif
   case (7)
#ifdef BIO_NPZD_FE
      call surface_fluxes_npzd_fe(nlev)
#endif
   case (8)
#ifndef NO_0D_BIO
   case (1000:)
      call surface_fluxes_0d(nlev)
      call update_sinking_rates_0d(numc,nlev,cc)
#endif
   end select

   do j=1,numc-vars_zero_d
         
!     do advection step due to settling or rising
      call adv_center(nlev,dt,h,h,ws(j,:),flux,                   &
           flux,_ZERO_,_ZERO_,w_adv_discr,adv_mode_1,cc(j,:))
         
!     do advection step due to vertical velocity
      if(w_adv_ctr .ne. 0) then
         call adv_center(nlev,dt,h,h,w,flux,                   &
              flux,_ZERO_,_ZERO_,w_adv_ctr,adv_mode_0,cc(j,:))
      end if
      
!     do diffusion step
      call diff_center(nlev,dt,cnpar,posconc(j),h,Neumann,Neumann,&
           sfl(j),bfl(j),nuh,Lsour,Qsour,RelaxTau,cc(j,:),cc(j,:))

   end do

   do split=1,split_factor
      dt_eff=dt/float(split_factor)

!     very important for 3D models to save extra 3D field:
      bioshade_=_ONE_

      select case (bio_model)
      case (-1)
#ifdef BIO_TEMPLATE
         call light_template(nlev,bioshade_feedback)
         !               call ode_solver(ode_method,numc,nlev,dt_eff,cc,do_bio_template)
#endif
      case (1)
#ifdef BIO_NPZD
         call light_npzd(nlev,bioshade_feedback)
#endif
      case (2)
#ifdef BIO_IOW
         call light_iow(nlev,bioshade_feedback)
#endif
      case (3)
      case (4)
#ifdef BIO_FASHAM
         call light_fasham(nlev,bioshade_feedback)
#endif
      case (5)
#ifdef BIO_MAB
         call light_mab(nlev,bioshade_feedback)
#endif
      case (6)
#ifdef BIO_ROLM
         call light_rolm(nlev,bioshade_feedback)
#endif
      case (7)
#ifdef BIO_NPZD_FE
         call light_npzd_fe(nlev,bioshade_feedback)
#endif
      case (8)
#ifndef NO_0D_BIO
      case (1000:)
         call light_0d(nlev,bioshade_feedback)
#endif
      end select
      
      call ode_solver(ode_method,numc,nlev,dt_eff,cc,right_hand_side_rhs,right_hand_side_ppdd)
      
   end do

   return
   end subroutine do_bio_eul
!EOC


   subroutine right_hand_side_ppdd(first,numc,nlev,cc,pp,dd)
      logical, intent(in)                  :: first
      integer, intent(in)                  :: numc,nlev
      REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
      
      REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
      REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
      
      REALTYPE,allocatable                 :: rhs(:,:)

      select case (bio_model)
      case (-1)
      case (1)
#ifdef BIO_NPZD
         call do_bio_npzd(first,numc,nlev,cc,pp,dd)
#endif
      case (2)
#ifdef BIO_IOW
         call do_bio_iow(first,numc,nlev,cc,pp,dd)
#endif
      case (3)
      case (4)
#ifdef BIO_FASHAM
         call do_bio_fasham(first,numc,nlev,cc,pp,dd)
#endif
      case (5)
#ifdef BIO_MAB
         call do_bio_mab(first,numc,nlev,cc,pp,dd)
#endif
      case (6)
#ifdef BIO_ROLM
         call do_bio_rolm(first,numc,nlev,cc,pp,dd)
#endif
      case (7)
#ifdef BIO_NPZD_FE
         call do_bio_npzd_fe(first,numc,nlev,cc,pp,dd)
#endif
      case (8)
#ifdef BIO_CL
         call do_bio_cl(first,numc,nlev,cc,pp,dd)
#endif
#ifndef NO_0D_BIO
      case (1000:)
         call do_bio_0d_eul(first,numc,nlev,cc,pp,dd)
#endif
      case default
         allocate(rhs(1:numc,0:nlev))
         call right_hand_side_rhs(first,numc,nlev,cc,rhs)
         pp = _ZERO_
         dd = _ZERO_
         do i=1,numc
            pp(i,i,:) = rhs(i,:)
         end do
         deallocate(rhs)
      end select
   end subroutine

   subroutine right_hand_side_rhs(first,numc,nlev,cc,rhs)
      logical, intent(in)                  :: first
      integer, intent(in)                  :: numc,nlev
      REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
      
      REALTYPE, intent(out)                :: rhs(1:numc,0:nlev)

      REALTYPE, allocatable                :: pp(:,:,:),dd(:,:,:)

      select case (bio_model)
#ifndef NO_0D_BIO
      case (1000:)
         call do_bio_0d_eul_rhs(first,numc,nlev,cc,rhs)
#endif
      case default
         allocate(pp(1:numc,1:numc,0:nlev))
         allocate(dd(1:numc,1:numc,0:nlev))
         call right_hand_side_ppdd(first,numc,nlev,cc,pp,dd)
         do i=1,numc
            rhs(i,:) = _ZERO_
            do j=1,numc
               rhs(i,:) = rhs(i,:) + pp(i,j,:)-dd(i,j,:)
            end do
         end do
         deallocate(pp)
         deallocate(dd)
      end select
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update the bio model for particles
!
! !INTERFACE:
   subroutine do_bio_par()
!
! !DESCRIPTION:
! Here, all particle properties are updated. This includes a
! random-walk step in order to update the particle positions, 
! as well as a modification of bio-geochemical properties 
! according to the particle model. 
!
! At the moment, only two particle models are available: model 3 for
! simple sedimentation and model 20 corresponding to a simple
! photo-adaption model.
!
! !USES:
   IMPLICIT NONE
!
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer                        :: nt

!-----------------------------------------------------------------------
!BOC

   
   if (par_allocation) then 
      call dealloc_par
      call alloc_par
   endif



   ! stochastic dispersion of all particle types
   do nt=1,ntype

      call lagrange(nlev,dt,zlev,nuh,ws(1,1),npar,                      &
           par_act(:,nt),                                               &
           par_ind(:,nt),                                               &
           par_z  (:,nt)  )


   end do


   ! update particle properties according to model
   select case (bio_model)
   case (3)
      call do_bio_sed_par
   case (20)
      call do_bio_photo_par
#ifndef NO_0D_BIO
   case (1000:)
      call light_0d_par(model%models(1),nlev,bioshade_feedback)
      call do_bio_0d_par(ode_method,dt)
#endif
   case default
      FATAL 'bio_model=', bio_model, ' is not a valid particle model.'
      stop 'do_bio_par()'
   end select


   return
   end subroutine do_bio_par
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Return updated bio variable
!
! !INTERFACE: 
   subroutine get_bio_updates(nlev,bioshade)
!
! !DESCRIPTION:
! Return variables that have been updated by the bio module. At the 
! moment, only the variable bioshade is returned in order to re-use
! it in the calling code for self-shading effects.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: bioshade(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   if (bioshade_feedback) then
      bioshade = bioshade_
   end if

   return
   end subroutine get_bio_updates
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory (Eulerian arrays)
!
! !INTERFACE:
   subroutine alloc_eul
!
! !DESCRIPTION:
!  Allocates memory for the arrays related to the Eulerian model. This
!  routine is only called once during the initialization procedure, 
!  with {\tt nmax} being the maximum number of layers that is expected
!  for the run. 
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard, Lars Umlauf, Karsten Bolding!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC

!  internal arrays
   allocate(par(0:nmax),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (par)'

   allocate(cc(1:numc,0:nmax),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (cc)'
   cc=_ZERO_

   allocate(ws(1:numc,0:nmax),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   ws=_ZERO_

   allocate(sfl(1:numc),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (sfl)'
   sfl=_ZERO_

   allocate(bfl(1:numc),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (bfl)'
   bfl=_ZERO_

   allocate(posconc(1:numc),stat=rc)
   if (rc /= 0) STOP 'allocate_memory(): Error allocating (posconc)'
   posconc=1

   allocate(zlev(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory_l: Error allocating (zlev)'
   zlev=_ZERO_

!  externally provided arrays
   allocate(h(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (h)'

   allocate(nuh(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (nuh)'

   allocate(t(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (t)'

   allocate(s(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (s)'

   allocate(rho(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (rho)'

   allocate(rad(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (rad)'

   allocate(w(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (w)'

   allocate(bioshade_(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (bioshade)'

   allocate(abioshade_(0:nmax),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating (abioshade)'


   return
   end subroutine alloc_eul
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory (particles)
!
! !INTERFACE:
   subroutine alloc_par
!
! !DESCRIPTION:
!  Allocates memory for the arrays related to the particle model.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf, Karsten Bolding
!EOP
!-----------------------------------------------------------------------!
! !LOCAL VARIABLES:
   integer                   :: rc

!-----------------------------------------------------------------------!
!BOC


   allocate(par_act(npar,ntype),stat=rc)
   if (rc /= 0) &
        STOP 'allocate_memory_par(): Error allocating (par_act)'

   allocate(par_ind(npar,ntype),stat=rc)
   if (rc /= 0) &
        STOP 'allocate_memory_par(): Error allocating (par_ind)'

   allocate(par_z(npar,ntype),stat=rc)
   if (rc /= 0) &
        STOP 'allocate_memory_par(): Error allocating (par_z)'

   allocate(par_prop(npar,nprop,ntype),stat=rc)
   if (rc /= 0) &
        STOP 'allocate_memory_par(): Error allocating (par_prop)'



   return
   end subroutine alloc_par
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Deallocate memory (Eulerian arrays)
!
! !INTERFACE:
   subroutine dealloc_eul()
!  Deallocates arrays related to the Eulerian model.
!
! !DESCRIPTION:
!  Deallocates arrays related to the Eulerian model.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Lars Umlauf, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

!  internal arrays
   if (allocated(par))            deallocate(par)
   if (allocated(cc))             deallocate(cc)
   if (allocated(ws))             deallocate(ws)
   if (allocated(sfl))            deallocate(sfl)
   if (allocated(bfl))            deallocate(bfl)
   if (allocated(posconc))        deallocate(posconc)
   if (allocated(var_ids))        deallocate(var_ids)
   if (allocated(var_names))      deallocate(var_names)
   if (allocated(var_units))      deallocate(var_units)
   if (allocated(var_long))       deallocate(var_long)
   if (allocated(zlev))           deallocate(zlev)

!  externally provided arrays
   if (allocated(h))              deallocate(h)
   if (allocated(nuh))            deallocate(nuh)
   if (allocated(t))              deallocate(t)
   if (allocated(s))              deallocate(s)
   if (allocated(rho))            deallocate(rho)
   if (allocated(rad))            deallocate(rad)
   if (allocated(w))              deallocate(w)
   if (allocated(bioshade_))      deallocate(bioshade_)
   if (allocated(abioshade_))     deallocate(abioshade_)

   return
   end subroutine dealloc_eul
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Deallocate memory (particles)
!
! !INTERFACE:
   subroutine dealloc_par
!
! !DESCRIPTION:
!  Deallocates arrays related to the particle model.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY: Lars Umlauf, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   if (allocated(par_act))        deallocate(par_act)
   if (allocated(par_ind))        deallocate(par_ind)
   if (allocated(par_z))          deallocate(par_z)
   if (allocated(par_prop))       deallocate(par_prop)


   return
   end subroutine dealloc_par
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine clean_bio
!
! !DESCRIPTION:
!  Deallocate memory and clean up some thins.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'clean_bio'

   ! clean up memory
   call dealloc_eul()

   if (.not. bio_eulerian) then

      call dealloc_par()

   end if

   ! clean up mussles
#if 0
   if (mussels_calc) then
      call end_mussels()
   end if
#endif

   
   LEVEL1 'done.'

   return
   end subroutine clean_bio
!EOC


   end module bio

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------

