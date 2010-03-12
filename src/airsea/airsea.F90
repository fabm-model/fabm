!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: airsea --- atmospheric fluxes \label{sec:airsea}
!
! !INTERFACE:
   module airsea
!
! !DESCRIPTION:
!  This module calculates the heat, momentum
!  and freshwater fluxes between the ocean and the atmosphere as well as
!  the incoming solar radiation. Fluxes and solar radiation may be
!  prescribed. Alternatively, they may be calculated by means
!  of bulk formulae from observed or modelled meteorological
!  parameters and the solar radiation may be calculated
!  from longitude, latitude,
!  time and cloudiness. For the prescibed fluxes and solar radiation,
!  values may be constant or read in from files. All necessary
!  setting have to be made in the namelist file {\tt airsea.nml}.
!
! !USES:
   use airsea_variables
   use time,         only: julian_day, time_diff
   use observations, only: read_obs
!
   IMPLICIT NONE

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public                              :: init_air_sea
   public                              :: do_air_sea
   public                              :: clean_air_sea
   public                              :: set_sst
   public                              :: integrated_fluxes
!
! !PUBLIC DATA MEMBERS:
   logical,  public                    :: calc_fluxes=.false.
!
!  wind speed (m/s)
   REALTYPE, public, target            :: w=_ZERO_
!
!  surface short-wave radiation
!  and surface heat flux (W/m^2)
   REALTYPE, public, target            :: I_0=_ZERO_
   REALTYPE, public                    :: heat=_ZERO_

!  surface stress components (Pa)
   REALTYPE, public                    :: tx=_ZERO_,ty=_ZERO_

!  precipitation and  evaporation
!  (m/s)
   REALTYPE, public, target            :: precip=_ZERO_
   REALTYPE, public, target            :: evap=_ZERO_

!  sea surface temperature (degC) and
!  sea surface salinity (psu)
   REALTYPE, public                    :: sst=_ZERO_
   REALTYPE, public                    :: sss=_ZERO_
   REALTYPE, public                    :: ssu=_ZERO_
   REALTYPE, public                    :: ssv=_ZERO_

!  integrated short-wave radiation,
!  surface heat flux (J/m^2)
   REALTYPE, public                    :: int_swr=_ZERO_,int_heat=_ZERO_

!  sum of short wave radiation
!  and surface heat flux (J/m^2)
   REALTYPE, public                    :: int_total=_ZERO_
!
! !DEFINED PARAMETERS:
   integer,  parameter                 :: meteo_unit=20
   integer,  parameter                 :: swr_unit=21
   integer,  parameter                 :: heatflux_unit=22
   integer,  parameter                 :: momentum_unit=23
   integer,  parameter                 :: precip_unit=24
   integer,  parameter                 :: sst_unit=25
   integer,  parameter                 :: sss_unit=26
   integer,  parameter                 :: CONSTVAL=1
   integer,  parameter                 :: FROMFILE=2
!
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Hans Burchard
!
!  $Log: airsea.F90,v $
!  Revision 1.30  2009-03-22 07:32:39  jorn
!  made swr clean-up independent of use of meteo data
!
!  Revision 1.29  2008-06-10 16:51:26  hb
!  bug concerning I_0 removed
!
!  Revision 1.28  2008-04-09 12:01:05  kb
!  initialise public variables to 0
!
!  Revision 1.27  2008-04-08 16:09:00  kb
!  assure valid qh and qe under all circumstance - Bruggeman, Stips
!
!  Revision 1.26  2008-01-02 15:11:48  kb
!  fixed w calculation - Kreuz
!
!  Revision 1.25  2008-01-02 14:36:37  kb
!  fixed indexing
!
!  Revision 1.24  2007-12-11 11:05:48  lars
!  corrected underscores in documentation
!
!  Revision 1.23  2007-12-09 10:25:48  kb
!  introduced swr_method: 1->constant, 2->from file, 3->equation
!
!  Revision 1.22  2007-12-07 10:12:20  kb
!  replaced p_e with precip and included evap
!
!  Revision 1.21  2007-09-25 10:06:10  kbk
!  modularized the airsea module - added Fairall method
!
!  Revision 1.20  2007-09-13 12:06:44  hb
!  fixed sign in momentum flux calculation
!
!  Revision 1.19  2007-06-26 18:24:29  jorn
!  made precipitation-evaporation clean-up independent of use of meteo data
!
!  Revision 1.18  2007-05-21 14:08:08  kbk
!  short wave radiation limitation suggested by Adolf Stips
!
!  Revision 1.17  2007-05-18 18:05:06  hb
!  Bug in short-wave radiation removed
!
!  Revision 1.16  2007-01-07 13:21:27  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.15  2006-12-08 06:50:37  kbk
!  fixed September in yday - Chris Locke
!
!  Revision 1.14  2006-11-27 10:08:33  kbk
!  use var init_saved_vars to initialise saved variables - air_sea_interaction -> do_air_sea
!
!  Revision 1.13  2006-11-17 07:13:17  kbk
!  rho amd wind-speed available via bio_var
!
!  Revision 1.12  2005/11/15 11:42:33  lars
!  documentation finish for print
!
!  Revision 1.11  2005/07/06 13:58:07  kbk
!  added fresh water, updated documentation
!
!  Revision 1.10  2004/07/30 09:19:03  hb
!  wet_mode now red from namelist
!
!  Revision 1.9  2004/06/25 07:50:29  hb
!  Preliminary wet mode choices improved
!
!  Revision 1.8  2004/05/28 13:14:14  hb
!  airsea.F90 extended for dew point temperature
!
!  Revision 1.7  2003/06/13 09:27:16  hb
!  Implemented freshwater fluxes
!
!  Revision 1.6  2003/03/28 09:20:34  kbk
!  added new copyright to files
!
!  Revision 1.5  2003/03/28 08:13:47  kbk
!  removed tabs
!
!  Revision 1.4  2003/03/10 08:37:56  gotm
!  HB fixed the Kondo calculations
!
!  Revision 1.3  2001/11/18 11:43:48  gotm
!  Cleaned
!
!  Revision 1.2  2001/06/13 07:40:39  gotm
!  Lon, lat was hardcoded in meteo.F90 - now passed via init_meteo()
!
!  Revision 1.1.1.1  2001/02/12 15:55:57  gotm
!  initial import into CVS
!
! !LOCAL VARIABLES:
   logical                   :: init_saved_vars=.true.
   integer                   :: swr_method=0
   integer                   :: fluxes_method=1
   integer                   :: back_radiation_method=1
   integer                   :: heat_method
   integer                   :: momentum_method
   integer                   :: precip_method=0
   integer                   :: sst_method=0
   integer                   :: sss_method=0
   integer                   :: hum_method
   logical, public           :: rain_impact=.false.
   logical, public           :: calc_evaporation=.false.

   character(len=PATH_MAX)   :: meteo_file
   character(len=PATH_MAX)   :: swr_file
   character(len=PATH_MAX)   :: heatflux_file
   character(len=PATH_MAX)   :: momentumflux_file
   character(len=PATH_MAX)   :: precip_file
   character(len=PATH_MAX)   :: sss_file
   character(len=PATH_MAX)   :: sst_file

   REALTYPE                  :: u10,v10
   REALTYPE                  :: airp
   REALTYPE                  :: airt,twet,tdew
   REALTYPE                  :: cloud
   REALTYPE                  :: rh

   REALTYPE                  :: const_swr=_ZERO_
   REALTYPE                  :: swr_factor=_ONE_
   REALTYPE                  :: const_heat
   REALTYPE                  :: const_tx,const_ty
   REALTYPE                  :: const_precip=_ZERO_
   REALTYPE                  :: precip_factor=_ONE_
   REALTYPE                  :: dlon,dlat
!
! !BUGS:
!  Wind speed - w - is not entirely correct. 
!  No temporal interpolation is done. If the momentum fluxes tx,ty are 
!  specified w=0. 
!  The Fairall and Kondo methods calculate their own w internally.
!  w is used by e.g. bio.F90
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the air--sea interaction module \label{sec:init-air-sea}
!
! !INTERFACE:
   subroutine init_air_sea(namlst,lat,lon)
!
! !DESCRIPTION:
!  This routine initialises the air-sea module by reading various variables
!  from the namelist {\tt airsea.nml} and opens relevant files.
!  These parameters are:
!
!  \begin{tabular}{ll}
! {\tt calc\_fluxes}     & {\tt .true.}: Sensible, latent and back-radiation are calculated by    \\
!                          & means of bulk formulae. In this case, {\tt meteo\_file} must be      \\
!                          & given and {\tt hum\_method} must be specified.                         \\
!                          & {\tt .false.}: Surface fluxes and solar radiation are prescribed.      \\
! {\tt fluxes\_method}   & Select which parameterisation to use for latent and sensible fluxes:   \\
!                          & 1: Kondo (1975)                                                        \\
!                          & 2: Fairall et al. (1996)                                               \\
! {\tt back\_radiation\_method}   & Select which parameterisation to use:                        \\
!                          & 1: Clark et al. (1974)                                                 \\
!                          & 2: Hastenrath and Lamb (1978)                                          \\
!                          & 3: Bignami et al. (1995)                                               \\
!                          & 4: Berliandand Berliand (1952)                                         \\
! {\tt meteo\_file}      & file with meteo data (for {\tt calc\_fluxes=.true.}) with            \\
!                          & date: {\tt yyyy-mm-dd hh:mm:ss}                                        \\
!                          & $x$-component of wind (10 m) in m\,s$^{-1}$                            \\
!                          & $y$-component of wind (10 m) in m\,s$^{-1}$                            \\
!                          & air pressure (2 m) in hectopascal                                      \\
!                          & dry air temperature (2 m) in Celsius                                   \\
!                          & rel. hum. in \% or wet bulb temp. in C or dew point temp. in C         \\
!                          & cloud cover in 1/10                                                    \\
!                          & Example:                                                               \\
!                          & {\tt 1998-01-01 00:00:00    6.87  10.95 1013.0   6.80   73.2   0.91}   \\
! {\tt hum\_method}        & 1: relative humidity given as 7.\ column in {\tt meteo\_file}        \\
!                          & 2: wet bulb temperature given as 7. column in {\tt meteo\_file}      \\
!                          & 3: dew point temperature given as 7. column in {\tt meteo\_file}     \\
! {\tt heat\_method}     & 0: heat flux not prescribed                                            \\
!                          & 1: constant value for short wave radiation ({\tt const\_swr})        \\
!                          &    and surface heat flux ({\tt const\_qh})                           \\
!                          & 2: {\tt swr}, {\tt heat} are read from {\tt heatflux\_file}          \\
! {\tt rain\_impact}     & {\tt .true.}: include sensible- and momentum-flux from precipitation   \\
! {\tt calc\_evaporation}& {\tt .true.}: calculate evaporation/condensation (m/s)                 \\
! {\tt swr\_method}      & 1: constant value for short wave radiation ({\tt const\_swr})        \\
!                          & 2: read short wave radiation from file                                 \\
!                          & 3: Solar radiation is calculated from time, longitude, latitude,       \\
!                          & and cloud cover.                                                       \\
! {\tt const\_swr}       & constant value for short wave radiation in W\,m$^{-2}$                 \\
!                          & (always positive)                                                      \\
! {\tt swr\_file}        & file with short wave radiation in W\,m$^{-2}$                          \\
! {\tt swr\_factor}      & scales data read from file to  W\,m$^{-2}$ - defaults to 1             \\
! {\tt const\_heat }     & constant value for surface heat flux in  W\,m$^{-2}$                   \\
!                          & (negative for heat loss)                                               \\
! {\tt heatflux\_file}   & file with date and {\tt heat} in W\,m$^{-2}$                           \\
! {\tt momentum\_method} & 0: momentum flux not prescribed                                        \\
!                          & 1: constant momentum fluxes {\tt const\_tx}, {\tt const\_tx} given \\
!                          & 2: surface momentum fluxes given from file {\tt momentumflux\_file}  \\
! {\tt const\_tx}        & $x$-component of constant surface momentum flux in N\,m$^{-2}$         \\
! {\tt const\_ty}        & $y$-component of constant surface momentum flux in N\,m$^{-2}$         \\
! {\tt momentumflux\_file} & File with date, {\tt tx} and {\tt ty} given                          \\
! {\tt precip\_method}   & 0: precipitation not included == precip=0.                             \\
!                          & 1: constant value for precipitation in  m\,s$^{-1}$                    \\
!                          & 2: values for precipitation read from file {\tt precip\_file}         \\
! {\tt const\_precip}    & value for precip in m\,s$^{-1}$                                        \\
! {\tt precip\_file}     & file with date and {\tt precip} in m\,s$^{-1}$                         \\
! {\tt precip\_factor}   & scales data read from file to  m\,s$^{-1}$ - defaults to 1             \\
! {\tt sst\_method}      & 0: no independent SST observation is read from file                    \\
!                          & 2: independent SST observation is read from file, only for output      \\
! {\tt sst\_file}        & file with date and SST (sea surface temperature) in Celsius            \\
! {\tt sss\_method}      & 0: no independent SSS observation is read from file                    \\
!                          & 2: independent SSS observation is read from file, only for output      \\
! {\tt sss\_file}        & file with date and SSS (sea surface salinity) in psu                   \\
!  \end{tabular}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
   REALTYPE, intent(in)                :: lat,lon
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   namelist /airsea/ calc_fluxes, &
                     fluxes_method, &
                     back_radiation_method, &
                     meteo_file, &
                     hum_method, &
                     heat_method, &
                     rain_impact, &
                     calc_evaporation, &
                     swr_method,const_swr,swr_file,swr_factor, &
                     const_heat, &
                     heatflux_file, &
                     momentum_method, &
                     const_tx,const_ty, &
                     momentumflux_file, &
                     precip_method,const_precip,precip_file,precip_factor,&
                     sst_method, sst_file, &
                     sss_method, sss_file
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_air_sea'

   open(namlst,file='airsea.nml',action='read',status='old',err=90)
   read(namlst,nml=airsea,err=91)
   close(namlst)

!  The short wave radiation
   select case (swr_method)
      case (1)
         LEVEL2 'Using constant short wave radiation= ',const_swr
      case (FROMFILE)
         open(swr_unit,file=swr_file,action='read', &
              status='old',err=92)
         LEVEL2 'Reading short wave radiation data from:'
         LEVEL3 trim(swr_file)
         if (swr_factor .ne. _ONE_) then
            LEVEL3 'applying swr factor= ',swr_factor
         end if
      case (3)
         if (.not. calc_fluxes) then
            LEVEL2 'Not possible to calculate swr - when calc_fluxes=.false.'
            stop 'init_airsea'
         else
            LEVEL2 'Calculating swr=swr(t(lon),lat,cloud)'
         end if
      case default
   end select

   if (calc_fluxes) then

      open(meteo_unit,file=meteo_file,action='read',status='old',err=93)
      LEVEL2 'Air-sea exchanges will be calculated'
      LEVEL2 'Reading meteo data from:'
      LEVEL3 trim(meteo_file)
      LEVEL3 'heat- and momentum-fluxes:'
      select case (fluxes_method)
         case(1)
            LEVEL4 'using Kondo formulation'
         case(2)
            LEVEL4 'using Fairall et. all formulation'
         case default
      end select
      LEVEL3 'long-wave back radiation:'
      select case (back_radiation_method)
         case(1)
            LEVEL4 'using Clark formulation'
         case(2)
            LEVEL4 'using Hastenrath formulation'
         case(3)
            LEVEL4 'using Bignami formulation'
         case(4)
            LEVEL4 'using Berliand formulation'
         case default
      end select

   else

!     The heat fluxes
      select case (heat_method)
         case (FROMFILE)
            open(heatflux_unit,file=heatflux_file,action='read', &
                 status='old',err=94)
            LEVEL2 'Reading heat fluxes from:'
            LEVEL3 trim(heatflux_file)
         case default
      end select

!     The momentum fluxes
      select case (momentum_method)
         case (FROMFILE)
            open(momentum_unit,file=momentumflux_file,action='read', &
                 status='old',err=95)
            LEVEL2 'Reading momentum fluxes from:'
            LEVEL3 trim(momentumflux_file)
         case default
      end select

!     The sea surface temperature
      select case (sst_method)
         case (FROMFILE)
            open(sst_unit,file=sst_file,action='read',status='old',err=96)
            LEVEL2 'Reading sea surface temperature from:'
            LEVEL3 trim(sst_file)
         case default
      end select

!     The sea surface salinity
      select case (sss_method)
         case (FROMFILE)
            open(sss_unit,file=sss_file,action='read',status='old',err=97)
            LEVEL2 'Reading sea surface salinity from:'
            LEVEL3 trim(sss_file)
         case default
      end select

   end if

!  The fresh water fluxes (used for calc_fluxes=.true. and calc_fluxes=.false.)
   select case (precip_method)
      case (1)
         LEVEL2 'Using constant precipitation= ',const_precip
         LEVEL2 'rain_impact=      ',rain_impact
         LEVEL2 'calc_evaporation= ',calc_evaporation
      case (FROMFILE)
         open(precip_unit,file=precip_file,action='read', &
              status='old',err=98)
         LEVEL2 'Reading precipitation data from:'
         LEVEL3 trim(precip_file)
         if (precip_factor .ne. _ONE_) then
            LEVEL3 'applying factor= ',precip_factor
         end if
         LEVEL2 'rain_impact=      ',rain_impact
         LEVEL2 'calc_evaporation= ',calc_evaporation
      case default
   end select

   twet=0.
   tdew=0.
   rh=0.
   cloud=0.
   sss=0.
   airt=0.

   dlon = lon
   dlat = lat

   return

90 FATAL 'I could not open airsea.nml'
   stop 'init_airsea'
91 FATAL 'I could not read airsea namelist'
   stop 'init_airsea'
92 FATAL 'I could not open ',trim(swr_file)
   stop 'init_airsea'
93 FATAL 'I could not open ',trim(meteo_file)
   stop 'init_airsea'
94 FATAL 'I could not open ',trim(heatflux_file)
   stop 'init_airsea'
95 FATAL 'I could not open ',trim(momentumflux_file)
   stop 'init_airsea'
96 FATAL 'I could not open ',trim(sst_file)
   stop 'init_airsea'
97 FATAL 'I could not open ',trim(sss_file)
   stop 'init_airsea'
98 FATAL 'I could not open ',trim(precip_file)
   stop 'init_airsea'

   end subroutine init_air_sea
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the air--sea fluxes
!
! !INTERFACE:
   subroutine do_air_sea(jul,secs)
!
! !DESCRIPTION:
!
!  Depending on the value of the boolean variable {\tt calc\_fluxes},
!  the subroutines for the calculation of the fluxes
!  and the short wave radiation are
!  called or the fluxes are directly read in from the namelist
!  {\tt airsea.nml} as constants or read in from files.
!  Furthermore, the surface freshwater flux is set to a constant
!  value or is read in from a file.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!-----------------------------------------------------------------------
!BOC

!  The short wave radiation
   select case (swr_method)
      case (CONSTVAL)
         I_0=const_swr
      case (FROMFILE)
         call read_swr(jul,secs,I_0)
      case default
   end select

!  The freshwater flux (used for calc_fluxes=.true. and calc_fluxes=.false.)
   select case (precip_method)
      case (CONSTVAL)
         precip=const_precip
      case (FROMFILE)
         call read_precip(jul,secs,precip)
      case default
   end select

   if (calc_fluxes) then
      call flux_from_meteo(jul,secs)
      if (swr_method .eq. 3) then
         call short_wave_radiation(jul,secs,dlon,dlat,cloud,I_0)
      end if
   else
!     The heat fluxes
      select case (heat_method)
         case (CONSTVAL)
            heat=const_heat
         case (FROMFILE)
            call read_heat_flux(jul,secs,heat)
         case default
      end select
!     The momentum fluxes
      select case (momentum_method)
         case (CONSTVAL)
            tx=const_tx
            ty=const_ty
         case (FROMFILE)
            call read_momentum_flux(jul,secs,tx,ty)
         case default
      end select
!     The sea surface temperature
      select case (sst_method)
         case (FROMFILE)
            call read_sst(jul,secs,sst)
         case default
      end select
!     The sea surface salinity
      select case (sss_method)
         case (FROMFILE)
            call read_sss(jul,secs,sss)
         case default
      end select
   end if

   if (init_saved_vars) init_saved_vars=.false.

   return
   end subroutine do_air_sea
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the air--sea interactions
!
! !INTERFACE:
   subroutine clean_air_sea
!
! !DESCRIPTION:
!  All files related to air-sea interaction which have been opened
!  are now closed by this routine.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (calc_fluxes) then
      close(meteo_unit)
   else
      if (heat_method     .eq. FROMFILE) close(heatflux_unit)
      if (momentum_method .eq. FROMFILE) close(momentum_unit)
      if (precip_method   .eq. FROMFILE) close(precip_unit)
      if (sst_method      .eq. FROMFILE) close(sst_unit)
      if (sss_method      .eq. FROMFILE) close(sss_unit)
   end if
   if (swr_method    .eq. FROMFILE) close(swr_unit)
   if (precip_method .eq. FROMFILE) close(precip_unit)
   init_saved_vars=.true.
   return
   end subroutine clean_air_sea
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read meteo data, interpolate in time
!
! !INTERFACE:
   subroutine flux_from_meteo(jul,secs)
!
! !DESCRIPTION:
!  For {\tt calc\_fluxes=.true.}, this routine reads meteo data
!  from {\tt meteo\_file} and calculates for each time step
!  in  {\tt meteo\_file} the
!  fluxes of heat and momentum, and the
!  long-wave back radiation by calling the routines {\tt humidity},
!  {\tt back\_radiation} and {\tt airsea\_fluxes}, see sections
!  \sect{sec:humidity}, \sect{sec:back-rad}, and \sect{sec:airsea-fluxes},
!  a wrapper routine for using the bulk fomulae from either \cite{Kondo75}
!  or \cite{Fairalletal96}. Afterwards, the airsea fluxes
!  are interpolated to the actual time step of GOTM. Finally, the
!  incoming short-wave radiation is calculated by using the interpolated
!  cloud cover and the actual UTC time of GOTM, see the routine
!  {\tt short\_wave\_radiation} described in \sect{sec:swr}.

!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t
   REALTYPE, SAVE            :: dt
   integer, save             :: meteo_jul1,meteo_secs1
   integer, save             :: meteo_jul2=0,meteo_secs2=0
   REALTYPE, save            :: obs(7)
   REALTYPE, save            :: alpha(4)
   REALTYPE, save            :: h1,tx1,ty1
   REALTYPE, save            :: h2=0.,tx2=0.,ty2=0.
   logical, save             :: first=.true.
   integer                   :: rc
   REALTYPE                  :: ta,ta_k,tw,tw_k 
   REALTYPE                  :: qe=_ZERO_,qh=_ZERO_,qb=_ZERO_
!
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      meteo_jul2=0
      meteo_secs2=0
      h2=0.
      tx2=0.
      ty2=0.
      first=.true.
   end if
!  This part initialises and reads in new values if necessary.
   if(time_diff(meteo_jul2,meteo_secs2,jul,secs) .lt. 0) then
      do
         meteo_jul1 = meteo_jul2
         meteo_secs1 = meteo_secs2
         call read_obs(meteo_unit,yy,mm,dd,hh,min,ss,6,obs,rc)
         call julian_day(yy,mm,dd,meteo_jul2)
         meteo_secs2 = hh*3600 + min*60 + ss
         if(time_diff(meteo_jul2,meteo_secs2,jul,secs) .gt. 0) EXIT
      end do
      u10    = obs(1)
      v10    = obs(2)
      airp  = obs(3)*100. !kbk mbar/hPa --> Pa
      airt  = obs(4)
      rh    = obs(5)
      cloud = obs(6)

      if (sst .lt. 100.) then
         tw  = sst
         tw_k= sst+KELVIN
      else
         tw  = sst-KELVIN
         tw_k= sst
      end if

      if (airt .lt. 100.) then
         ta_k  = airt + KELVIN
         ta = airt
      else
         ta  = airt - KELVIN
         ta_k = airt
      end if

      if (first) then
         call humidity(hum_method,rh,airp,tw,ta)
         call back_radiation(back_radiation_method, &
                             dlat,tw_k,ta_k,cloud,qb)
         call airsea_fluxes(fluxes_method, &
                            tw,ta,u10-ssu,v10-ssv,precip,evap,tx1,ty1,qe,qh)
         h1=(qb+qe+qh)
         h2  = h1
         tx2 = tx1
         ty2 = ty1
         first = .false.
      else
         h1  = h2
         tx1 = tx2
         ty1 = ty2
         call humidity(hum_method,rh,airp,tw,ta)
         call back_radiation(back_radiation_method, &
                             dlat,tw_k,ta_k,cloud,qb)
         call airsea_fluxes(fluxes_method, &
                            tw,ta,u10-ssu,v10-ssv,precip,evap,tx2,ty2,qe,qh)
         h2=(qb+qe+qh)
      end if
      dt = time_diff(meteo_jul2,meteo_secs2,meteo_jul1,meteo_secs1)
      alpha(2) = (h2-h1)/dt
      alpha(3) = (tx2-tx1)/dt
      alpha(4) = (ty2-ty1)/dt
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,meteo_jul1,meteo_secs1)
   heat = h1  + t*alpha(2)
   tx   = tx1 + t*alpha(3)
   ty   = ty1 + t*alpha(4)

   w = sqrt(u10*u10+v10*v10)

   return
   end subroutine flux_from_meteo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read short wave radiation, interpolate in time
!
! !INTERFACE:
   subroutine read_swr(jul,secs,swr)
!
! !DESCRIPTION:
!  This routine reads the short wave radiation (in W\,s$^{-2}$) from
!  {\tt swr\_file} and interpolates in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: swr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: swr_jul1,swr_secs1
   integer, save             :: swr_jul2=0,swr_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      swr_jul2=0
      swr_secs2=0
      obs2(1)=0.
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(swr_jul2,swr_secs2,jul,secs) .lt. 0) then
      do
         swr_jul1 = swr_jul2
         swr_secs1 = swr_secs2
         obs1 = obs2
         call read_obs(swr_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,swr_jul2)
         swr_secs2 = hh*3600 + min*60 + ss
         if(time_diff(swr_jul2,swr_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(swr_jul2,swr_secs2,swr_jul1,swr_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,swr_jul1,swr_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   swr = (obs1(1) + t*alpha)*swr_factor

   return
   end subroutine read_swr

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read heat flux data, interpolate in time
!
! !INTERFACE:
   subroutine read_heat_flux(jul,secs,heat)
!
! !DESCRIPTION:
!   For {\tt calc\_fluxes=.false.}, this routine reads in
!   the sum of sensible, latent and long-wave back-radiation flux 
!   in W\,m$^{-2}$ from {\tt heatflux\_file} and interpolate in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: heat
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, SAVE            :: dt
   integer, save             :: heat_jul1,heat_secs1
   integer, save             :: heat_jul2=0,heat_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      heat_jul2=0
      heat_secs2=0
      obs2(1)=0.
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(heat_jul2,heat_secs2,jul,secs) .lt. 0) then
      do
         heat_jul1 = heat_jul2
         heat_secs1 = heat_secs2
         obs1 = obs2
         call read_obs(heatflux_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,heat_jul2)
         heat_secs2 = hh*3600 + min*60 + ss
         if(time_diff(heat_jul2,heat_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(heat_jul2,heat_secs2,heat_jul1,heat_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,heat_jul1,heat_secs1)

   alpha = (obs2(1)-obs1(1))/dt
   heat = obs1(1) + t*alpha

   return
   end subroutine read_heat_flux
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read momentum flux data, interpolate in time
!
! !INTERFACE:
   subroutine read_momentum_flux(jul,secs,tx,ty)
!
! !DESCRIPTION:
!   For {\tt calc\_fluxes=.false.}, this routine reads momentum fluxes
!   in N\,m$^{-2}$ from \linebreak
!   {\tt momentumflux\_file} and interpolates them in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                  :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: tx,ty
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: mom_jul1,mom_secs1
   integer, save             :: mom_jul2=0,mom_secs2=0
   REALTYPE, save            :: obs1(2),obs2(2)=0.
   integer                   :: rc
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      mom_jul2=0
      mom_secs2=0
      obs2(2)=0.
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(mom_jul2,mom_secs2,jul,secs) .lt. 0) then
      do
         mom_jul1 = mom_jul2
         mom_secs1 = mom_secs2
         obs1 = obs2
         call read_obs(momentum_unit,yy,mm,dd,hh,min,ss,2,obs2,rc)
         call julian_day(yy,mm,dd,mom_jul2)
         mom_secs2 = hh*3600 + min*60 + ss
         if(time_diff(mom_jul2,mom_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(mom_jul2,mom_secs2,mom_jul1,mom_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,mom_jul1,mom_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   tx = obs1(1) + t*alpha
   alpha = (obs2(2)-obs1(2))/dt
   ty = obs1(2) + t*alpha

   return
   end subroutine read_momentum_flux
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read precipitation, interpolate in time
!
! !INTERFACE:
   subroutine read_precip(jul,secs,precip)
!
! !DESCRIPTION:
!  This routine reads the precipitation (in m\,s$^{-1}$) from
!  {\tt precip\_file} and interpolates in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: precip
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: precip_jul1,precip_secs1
   integer, save             :: precip_jul2=0,precip_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      precip_jul2=0
      precip_secs2=0
      obs2(1)=0.
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(precip_jul2,precip_secs2,jul,secs) .lt. 0) then
      do
         precip_jul1 = precip_jul2
         precip_secs1 = precip_secs2
         obs1 = obs2
         call read_obs(precip_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,precip_jul2)
         precip_secs2 = hh*3600 + min*60 + ss
         if(time_diff(precip_jul2,precip_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(precip_jul2,precip_secs2,precip_jul1,precip_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,precip_jul1,precip_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   precip = (obs1(1) + t*alpha)*precip_factor

   return
   end subroutine read_precip

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read SST, interpolate in time
!
! !INTERFACE:
   subroutine read_sst(jul,secs,sst)
!
! !DESCRIPTION:
!   For {\tt calc\_fluxes=.false.}, this routine reads sea surface
!   temperature (SST) from {\tt sst\_file}
!  and interpolates in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: sst
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: sst_jul1,sst_secs1
   integer, save             :: sst_jul2=0,sst_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
     sst_jul2=0
     sst_secs2=0
     obs2(1)=0.
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(sst_jul2,sst_secs2,jul,secs) .lt. 0) then
      do
         sst_jul1 = sst_jul2
         sst_secs1 = sst_secs2
         obs1 = obs2
         call read_obs(sst_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,sst_jul2)
         sst_secs2 = hh*3600 + min*60 + ss
         if(time_diff(sst_jul2,sst_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(sst_jul2,sst_secs2,sst_jul1,sst_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,sst_jul1,sst_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   sst = obs1(1) + t*alpha

   return
   end subroutine read_sst
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read SSS, interpolate in time
!
! !INTERFACE:
   subroutine read_sss(jul,secs,sss)
!
! !DESCRIPTION:
!   For {\tt calc\_fluxes=.false.}, this routine reads sea surface
!   salinity (SSS) from {\tt sss\_file}
!  and interpolates in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   REALTYPE,intent(out)                :: sss
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: sss_jul1,sss_secs1
   integer, save             :: sss_jul2=0,sss_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      sss_jul2=0
      sss_secs2=0
      obs2(1)=0.
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(sss_jul2,sss_secs2,jul,secs) .lt. 0) then
      do
         sss_jul1 = sss_jul2
         sss_secs1 = sss_secs2
         obs1 = obs2
         call read_obs(sss_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,sss_jul2)
         sss_secs2 = hh*3600 + min*60 + ss
         if(time_diff(sss_jul2,sss_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(sss_jul2,sss_secs2,sss_jul1,sss_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(jul,secs,sss_jul1,sss_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   sss = obs1(1) + t*alpha

   return
   end subroutine read_sss
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Integrate short--wave and sea surface fluxes
!
! !INTERFACE:
   subroutine integrated_fluxes(dt)
!
! !DESCRIPTION:
!  This utility routine integrates the short--wave radiation
!  and heat--fluxes over time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: dt
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   int_swr   = int_swr + dt*I_0
   int_heat  = int_heat + dt*heat
   int_total = int_swr + int_heat
   return
   end subroutine integrated_fluxes
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set the SST to be used from model.
!
! !INTERFACE:
   subroutine set_sst(temp,u,v)
!
! !DESCRIPTION:
!  This routine sets the simulated
!  sea surface temperature (SST) to be used for
!  the surface flux calculations.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: temp,u,v
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  See log for airsea module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   sst = temp
   !ssu = -u
   !ssv = -v
   return
   end subroutine set_sst
!EOC

!-----------------------------------------------------------------------

   end module airsea

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
