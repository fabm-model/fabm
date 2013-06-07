#ifdef _FABM_F2003_
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_au_sediment --- Benthic biogeochemical model
! Aarhus university, department of BioScience
!
! !INTERFACE:
   module fabm_au_sediment
!
! !DESCRIPTION:
!  The benthic biogeochemical model .....
!
!
! !USES:
   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Karen Timmermann
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model), public :: type_au_sediment
!     Variable identifiers
      type (type_state_variable_id)        :: id_C, id_W, id_M !Pelagic Carbon, phosphor and Nitrogen (if model is coupled to pelagic)
      type (type_bottom_state_variable_id) :: id_sediC, id_sediW, id_sediM !Carbon, phosphor and Nitrogen pools
      type (type_bottom_state_variable_id) :: id_sediPO4, id_mePO4 !phosphor state variables
      type (type_bottom_state_variable_id) :: id_sediNH4, id_sediNO3, id_sediN2 ! Nitrogen state variables
      type (type_state_variable_id)        :: id_PO4, id_NH4, id_NO3 !pelagic bottom varaiables which this model will modify
      type (type_state_variable_id)        :: id_O2
      type (type_dependency_id)            :: id_temp

      type (type_conserved_quantity_id)    :: id_totN, id_totP

!     Model parameters
      real(rk) :: sorpRconst, DPO4, DNH4, DNO3, kDenit, sediNHrmax, KO2NH, KO2C, diffLength
      real(rk) :: Tref, qT, minQWC, Wrmax, minQMC, Mrmax, Crmax, sediTopDZ, w_d !_KB_, dt
      logical  :: couple2Pelagic

      contains

      procedure :: do_benthos
      procedure :: get_conserved_quantities
   end type type_au_sediment
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the benthic model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the au namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_au_sediment), intent(inout), target :: self
   integer,                  intent(in)            :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   real(rk) :: sediC_initial
   real(rk) :: sediW_initial
   real(rk) :: sediM_initial
   real(rk) :: sediPO4_initial
   real(rk) :: mePO4_initial
   real(rk) :: sediNH4_initial
   real(rk) :: sediNO3_initial
   real(rk) :: sediN2_initial

   real(rk) :: sorpRconst
   real(rk) :: DPO4
   real(rk) :: DNH4
   real(rk) :: DNO3
   real(rk) :: kDenit
   real(rk) :: sediNHrmax !sec-1, max nitrification rate
   real(rk) :: KO2NH      !mmol O2/m-3, half saturation for O2 used in nitrification
   real(rk) :: KO2C
   real(rk) :: diffLength
   real(rk) :: Tref
   real(rk) :: qT
   real(rk) :: minQWC
   real(rk) :: Wrmax
   real(rk) :: minQMC
   real(rk) :: Mrmax
   real(rk) :: Crmax
   real(rk) :: sediTopDZ
   real(rk) :: w_d       !detritus settling rate (m/s)
!_KB_
!_KB_   real(rk) :: dt
!_KB_

   character(len=64) :: C_var
   character(len=64) :: M_var
   character(len=64) :: W_var
   character(len=64) :: PO4_var
   character(len=64) :: NH4_var
   character(len=64) :: NO3_var
   character(len=64) :: O2_var
   logical  :: couple2Pelagic
   real(rk) :: PO4_initial
   real(rk) :: NH4_initial
   real(rk) :: NO3_initial
   real(rk) :: O2_initial

   namelist /au_sediment/ sediC_initial, sediW_initial, sediM_initial, sediPO4_initial, mePO4_initial, sediNH4_initial, sediNO3_initial, sediN2_initial, &
                            sorpRconst, DPO4, DNH4, DNO3, kDenit, sediNHrmax, KO2NH, KO2C, diffLength, Tref, qT, minQWC, Wrmax, minQMC, Mrmax, Crmax, &
                            C_var, M_var, W_var, PO4_var, NH4_var, NO3_var, O2_var, couple2Pelagic, PO4_initial, NH4_initial, NO3_initial, O2_initial, &
                            sediTopDZ, w_d !, _KB_ dt
!EOP
!-----------------------------------------------------------------------
!BOC
   sediC_initial=1_rk
   sediW_initial=1_rk
   sediM_initial=1_rk
   sediPO4_initial=1_rk
   mePO4_initial=1_rk
   sediNH4_initial=1_rk
   sediNO3_initial=1_rk
   sediN2_initial=1_rk

   sorpRconst=1e-8_rk
   DPO4=1e-8_rk
   DNH4=1e-5_rk
   DNO3=1e-5_rk
   kDenit=1.16e-7_rk
   sediNHrmax=1.16e-7_rk
   KO2NH=30_rk
   KO2C=10_rk
   diffLength=0.1_rk
   Tref=20_rk
   qT=0.07_rk
   minQWC=0.005_rk
   Wrmax=9.26e-7_rk
   minQMC=0.09_rk
   Mrmax=9.26e-7_rk
   Crmax=6.94e-7_rk
   sediTopDZ=0.1_rk
   w_d=0.0001_rk
!_KB_   dt=60_rk

   C_var='au_daneco_C'
   M_var='au_daneco_M'
   W_var='au_daneco_W'
   PO4_var='au_daneco_PO4'
   NH4_var='au_daneco_NH4'
   NO3_var='au_daneco_NO3'
   O2_var='au_daneco_O2'
   couple2Pelagic=.false.
   PO4_initial=1
   NH4_initial=1
   NO3_initial=1
   O2_initial=1

!  Read the namelist
   read(configunit,nml=au_sediment,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%sorpRconst=sorpRconst
   self%DPO4=DPO4
   self%DNH4=DNH4
   self%DNO3=DNO3
   self%kDenit=kDenit
   self%sediNHrmax=sediNHrmax
   self%KO2NH=KO2NH
   self%KO2C=KO2C
   self%diffLength=diffLength
   self%Tref=Tref
   self%qT=qT
   self%minQWC=minQWC
   self%Wrmax=Wrmax
   self%minQMC=minQMC
   self%Mrmax=Mrmax
   self%Crmax=Crmax
   self%sediTopDZ=sediTopDZ
   self%w_d = w_d
!_KB_   self%dt = dt
   self%couple2Pelagic = couple2Pelagic

!  Register state variables

   call self%register_state_variable(self%id_sediC,'sediC','mmol-C/m**2','benthic carbon',sediC_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sediW,'sediW','mmol/m**2','benthic phospor',sediW_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sediM,'sediM','mmol/m**2','benthic nitrogen',sediM_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sediPO4,'sediPO4','mmol/m**2','benthic PO4',sediPO4_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_mePO4,'mePO4','mmol/m**2','benthic metal PO4',mePO4_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sediNH4,'sediNH4','mmol/m**2','benthic NH4',sediNH4_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sediNO3,'sediNO3','mmol/m**2','benthic NO3',sediNO3_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sediN2,'sediN2','mmol/m**2','outgassed N2',sediN2_initial,minimum=_ZERO_)

!  Register environmental dependencies
   call self%register_dependency(self%id_temp,varname_temp)
   if (couple2Pelagic) then
      call self%register_state_dependency(self%id_C,C_var)
      call self%register_state_dependency(self%id_W,M_var)
      call self%register_state_dependency(self%id_M,W_var)
      call self%register_state_dependency(self%id_PO4,PO4_var)
      call self%register_state_dependency(self%id_NH4,NH4_var)
      call self%register_state_dependency(self%id_NO3,NO3_var)
      call self%register_state_dependency(self%id_O2,O2_var)
   else
      call self%register_state_variable(self%id_PO4,PO4_var,'mmol/m**2','pelagic PO4',PO4_initial,minimum=_ZERO_,no_river_dilution=.true.)
      call self%register_state_variable(self%id_NH4,NH4_var,'mmol/m**2','pelagic NH4',NH4_initial,minimum=_ZERO_,no_river_dilution=.true.)
      call self%register_state_variable(self%id_NO3,NO3_var,'mmol/m**2','pelagic NO3',NO3_initial,minimum=_ZERO_,no_river_dilution=.true.)
      call self%register_state_variable(self%id_O2,O2_var,'mmol/m**2','pelagic O2',O2_initial,minimum=_ZERO_,no_river_dilution=.true.)
   end if

!  Register conserved quantities
   call self%register_conserved_quantity(self%id_totN,'N','mmol/m**3','totN')
   call self%register_conserved_quantity(self%id_totP,'P','mmol/m**3','totP')

   return

99 call fatal_error('initialize','Error reading namelist au_sediment')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   subroutine get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   class (type_au_sediment), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk) :: sediNsum, sediPsum
   real(rk) :: botM, botW !pelagic bottom concentrations
   real(rk) :: sediW, sediM, sediPO4, mePO4, sediNH4, sediNO3, sediN2 ! local benthic state variables
   real(rk) :: botPO4, botNH4, botNO3 ! local pelagic state variables
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   ! Benthic
   _GET_HORIZONTAL_(self%id_sediW, sediW)!
   _GET_HORIZONTAL_(self%id_sediM, sediM)!
   _GET_HORIZONTAL_(self%id_sediPO4, sediPO4)!
   _GET_HORIZONTAL_(self%id_mePO4, mePO4)!
   _GET_HORIZONTAL_(self%id_sediNH4, sediNH4)!
   _GET_HORIZONTAL_(self%id_sediNO3, sediNO3)!
   _GET_HORIZONTAL_(self%id_sediN2, sediN2)!

   sediNsum = sediM+sediNH4+sediNO3+sediN2
   sediPsum = sediW+sediPO4+mePO4

   if (.not.self%couple2Pelagic) then
      ! Pelagic
      _GET_(self%id_PO4, botPO4)!
      _GET_(self%id_NH4, botNH4)!
      _GET_(self%id_NO3, botNO3)!
      sediNsum = sediNsum + botNH4 +botNO3
      sediPsum = sediPsum + botPO4
   end if

   _SET_CONSERVED_QUANTITY_(self%id_totN,sediNsum)
   _SET_CONSERVED_QUANTITY_(self%id_totP,sediPsum)

   ! Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_

   end subroutine get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!
! !DESCRIPTION:
! This routine calculates the benthic sink and source terms, as well as
! (matching) bottom fluxes for pelagic variables. Both have units mmol/m**2/s.
!
! !INPUT PARAMETERS:
   class (type_au_sediment),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Karen Timmermann & Janus Larsen
!
! !LOCAL VARIABLES:
   real(rk) :: botC, botM, botW !pelagic bottom concentrations
   real(rk) :: sediC, sediW, sediM, sediPO4, mePO4, sediNH4, sediNO3, sediN2 ! local benthic state variables
   real(rk) :: botPO4, botNH4, botNO3, botO2, botTemp ! local pelagic state variables
   real(rk) :: sediFT !common
   real(rk) :: sediWr, sediQWC, sediRsorp, sediRdesorp, JPO4 !phosphor
   real(rk) :: sediQMC, sediMr, JNH4, JNO3, sediNHr  !nitrogen
   real(rk) :: sediCr ! Carbon
   real(rk) :: botDZ

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   ! Benthic
   _GET_HORIZONTAL_(self%id_sediC, sediC)
   _GET_HORIZONTAL_(self%id_sediW, sediW)
   _GET_HORIZONTAL_(self%id_sediM, sediM)
   _GET_HORIZONTAL_(self%id_sediPO4, sediPO4)
   _GET_HORIZONTAL_(self%id_mePO4, mePO4)
   _GET_HORIZONTAL_(self%id_sediNH4, sediNH4)
   _GET_HORIZONTAL_(self%id_sediNO3, sediNO3)
   _GET_HORIZONTAL_(self%id_sediN2, sediN2)
   ! Pelagic
   _GET_(self%id_PO4, botPO4)
   _GET_(self%id_NH4, botNH4)
   _GET_(self%id_NO3, botNO3)
   _GET_(self%id_O2, botO2)
   _GET_(self%id_temp, botTemp)

   !calc sedimentation from pelagic bottom layer to sedi pools
   if (self%couple2Pelagic) then
      _GET_(self%id_C, botC)
      _GET_(self%id_M, botM)
      _GET_(self%id_W, botW)
      !mass that falls out of bottom pelagic layer = fallRate[m/s] * CellArea[m2] * timestep[s] * concentration [kg/m3]
      !new concentration in sediment = [old conc] + [new mass]/[sedi cell volume] because = ([old conc]*[cell volume]+[new mass])/[cell volume]
      !i.e. new concentration in sediment = [old conc] + (fallRate[m/s] * CellArea[m2] * timestep[s] * concentration [kg/m3])/[sedi cell volume] =
      !   [old conc] + (fallRate[m/s] * timestep[s] * concentration [kg/m3])/[sedi cell thickness]
!_KB_      sediC = sediC + self%w_d*self%dt*botC/self%sediTopDZ
!_KB_      sediM = sediM + self%w_d*self%dt*botM/self%sediTopDZ
!_KB_      sediW = sediW + self%w_d*self%dt*botW/self%sediTopDZ
   end if

   !functions
   botDZ = 10_rk
   sediFT = exp(self%qT*(botTemp-self%Tref))
   sediQWC = sediW/sediC
   if (sediQWC.ge.self%minQWC) then
      sediWr = sediFT*self%Wrmax*(1-self%minQWC/sediQWC)**2
   else
      sediWr = 0_rk
   end if
   if (botO2.gt.0) then
      sediRsorp = self%sorpRconst
      sediRdesorp = 0_rk
   else
      sediRsorp = 0_rk
      sediRdesorp = self%sorpRconst
   end if
   JPO4 = self%DPO4*(botPO4-sediPO4)/self%diffLength
   !ensure that diffusion does not become too large
   !JPO4=if(JPO4.gt.0,if(JPO4/botDZ*sediDt.gt.botPO4,botPO4/sediDt*botDZ,JPO4),if(sediPO4.lt.(sediWr*sediQWC-sediRsorp*sediPO4+sediRdesorp*mePO4+JPO4/sediTopDZ)*sediDt,(sediPO4/sediDt-sediWr*sediQWC+sediRsorp*sediPO4-sediRdesorp*mePO4)*sediTopDZ,JPO4))
   if (JPO4.gt.0) then
!_KB_      if (JPO4/botDZ*self%dt.gt.botPO4) then
!_KB_         JPO4 = botPO4/self%dt*botDZ
!_KB_      end if
   else
!_KB_      if (sediPO4.lt.(sediWr*sediQWC-sediRsorp*sediPO4+sediRdesorp*mePO4+JPO4/self%sediTopDZ)*self%dt) then
!_KB_         JPO4 = (sediPO4/self%dt-sediWr*sediQWC+sediRsorp*sediPO4-sediRdesorp*mePO4)*self%sediTopDZ
!_KB_      end if
   end if

   sediQMC = sediM/sediC
   sediNHr = sediFT*self%sediNHrmax*(botO2/(self%KO2NH+botO2))
   if (sediQMC.ge.self%minQMC) then
      sediMr = sediFT*self%Mrmax*(1-self%minQMC/sediQMC)**2
   else
      sediMr = 0_rk
   end if
   JNH4 = self%DNH4*(botNH4-sediNH4)/self%diffLength
   !ensure that diffusion does not become too large
   !JNH4=if(JNH4.gt.0,if(botNH4.lt.JNH4*sediDt/botDZ,botNH4*botDZ/sediDt,JNH4),if(sediNH4.lt.(sediNHr*sediNH4-sediMr*sediQMC-JNH4/sediTopDZ)*sediDt,sediNHr*sediNH4*sediTopDZ-sediNH4/sediDt*sediTopDZ-sediMr*sediQMC*sediTopDZ,JNH4))
!_KB_   if (JNH4.gt.0) then
!_KB_      if(botNH4.lt.JNH4*self%dt/botDZ) then
!_KB_         JNH4 = botNH4*botDZ/self%dt
!_KB_      end if
!_KB_   else
!_KB_      if (sediNH4.lt.(sediNHr*sediNH4-sediMr*sediQMC-JNH4/self%sediTopDZ)*self%dt) then
!_KB_          JNH4 = sediNHr*sediNH4*self%sediTopDZ-sediNH4/self%dt*self%sediTopDZ-sediMr*sediQMC*self%sediTopDZ
!_KB_      end if
!_KB_   end if

   JNO3 = self%DNO3*(botNO3-sediNO3)/self%diffLength
   !ensure that diffusion does not become too large
   !JNO3=if(JNO3.gt.0,if(botNO3.lt.JNO3*sediDt/botDZ,botNO3*botDZ/sediDt,JNO3),if(sediNO3.lt.(kDenit*sediNO3-sediNHr*sediNH4-JNO3/sediTopDZ)*sediDt,sediTopDZ*(sediNHr*sediNH4-sediNO3/sediDt-kDenit*sediNO3),JNO3))
!_KB_   if (JNO3.gt.0) then
!_KB_      if (botNO3.lt.JNO3*self%dt/botDZ) then
!_KB_         JNO3 = botNO3*botDZ/self%dt
!_KB_      end if
!_KB_   else
!_KB_      if (sediNO3.lt.(self%kDenit*sediNO3-sediNHr*sediNH4-JNO3/self%sediTopDZ)*self%dt) then
!_KB_         JNO3 = self%sediTopDZ*(sediNHr*sediNH4-sediNO3/self%dt-self%kDenit*sediNO3)
!_KB_      end if
!_KB_   end if

   if (sediQMC.ge.self%minQMC) then
      sediCr = sediFT*(botO2/(self%KO2C+botO2))*self%Crmax*(1-self%minQMC/sediQMC)**2
   else
      sediCr = 0_rk
   end if

    !!sediment state variables. Concentrations
    !_SET_ODE_BEN_(self%id_sediW, -sediWr*sediQWC)
    !_SET_ODE_BEN_(self%id_sediPO4, sediWr*sediQWC-sediRsorp*sediPO4+sediRdesorp*mePO4+JPO4/self%sediTopDZ)
    !_SET_ODE_BEN_(self%id_mePO4, sediRsorp*sediPO4-sediRdesorp*mePO4)
    !_SET_ODE_BEN_(self%id_sediM, -sediMr*sediQMC)
    !_SET_ODE_BEN_(self%id_sediN2, self%kDenit*sediNO3)
    !_SET_ODE_BEN_(self%id_sediNO3, sediNHr*sediNH4-self%kDenit*sediNO3+JNO3/self%sediTopDZ)
    !_SET_ODE_BEN_(self%id_sediNH4, sediMr*sediQMC-sediNHr*sediNH4+JNH4/self%sediTopDZ)
    !_SET_ODE_BEN_(self%id_sediC, sediCr*sediQMC)
    !!pelagic state variables fluxes
    !_SET_BOTTOM_EXCHANGE_(self%id_PO4, -JPO4)
    !_SET_BOTTOM_EXCHANGE_(self%id_NO3, -JNO3)
    !_SET_BOTTOM_EXCHANGE_(self%id_NH4, -JNH4)

!sediment state variables. Concentrations
   _SET_ODE_BEN_(self%id_sediW, 0)
   _SET_ODE_BEN_(self%id_sediPO4,0)
   _SET_ODE_BEN_(self%id_mePO4, 0)
   _SET_ODE_BEN_(self%id_sediM, 0)
   _SET_ODE_BEN_(self%id_sediN2, 0)
   _SET_ODE_BEN_(self%id_sediNO3, 0)
   _SET_ODE_BEN_(self%id_sediNH4, 0)
   _SET_ODE_BEN_(self%id_sediC, 0)

    !pelagic state variables fluxes
   _SET_BOTTOM_EXCHANGE_(self%id_PO4,0)
   _SET_BOTTOM_EXCHANGE_(self%id_NO3, 0)
   _SET_BOTTOM_EXCHANGE_(self%id_NH4, 0)

   ! Leave spatial loops over the horizontal domain (if any).
  _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_benthos
!EOC

!-----------------------------------------------------------------------

   end module fabm_au_sediment

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#endif
