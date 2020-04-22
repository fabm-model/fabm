#include "fabm_driver.h"

!---------------------------------------
!---CRi PREDATION FUNCTION --------------
!---------------------------------------

module su_ersem_mesozooplankton

   use fabm_types
   use fabm_particle
   use fabm_expressions
   use fabm_builtin_models

   use ersem_shared
   use ersem_pelagic_base

   implicit none

   private

   type,extends(type_ersem_pelagic_base),public :: type_su_ersem_mesozooplankton
      ! Variables
      type (type_model_id),               allocatable,dimension(:) :: id_prey
      type (type_dependency_id),          allocatable,dimension(:) :: id_preyc,id_preyn,id_preyp,id_preys,id_preyf,id_preyl
      type (type_state_variable_id),      allocatable,dimension(:) :: id_preyf_target
      type (type_state_variable_id)      :: id_O3c,id_O2o,id_L2c,id_TA
      type (type_state_variable_id)      :: id_R1c,id_R1p,id_R1n
      type (type_state_variable_id)      :: id_R2c
      type (type_state_variable_id)      :: id_RPc,id_RPp,id_RPn,id_RPs
      type (type_state_variable_id)      :: id_N1p,id_N4n
      type (type_dependency_id)          :: id_ETW,id_eO2mO2
      type (type_dependency_id)          :: id_totprey
      type (type_horizontal_dependency_id) :: id_inttotprey

      type (type_diagnostic_variable_id) :: id_fZIO3c, id_SZIc, id_rug, id_fZIRDc, id_fZIRPc, id_sd, id_excess_n, id_excess_p, id_I
      type (type_diagnostic_variable_id),allocatable,dimension(:) :: id_sprey, id_fpreyc, id_Cp_PB, id_Ii

      ! Parameters
      integer  :: nprey
      real(rk) :: qpc,qnc
      real(rk) :: q10,chro,minfood,chuc
      real(rk),allocatable :: Cr(:),pu_ea(:)
      real(rk) :: sum
      real(rk) :: sdo, sd, srs
      real(rk) :: pu
      real(rk) :: pe_R1
      real(rk) :: gutdiss
      real(rk) :: xR1n,xR1p
      real(rk) :: Imax, KI   ! parameters for the new predation function

      real(rk) :: Minprey,repw,mort

      ! ERSEM global parameters
      real(rk) :: R1R2,urB1_O2
   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_su_ersem_mesozooplankton),intent(inout),target :: self
      integer,                              intent(in)           :: configunit
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
      integer           :: iprey
      character(len=16) :: index
      type (type_weighted_sum),pointer :: total_prey_calculator
      real(rk)          :: pu_ea,pu_eaR
      logical           :: preyispom
      real(rk)          :: c0
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%q10,    'q10',    '-',          'Q_10 temperature coefficient')
      call self%get_parameter(self%minfood,'minfood','mg C/m^3',   'Michaelis-Menten constant to perceive food')
      call self%get_parameter(self%chuc,   'chuc',   'mg C/m^3',   'Michaelis-Menten constant for food uptake')
      call self%get_parameter(self%sum,    'sum',    '1/d',        'maximum specific uptake at reference temperature')
      call self%get_parameter(self%pu,     'pu',     '-',          'assimilation efficiency')
      call self%get_parameter(pu_ea,       'pu_ea',  '-',          'fraction of unassimilated prey that is excreted (not respired)')
      call self%get_parameter(pu_eaR,      'pu_eaR', '-',          'fraction of unassimilated detritus that is excreted (not respired)')
      call self%get_parameter(self%pe_R1,  'pe_R1',  '-',          'dissolved fraction of excreted/dying matter')
      call self%get_parameter(self%srs,    'srs',    '1/d',        'specific rest respiration at reference temperature')
      call self%get_parameter(self%sd,     'sd',     '1/d',        'basal mortality')
      call self%get_parameter(self%sdo,    'sdo',    '1/d',        'maximum mortality due to oxygen limitation')
      call self%get_parameter(self%chro,   'chro',   '-',          'Michaelis-Menten constant for oxygen limitation')
      call self%get_parameter(self%qpc,    'qpc',    'mmol P/mg C','phosphorus to carbon ratio')
      call self%get_parameter(self%qnc,    'qnc',    'mmol N/mg C','nitrogen to carbon ratio')
      call self%get_parameter(self%Minprey,'Minprey','mg C/m^2',   'food threshold for overwintering state')
      call self%get_parameter(self%repw,   'repw',   '1/d',        'specific overwintering respiration')
      call self%get_parameter(self%mort,   'mort',   '1/d',        'specific overwintering mortality')
      call self%get_parameter(self%Imax,   'Imax',   'gC/gC/d',    'maximum ingestion rate')
      call self%get_parameter(self%KI,     'KI',     'dl',         'half saturation constant for satiation')

      call self%get_parameter(self%R1R2,   'R1R2','-','labile fraction of produced dissolved organic carbon')
      call self%get_parameter(self%xR1p,   'xR1p','-','transfer of phosphorus to DOM, relative to POM',default=0._rk)
      call self%get_parameter(self%xR1n,   'xR1n','-','transfer of nitrogen to DOM, relative to POM',default=0._rk)
      call self%get_parameter(self%urB1_O2,'urB1_O2','mmol O_2/mg C','oxygen consumed per carbon respired')

      call self%get_parameter(self%gutdiss,'gutdiss','-','fraction of prey calcite that dissolves after ingestion')

      call self%get_parameter(c0,'c0','mg C/m^3','background concentration')

      ! Register state variables
      call self%initialize_ersem_base(sedimentation=.false.)
      call self%add_constituent('c',1.e-4_rk,c0,qn=self%qnc,qp=self%qpc)

      ! Determine number of prey types.
      call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)

      ! Get prey-specific parameters.
      allocate(self%Cr(self%nprey))
      allocate(self%pu_ea(self%nprey))
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%get_parameter(self%Cr(iprey),'Cr'//trim(index),'-','clearance rate on prey i'//trim(index))
      end do
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%get_parameter(preyispom,'prey'//trim(index)//'ispom','','prey type '//trim(index)//' is detritus',default=.false.)
         if (preyispom) then
            self%pu_ea(iprey) = pu_eaR
         else
            self%pu_ea(iprey) = pu_ea
         end if
      end do

      ! Get prey-specific coupling links.
      allocate(self%id_prey(self%nprey))
      allocate(self%id_preyc(self%nprey))
      allocate(self%id_preyn(self%nprey))
      allocate(self%id_preyp(self%nprey))
      allocate(self%id_preyf(self%nprey))
      allocate(self%id_preys(self%nprey))
      allocate(self%id_preyl(self%nprey))
      allocate(self%id_preyf_target(self%nprey))
      allocate(self%id_sprey(self%nprey))
      allocate(self%id_fpreyc(self%nprey))
      allocate(self%id_Cp_PB(self%nprey))
      allocate(self%id_Ii(self%nprey))
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call self%register_dependency(self%id_preyc(iprey),'prey'//trim(index)//'c','mmol C/m^3', 'prey '//trim(index)//' carbon')
         call self%register_dependency(self%id_preyn(iprey),'prey'//trim(index)//'n','mmol N/m^3', 'prey '//trim(index)//' nitrogen')
         call self%register_dependency(self%id_preyp(iprey),'prey'//trim(index)//'p','mmol P/m^3', 'prey '//trim(index)//' phosphorus')
         call self%register_dependency(self%id_preys(iprey),'prey'//trim(index)//'s','mmol Si/m^3','prey '//trim(index)//' silicate')
         call self%register_dependency(self%id_preyl(iprey),'prey'//trim(index)//'l','mg C/m^3',   'prey '//trim(index)//' calcite')

         call self%register_diagnostic_variable(self%id_sprey(iprey), 'sprey'//trim(index), 'day-1',    'specific predation rate'//trim(index))
         call self%register_diagnostic_variable(self%id_fpreyc(iprey),'fpreyc'//trim(index),'mgC/m3/d', 'carbon predation rate'//trim(index))
         call self%register_diagnostic_variable(self%id_Cp_PB(iprey),'Cp_PB'//trim(index),'', 'akin to clearance rate rate'//trim(index))
         call self%register_diagnostic_variable(self%id_Ii(iprey),'Ii'//trim(index),'', 'carbon-specific ingestion rate'//trim(index))

         call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
         call self%request_coupling_to_model(self%id_preyc(iprey),self%id_prey(iprey),standard_variables%total_carbon)
         call self%request_coupling_to_model(self%id_preyn(iprey),self%id_prey(iprey),standard_variables%total_nitrogen)
         call self%request_coupling_to_model(self%id_preyp(iprey),self%id_prey(iprey),standard_variables%total_phosphorus)
         call self%request_coupling_to_model(self%id_preys(iprey),self%id_prey(iprey),standard_variables%total_silicate)
         call self%request_coupling_to_model(self%id_preyl(iprey),self%id_prey(iprey),total_calcite_in_biota)

         if (use_iron) then
            call self%register_dependency(self%id_preyf(iprey),'prey'//trim(index)//'f','umol Fe/m^3','prey '//trim(index)//' iron')
            call self%register_state_dependency(self%id_preyf_target(iprey),'prey'//trim(index)//'f_sink','umol Fe/m^3','sink for Fe of prey '//trim(index),required=.false.)
            call self%request_coupling_to_model(self%id_preyf(iprey),self%id_prey(iprey),standard_variables%total_iron)
         end if
      end do

      ! Create a submodel that will compute total prey for us, and create a variable that will contain its depth integral.
      ! This quantity will be depth integrated to determine whether we should be overwintering.
      allocate(total_prey_calculator)
      total_prey_calculator%units = 'mg C/m^3'
      total_prey_calculator%result_output = output_none
      do iprey=1,self%nprey
         write (index,'(i0)') iprey
         call total_prey_calculator%add_component('prey'//trim(index)//'c',self%Cr(iprey))
      end do
      call self%add_child(total_prey_calculator,'totprey_calculator',configunit=-1)

      ! Create a link to the total prey, and request its depth-integrated value.
      call self%register_dependency(self%id_totprey,'totprey','mg C/m^3','total carbon in prey')
      call self%request_coupling(self%id_totprey,'totprey_calculator/result')
      call self%register_expression_dependency(self%id_inttotprey,vertical_integral(self%id_totprey))

      ! Register links to external nutrient pools.
      call self%register_state_dependency(self%id_N1p,'N1p','mmol P/m^3','phosphate')
      call self%register_state_dependency(self%id_N4n,'N4n','mmol N/m^3','ammonium')

      ! Register links to external labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R1c,'R1c','mg C/m^3',  'dissolved organic carbon')
      call self%register_state_dependency(self%id_R1p,'R1p','mmol P/m^3','dissolved organic phosphorus')
      call self%register_state_dependency(self%id_R1n,'R1n','mmol N/m^3','dissolved organic nitrogen')

      ! Register links to external semi-labile dissolved organic matter pools.
      call self%register_state_dependency(self%id_R2c,'R2c','mg C/m^3','semi-labile dissolved organic carbon')

      ! Register links to external particulate organic matter pools.
      call self%register_state_dependency(self%id_RPc,'RPc','mg C/m^3',   'particulate organic carbon')
      call self%register_state_dependency(self%id_RPp,'RPp','mmol P/m^3', 'particulate organic phosphorus')
      call self%register_state_dependency(self%id_RPn,'RPn','mmol N/m^3', 'particulate organic nitrogen')
      call self%register_state_dependency(self%id_RPs,'RPs','mmol Si/m^3','particulate organic silicate')

      ! Allow coupling of all required particulate organic matter variables to a single source model.
      call self%request_coupling_to_model(self%id_RPc,'RP','c')
      call self%request_coupling_to_model(self%id_RPn,'RP','n')
      call self%request_coupling_to_model(self%id_RPp,'RP','p')
      call self%request_coupling_to_model(self%id_RPs,'RP','s')

      ! Register links to external total dissolved inorganic carbon, dissolved oxygen pools
      call self%register_state_dependency(self%id_O2o,'O2o','mmol O_2/m^3','oxygen source')
      call self%register_state_dependency(self%id_O3c,'O3c','mmol C/m^3','carbon dioxide sink')
      call self%register_state_dependency(self%id_TA,standard_variables%alkalinity_expressed_as_mole_equivalent)    

      call self%register_state_dependency(self%id_L2c,'L2c','mg C/m^3','calcite',required=.false.)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_ETW,standard_variables%temperature)
      call self%register_dependency(self%id_eO2mO2,standard_variables%fractional_saturation_of_oxygen)

      ! Register diagnostics
      call self%register_diagnostic_variable(self%id_fZIO3c,   'fZIO3c',   'mg C/m^3/d', 'respiration',             output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_SZIc,     'SZIc',     'mg C/m^3/d', 'growth',                  output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_rug,      'rug',      'mg C/m^3/d', 'ingestion',               output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_fZIRDc,   'fZIRDc',   'mg C/m^3/d', 'excretion dissolved',     output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_fZIRPc,   'fZIRPc',   'mg C/m^3/d', 'excretion particulated',  output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_sd,       'sd',       '1/d',        'specific mortality',      output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_excess_n, 'excess_n', 'mg C/m^3/d', 'excess nitrogen',         output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_excess_p, 'excess_p', 'mg C/m^3/d', 'excess phosphorus',       output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_I,        'I',        'g C/g C/d',  'total carbon ingestion',  output=output_time_step_averaged)

      ! Contribute to aggregate fluxes.
      call self%add_to_aggregate_variable(type_bulk_standard_variable(name='zooplankton_respiration_rate',units='mg C/m^3/d',aggregate_variable=.true.),self%id_fZIO3c)

   end subroutine

   subroutine do(self,_ARGUMENTS_DO_)

      class (type_su_ersem_mesozooplankton),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      integer  :: iprey,istate
      real(rk),dimension(self%nprey) :: Cp_PB, Ii     ! Implementing new predation function
      real(rk) :: I                         ! Implementing new predation function
      real(rk) :: ETW,eO2mO2
      real(rk) :: c,cP
      real(rk),dimension(self%nprey) :: preycP,preypP,preynP,preysP,preylP
      real(rk),dimension(self%nprey) :: sprey,rupreyc,fpreyc
      real(rk) :: preyP
      real(rk) :: et,CORROX,eO2
      real(rk) :: rum,put_u,rug
      real(rk) :: sd,rd
      real(rk) :: ineff
      real(rk) :: rrs,rra
      real(rk) :: fZIO3c
      real(rk) :: ret,fZIRDc,fZIRPc
      real(rk) :: fZIRIp,fZIRDp,fZIRPp
      real(rk) :: fZIRIn,fZIRDn,fZIRPn
      real(rk) :: SZIc,SZIn,SZIp
      real(rk) :: excess_c, excess_n, excess_p
      real(rk) :: intprey

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_inttotprey,intprey)
         intprey = intprey*CMass

         if (intprey>self%Minprey) then
            ! Enough prey available - not overwintering

            ! Get environment (temperature, oxygen saturation)
            _GET_(self%id_ETW,ETW)
            _GET_(self%id_eO2mO2,eO2mO2)
            eO2mO2 = min(1.0_rk,eO2mO2)

            ! Get own concentrations: c includes background concentration and is used in source terms,
            ! cP excludes background concentration and is used in sink terms.
            _GET_WITH_BACKGROUND_(self%id_c,c)
            _GET_(self%id_c,cP)

            ! Get prey concentrations
            do iprey=1,self%nprey
               _GET_(self%id_preyc(iprey), preycP(iprey))
               _GET_(self%id_preyp(iprey), preypP(iprey))
               _GET_(self%id_preyn(iprey), preynP(iprey))
               _GET_(self%id_preys(iprey), preysP(iprey))
               _GET_(self%id_preyl(iprey), preylP(iprey))
            end do

            ! Prey carbon was returned in mmol (due to units of standard_variables%total_carbon); convert to mg
            preycP = preycP*CMass

            ! Temperature effect:
            et = self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk)

            ! Oxygen limitation (based on oxygen saturation eO2mO2):
            CORROX = 1._rk + self%chro
            eO2 = MIN(1._rk,CORROX*(eO2mO2/(self%chro + eO2mO2)))


!-------------------------------------------------------------
!-----------Implementing the Cpi predation function-----------
!-------------------------------------------------------------

            ! Potential capture/ingestion rate of prey i (gCprey/gCpredator/day):
            Cp_PB = self%Cr * preycP *et

            ! Total ingestion rate (gCprey/gCpredator/day): 
            I = MIN(self%Imax * sum(Cp_PB)/(sum(Cp_PB) + self%KI), sum(Cp_PB))

            ! Ingestion rate of prey i (gCprey/gCpredator/day):
            Ii = I * Cp_PB/sum(Cp_PB)

            ! Specific ingestion rate on prey i (day-1):
            do iprey=1,self%nprey
               if (preycP(iprey) > 0._rk) then
                  sprey(iprey) = Ii(iprey) * c / preycP(iprey)
               else 
                  sprey(iprey) = 0._rk
               end if
            end do

            ! Total ingestion rate in ug C prey/L/day
            rug  =   sum(sprey * preycP) + 1e-60_rk

            ! Ingestion rate of prey i in ugC prey/L/day:
            fpreyc = sprey * preycP    

            ! Apply specific predation rates to all state variables of every prey.
            do iprey=1,self%nprey
               do istate=1,size(self%id_prey(iprey)%state)
                  _GET_(self%id_prey(iprey)%state(istate),preyP)
                  _ADD_SOURCE_(self%id_prey(iprey)%state(istate),-sprey(iprey)*preyP)
               end do
            end do

!--------------------------------------------------------------------------------------
!-------------End of implementation of the predation function Cpi----------------------
!--------------------------------------------------------------------------------------

            ! Specific mortality (1/d): background mortality + mortality due to oxygen limitation.
            sd = self%sd + (1._rk - eO2)*self%sdo

            ! Compute abolute mortality (mg C/m3/d) from specific mortality and biomass.
            rd = sd*cP

            ! Assimilation inefficiency (dimensionless):
            ineff = 1._rk - self%pu

            ! Excretion and egestion of organic matter (part dissolved, part particulate)
            ret = ineff * sum(fpreyc*self%pu_ea)
            fZIRDc = (ret + rd)*self%pe_R1
            fZIRPc = (ret + rd)*(1._rk - self%pe_R1)

            ! Rest respiration (mg C/m3/d), corrected for prevailing temperature
            rrs = self%srs*et*cP

            ! Activity respiration (mg C/m3/d)
            rra = ineff * rug - ret

            ! Total respiration (mg C/m3/d)
            fZIO3c = rrs + rra

            if (_AVAILABLE_(self%id_L2c)) then
               ! Send non-dissolved fraction of consumed calcite to external calcite pool.
               ! As the consumed calcite has not yet been taken away from the dissolved
               ! inorganic carbon pool (prey calcite is "virtual calcite", materializing only
               ! the moment the prey dies), do so now.
               _ADD_SOURCE_(self%id_L2c,  (1.0_rk-self%gutdiss)*ineff*sum(self%pu_ea*sprey*preylP))
               _ADD_SOURCE_(self%id_O3c, -(1.0_rk-self%gutdiss)*ineff*sum(self%pu_ea*sprey*preylP)/CMass)
               _ADD_SOURCE_(self%id_TA,-2*(1.0_rk-self%gutdiss)*ineff*sum(self%pu_ea*sprey*preylP)/CMass)   ! CaCO3 formation decreases alkalinity by 2 units
            end if

            ! Source equation for carbon in biomass (NB cannibalism is handled as part of predation formulation)
            SZIc = rug - fZIRDc - fZIRPc - fZIO3c

            ! Carbon flux to labile dissolved, refractory dissolved, and particulate organic matter.
            _ADD_SOURCE_(self%id_R1c, + fZIRDc * self%R1R2)
            _ADD_SOURCE_(self%id_R2c, + fZIRDc * (1._rk-self%R1R2))
            _ADD_SOURCE_(self%id_RPc, + fZIRPc)

            ! Account for CO2 production and oxygen consumption in respiration.
            _ADD_SOURCE_(self%id_O3c, + fZIO3c/CMass)
            _ADD_SOURCE_(self%id_O2o, - fZIO3c*self%urB1_O2)

            ! -------------------------------
            ! Phosphorus
            ! -------------------------------

            ! Phosphorus loss to dissolved and particulate organic matter.
            fZIRIp = (fZIRDc + fZIRPc) * self%qpc
            fZIRDp = min(fZIRIp, fZIRDc * self%qpc * self%xR1p)
            fZIRPp = fZIRIp - fZIRDp

            ! Source equation for phosphorus in biomass
            SZIp = sum(sprey*preypP) - fZIRPp - fZIRDp

            ! Phosphorus flux to dissolved and particulate organic matter.
            _ADD_SOURCE_(self%id_R1p, + fZIRDp)
            _ADD_SOURCE_(self%id_RPp, + fZIRPp)

            ! -------------------------------
            ! Nitrogen
            ! -------------------------------

            ! Nitrogen loss to dissolved and particulate organic matter.
            fZIRIn = (fZIRDc + fZIRPc) * self%qnc
            fZIRDn = min(fZIRIn, fZIRDc * self%qnc * self%xR1n)
            fZIRPn = fZIRIn - fZIRDn

            ! Source equation for nitrogen in biomass
            SZIn = sum(sprey*preynP) - fZIRPn - fZIRDn

            ! Nitrogen flux to dissolved and particulate organic matter.
            _ADD_SOURCE_(self%id_R1n, + fZIRDn)
            _ADD_SOURCE_(self%id_RPn, + fZIRPn)

            ! -------------------------------
            ! Silicate
            ! -------------------------------

            ! Send all consumed silicate to particulate organic matter pool.
            _ADD_SOURCE_(self%id_RPs,sum(sprey*preysP))

            ! -------------------------------
            ! Iron
            ! -------------------------------

            if (use_iron) then
               ! Iron dynamics:
               ! Following Vichi et al., 2007 it is assumed that the iron fraction of the ingested phytoplankton
               ! is egested as particulate detritus (Luca)
               do iprey=1,self%nprey
                  _GET_(self%id_preyf(iprey),preyP)
                  if (preyP/=0.0_rk) _ADD_SOURCE_(self%id_preyf_target(iprey),sprey(iprey)*preyP)
               end do
            end if

            ! -------------------------------
            ! Maintain constant stoichiometry
            ! -------------------------------
            
            ! Compute excess carbon flux, given that the maximum realizable carbon flux needs to be balanced
            ! by corresponding nitrogen and phosphorus fluxes to maintain constant stoichiometry.
            excess_c = max(max(SZIc - SZIp/self%qpc,SZIc - SZIn/self%qnc),0._rk)
            SZIc = SZIc - excess_c
            _ADD_SOURCE_(self%id_c,SZIc)
            _ADD_SOURCE_(self%id_RPc,excess_c)

            ! Compute excess nitrogen and phosphorus fluxes, based on final carbon flux.
            excess_n = max(SZIn - SZIc*self%qnc,0.0_rk)
            excess_p = max(SZIp - SZIc*self%qpc,0.0_rk)

            ! Send excess nitrogen and phosphorus to ammonium and phosphate, respectively.
            _ADD_SOURCE_(self%id_N4n,excess_n)
            _ADD_SOURCE_(self%id_N1p,excess_p)
            _ADD_SOURCE_(self%id_TA,excess_n-excess_p)  ! Alkalinity contributions: +1 for NH4, -1 for PO4

         else

            ! Insufficient prey - overwintering
            _GET_(self%id_c,cP)

            ! Respiration
            fZIO3c = cP * self%repw

            ! Mortality
            fZIRPc = cP * self%mort

            _ADD_SOURCE_(self%id_RPc,fZIRPc)
            _ADD_SOURCE_(self%id_RPn,fZIRPc*self%qnc)
            _ADD_SOURCE_(self%id_RPp,fZIRPc*self%qpc)

            _ADD_SOURCE_(self%id_O3c,fZIO3c/CMass)
            _ADD_SOURCE_(self%id_N4n,fZIO3c*self%qnc)
            _ADD_SOURCE_(self%id_N1p,fZIO3c*self%qpc)
            _ADD_SOURCE_(self%id_TA, fZIO3c*(self%qnc-self%qpc))  ! Alkalinity contributions: +1 for NH4, -1 for PO4

            _ADD_SOURCE_(self%id_c,- fZIRPc - fZIO3c)

            SZIc = 0
            rug = 0
            fZIRDc = 0
            sd = 0
            excess_n = 0
            excess_p = 0
            I = 0
            sprey = 0
            fpreyc = 0
            Cp_PB = 0
            Ii = 0
         end if

         _SET_DIAGNOSTIC_(self%id_fZIO3c,   fZIO3c)
         _SET_DIAGNOSTIC_(self%id_SZIc,     SZIc)
         _SET_DIAGNOSTIC_(self%id_rug,      rug)
         _SET_DIAGNOSTIC_(self%id_fZIRDc,   fZIRDc)
         _SET_DIAGNOSTIC_(self%id_fZIRPc,   fZIRPc)
         _SET_DIAGNOSTIC_(self%id_sd,       sd)
         _SET_DIAGNOSTIC_(self%id_excess_n, excess_n)
         _SET_DIAGNOSTIC_(self%id_excess_p, excess_p)
         _SET_DIAGNOSTIC_(self%id_I,        I)
         do iprey=1,self%nprey
            _SET_DIAGNOSTIC_(self%id_sprey(iprey), sprey(iprey))
            _SET_DIAGNOSTIC_(self%id_fpreyc(iprey), fpreyc(iprey))
            _SET_DIAGNOSTIC_(self%id_Cp_PB(iprey), Cp_PB(iprey))
            _SET_DIAGNOSTIC_(self%id_Ii(iprey), Ii(iprey))
         end do

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

end module
