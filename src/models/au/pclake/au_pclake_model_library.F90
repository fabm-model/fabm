!*******************************************************************************
!*                                                                             *
!* au_pclake_model_library.F90                                                 *
!*                                                                             *
!* Developed by :                                                              *
!*     Aquatic Ecology and Water Quality Management                            *
!*     Wageningen University,Netherland                                        *
!*     Bioscience-Silkeborg                                                    *
!*     Aarhus University,Denmark                                               *
!*     Center of Lake Restoration(CLEAR)                                       *
!*     Southern Denmark University,Denmark                                     *
!*     Copyright by the PCLake_FABM group @AU-silkeborg                        *
!*     under the GNU Public License - www.gnu.org                              *
!*                                                                             *
!*      Created March 2014                                                     *
!*-----------------------------------------------------------------------------*
!*-----------------------------------------------------------------------------*
!* If you have questions and suggestions regarding this model, please write to *
!*      Fenjuan Hu: fenjuan@bios.au.dk                                             *
!*      Dennis Trolle:dtr@bios.au.dk                                           *
!*      Karsten Bolding:bold@bios.au.dk                                        *
!*-----------------------------------------------------------------------------*
!*******************************************************************************

module au_pclake_model_library

   use fabm_types, only: type_base_model_factory,type_base_model
   
   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: au_pclake_model_factory

contains

   subroutine create(self,name,model)

      use au_pclake_abiotic_water
      use au_pclake_abiotic_sediment
      use au_pclake_phytoplankton_water
      use au_pclake_phytoplankton_sediment
      use au_pclake_macrophytes
      use au_pclake_foodweb_water
      use au_pclake_foodweb_sediment
      use au_pclake_auxilary
      use au_pclake_dummy_state
      

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

   
    !NULLIFY(model)
    
       select case (name)
       case ('au_pclake_abiotic_water');                allocate(type_au_pclake_abiotic_water::model);
       case ('au_pclake_abiotic_sediment');                allocate(type_au_pclake_abiotic_sediment::model);
       case ('au_pclake_phytoplankton_water');          allocate(type_au_pclake_phytoplankton_water::model);
       case ('au_pclake_phytoplankton_sediment');          allocate(type_au_pclake_phytoplankton_sediment::model);
       case ('au_pclake_macrophytes');                allocate(type_au_pclake_macrophytes::model);
       case ('au_pclake_foodweb_water');                allocate(type_au_pclake_foodweb_water::model);
       case ('au_pclake_foodweb_sediment');                allocate(type_au_pclake_foodweb_sediment::model);
       case ('au_pclake_auxilary');                   allocate(type_au_pclake_auxilary::model);
       case ('au_pclake_dummy_state');                   allocate(type_au_pclake_dummy_state::model);
       
       !case default
           
           
       end select

   end subroutine
!EOC

!------------------------------------------------------------------------------
!
!EOP
!------------------------------------------------------------------------------
end module
!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
