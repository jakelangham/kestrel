module deposition_types

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use messages_module, only: FatalErrorMessage, WarningMessage
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use deposition_interface, only: DepositionModel

   implicit none

   private
   public :: NoDeposition, SimpleHinderedSettling, SpearmanManningHinderedSettling
   public :: create_deposition_model

   ! -- Specific deposition model interfaces --
   !
   !    The 'validate' function does validation of model parameters
   !
   !    Implementations of the constructor, validator and hindering are in submodules.

   ! No deposition
   type, extends(DepositionModel) :: NoDeposition
      contains
         procedure :: validate => NoDeposition_validator
         procedure :: hindering => NoDeposition_hindering
   end type NoDeposition

   interface NoDeposition
      procedure NoDeposition_constructor
   end interface NoDeposition

   interface
      module function NoDeposition_constructor(params, RunParams) result(this)
         type(Dict), intent(in) :: params
         type(RunSet), intent(in) :: RunParams
         type(NoDeposition) :: this
      end function NoDeposition_constructor

      module subroutine NoDeposition_validator(this)
         class(NoDeposition), intent(in) :: this
      end subroutine NoDeposition_validator

      module pure function NoDeposition_hindering(this,psi) result(f)
         class(NoDeposition), intent(in) :: this
         real(kind=wp), intent(in) :: psi
         real(kind=wp) :: f
      end function NoDeposition_hindering
   end interface

   ! Simple quadratic hindered settling function, with zeros at 
   ! psi = 0 and psi = maximum packing
   type, extends(DepositionModel) :: SimpleHinderedSettling
      contains
         procedure :: validate => SimpleHinderedSettling_validator
         procedure :: hindering => SimpleHinderedSettling_hindering
   end type SimpleHinderedSettling

   interface SimpleHinderedSettling
      procedure SimpleHinderedSettling_constructor
   end interface SimpleHinderedSettling

   interface
      module function SimpleHinderedSettling_constructor(params, RunParams) result(this)
         type(Dict), intent(in) :: params
         type(RunSet), intent(in) :: RunParams
         type(SimpleHinderedSettling) :: this
      end function SimpleHinderedSettling_constructor

      module subroutine SimpleHinderedSettling_validator(this)
         class(SimpleHinderedSettling), intent(in) :: this
      end subroutine SimpleHinderedSettling_validator

      module pure function SimpleHinderedSettling_hindering(this,psi) result(f)
         class(SimpleHinderedSettling), intent(in) :: this
         real(kind=wp), intent(in) :: psi
         real(kind=wp) :: f
      end function SimpleHinderedSettling_hindering
   end interface

   ! Hindered settling function from Spearman & Manning (2017)
   ! doi:10.1007/s10236-017-1034-7.
   type, extends(DepositionModel) :: SpearmanManningHinderedSettling
      real(kind=wp) :: nsettling
      real(kind=wp) :: a
      real(kind=wp) :: b
      contains
         procedure :: validate => SpearmanManningHinderedSettling_validator
         procedure :: hindering => SpearmanManningHinderedSettling_hindering
   end type SpearmanManningHinderedSettling

   interface SpearmanManningHinderedSettling
      procedure SpearmanManningHinderedSettling_constructor
   end interface SpearmanManningHinderedSettling

   interface
      module function SpearmanManningHinderedSettling_constructor(params, RunParams) result(this)
         type(Dict), intent(in) :: params
         type(RunSet), intent(in) :: RunParams
         type(SpearmanManningHinderedSettling) :: this
      end function SpearmanManningHinderedSettling_constructor

      module subroutine SpearmanManningHinderedSettling_validator(this)
         class(SpearmanManningHinderedSettling), intent(in) :: this
      end subroutine SpearmanManningHinderedSettling_validator

      module pure function SpearmanManningHinderedSettling_hindering(this,psi) result(f)
         class(SpearmanManningHinderedSettling), intent(in) :: this
         real(kind=wp), intent(in) :: psi
         real(kind=wp) :: f
      end function SpearmanManningHinderedSettling_hindering
   end interface

contains
   
   subroutine create_deposition_model(choice, ParamsDict, RunParams, deposition_model)
      character(len=*), intent(in) :: choice
      type(Dict), intent(in) :: ParamsDict
      type(RunSet), intent(in) :: RunParams
      class(DepositionModel), allocatable, intent(out) :: deposition_model

      if (len_trim(choice) == 0) then
         call WarningMessage("In the 'Parameters' block 'deposition' is not given." // new_line('A') // " Using default 'Spearman Manning' deposition model.")
         allocate(deposition_model, source =  SpearmanManningHinderedSettling(ParamsDict, RunParams))
         return
      end if

      select case (trim(choice))
         case ('none')
            allocate(deposition_model, source = NoDeposition(ParamsDict, RunParams))
         case ('simple')
            allocate(deposition_model, source = SimpleHinderedSettling(ParamsDict, RunParams))
         case ('spearman manning')
            allocate(deposition_model, source =  SpearmanManningHinderedSettling(ParamsDict, RunParams))
         case default
            call WarningMessage("In the 'Parameters' block 'deposition = " // choice // "' is not recognized." // new_line('A') // " Using default 'Spearman Manning' deposition model.")
            allocate(deposition_model, source =  SpearmanManningHinderedSettling(ParamsDict, RunParams))
      end select
   end subroutine create_deposition_model
  
end module deposition_types