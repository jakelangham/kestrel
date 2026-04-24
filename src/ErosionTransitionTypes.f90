module erosion_transition_types

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use messages_module, only: FatalErrorMessage, WarningMessage
   use dict_module, only: Dict
   use erosion_transition_interface, only: ErosionTransition

   implicit none

   private
   public :: NoErosionTransition, SmoothErosionTransition, StepErosionTransition
   public :: create_erosion_transition

   type, extends(ErosionTransition) :: SmoothErosionTransition
      real(kind=wp) :: rate
   contains
      procedure :: validate => SmoothErosionTransition_validate
      procedure :: transition => SmoothErosionTransition_transition
   end type

   interface SmoothErosionTransition
      procedure SmoothErosionTransition_constructor
   end interface SmoothErosionTransition

   interface
      module function SmoothErosionTransition_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(SmoothErosionTransition) :: this
      end function SmoothErosionTransition_constructor

      module subroutine SmoothErosionTransition_validate(this)
         class(SmoothErosionTransition), intent(in) :: this
      end subroutine SmoothErosionTransition_validate

      module pure function SmoothErosionTransition_transition(this, Bt) result(factor)
         class(SmoothErosionTransition), intent(in) :: this
         real(kind=wp), intent(in) :: Bt
         real(kind=wp) :: factor
      end function SmoothErosionTransition_transition
   end interface

   type, extends(ErosionTransition) :: StepErosionTransition
   contains
      procedure :: validate => StepErosionTransition_validate
      procedure :: transition => StepErosionTransition_transition
   end type

   interface StepErosionTransition
      procedure StepErosionTransition_constructor
   end interface StepErosionTransition

   interface
      module function StepErosionTransition_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(StepErosionTransition) :: this
      end function StepErosionTransition_constructor

      module subroutine StepErosionTransition_validate(this)
         class(StepErosionTransition), intent(in) :: this
      end subroutine StepErosionTransition_validate

      module pure function StepErosionTransition_transition(this, Bt) result(factor)
         class(StepErosionTransition), intent(in) :: this
         real(kind=wp), intent(in) :: Bt
         real(kind=wp) :: factor
      end function StepErosionTransition_transition
   end interface

   type, extends(ErosionTransition) :: NoErosionTransition
   contains
      procedure :: validate => NoErosionTransition_validate
      procedure :: transition => NoErosionTransition_transition
   end type

   interface NoErosionTransition
      procedure NoErosionTransition_constructor
   end interface NoErosionTransition

   interface
      module function NoErosionTransition_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(NoErosionTransition) :: this
      end function NoErosionTransition_constructor

      module subroutine NoErosionTransition_validate(this)
         class(NoErosionTransition), intent(in) :: this
      end subroutine NoErosionTransition_validate

      module pure function NoErosionTransition_transition(this, Bt) result(factor)
         class(NoErosionTransition), intent(in) :: this
         real(kind=wp), intent(in) :: Bt
         real(kind=wp) :: factor
      end function NoErosionTransition_transition
   end interface

contains

   function create_erosion_transition(choice, ParamsDict) result(erosion_transition)
      type(varString), intent(in) :: choice
      type(Dict), intent(in) :: ParamsDict
      class(ErosionTransition), allocatable :: erosion_transition

      if (choice%len()==0) then
         call WarningMessage("In the 'Parameters' block 'erosion transition' is not given." // new_line('a') // &
                                      " Using default 'erosion transition = smooth'")
         allocate(erosion_transition, source = SmoothErosionTransition(ParamsDict))
         return
      end if

      select case (choice%s)
         case ('smooth')
            allocate(erosion_transition, source=SmoothErosionTransition(ParamsDict))
         case ('step')
            allocate(erosion_transition, source=StepErosionTransition(ParamsDict))
         case ('none')
            allocate(erosion_transition, source=NoErosionTransition(ParamsDict))
         case default
            call WarningMessage("In the 'Parameters' block 'erosion transition = " // choice%s // "' is not recognized." // new_line('a') // &
                                " Using default 'erosion transition = Smooth'")
            allocate(erosion_transition, source=SmoothErosionTransition(ParamsDict))
      end select
      
   end function create_erosion_transition

end module erosion_transition_types