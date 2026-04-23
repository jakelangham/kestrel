module erosion_transition_interface
! -- Transition functions from erodible to inerodible at the maximum erosion
!    depth (erosion_depth). --
!
!    The ErosionTransition class specifies how to transition when the elevation decrease
!    approaches the erosion depth.
!
!    The 'transition' functions have a common signature, taking the maximum erosion depth (erosion_depth)
!    and current change in elevation (Bt) as inputs and returning a factor (0<= factor <=1)
!    which is used in other places to transition the erosion rate etc. as erosion_depth is approached.
!
!    Specific transition functions extend the ErosionTransition class
!    and are defined in ErosionTransitionTypes.f90
!    with implementations in submodules in src/erosion_transition_models

   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage
   use varstring_module, only: varString
   use dict_module, only: Dict

   implicit none

   private
   public :: ErosionTransition, erosion_transition

   type, abstract :: ErosionTransition
      type(varString) :: name
      real(kind=wp) :: erosion_depth
   contains
      procedure, non_overridable :: initialize => ErosionTransition_constructor
      procedure, non_overridable :: base_validate => ErosionTransition_validate
      procedure(validate_transition), deferred :: validate
      procedure(transition_func), deferred :: transition
   end type ErosionTransition

   abstract interface
      subroutine validate_transition(this)
         import :: ErosionTransition
         class(ErosionTransition), intent(in) :: this
      end subroutine validate_transition
   end interface

   abstract interface
      pure function transition_func(this, Bt) result(factor)
         import :: ErosionTransition
         import :: wp
         class(ErosionTransition), intent(in) :: this
         real(kind=wp), intent(in) :: Bt
         real(kind=wp) :: factor
      end function transition_func
   end interface

   class(ErosionTransition), allocatable :: erosion_transition

contains

   ! initialize ErosionTransition with shared parameters
   module subroutine ErosionTransition_constructor(this, params)
      class(ErosionTransition), intent(inout) :: this
      type(Dict), intent(in) :: params
      
      if (params%has_key("erosion depth")) then
         call params%get("erosion depth", this%erosion_depth)
      else
         call FatalErrorMessage("Required parameter Erosion Depth not given.")
      end if

    end subroutine ErosionTransition_constructor

   module subroutine ErosionTransition_validate(this)
      class(ErosionTransition), intent(in) :: this
      if (this%erosion_depth < 0) then
         call FatalErrorMessage("Erosion depth must be non-negative.  Received ", this%erosion_depth)
      end if
   end subroutine ErosionTransition_validate   

end module erosion_transition_interface