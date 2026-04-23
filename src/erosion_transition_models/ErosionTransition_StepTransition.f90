submodule (erosion_transition_types) erosion_step_transition

   use set_precision_module, only: wp
   use dict_module, only: Dict

   implicit none

contains

   module function StepErosionTransition_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(StepErosionTransition) :: this
      call this%initialize(params)
      this%name = varString('step transition')
   end function StepErosionTransition_constructor

   module subroutine StepErosionTransition_validate(this)
      class(StepErosionTransition), intent(in) :: this
      call this%base_validate()
   end subroutine StepErosionTransition_validate

   module pure function StepErosionTransition_transition(this, Bt) result(factor)
      class(StepErosionTransition), intent(in) :: this
      real(kind=wp), intent(in) :: Bt
      real(kind=wp) :: factor

      if (Bt < -this%erosion_depth) then
        factor = 0.0_wp
      else
        factor = 1.0_wp
      end if
   end function StepErosionTransition_transition

end submodule erosion_step_transition