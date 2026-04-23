submodule (erosion_transition_types) erosion_no_transition

   use set_precision_module, only: wp
   use dict_module, only: Dict

   implicit none

contains

   module function NoErosionTransition_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(NoErosionTransition) :: this
      call this%initialize(params)
      this%name = varString('no transition')
   end function NoErosionTransition_constructor

   module subroutine NoErosionTransition_validate(this)
      class(NoErosionTransition), intent(in) :: this
      call this%base_validate()
   end subroutine NoErosionTransition_validate

   module pure function NoErosionTransition_transition(this, Bt) result(factor)
      class(NoErosionTransition), intent(in) :: this
      real(kind=wp), intent(in) :: Bt
      real(kind=wp) :: factor

      factor = 1.0_wp
   end function NoErosionTransition_transition

end submodule erosion_no_transition