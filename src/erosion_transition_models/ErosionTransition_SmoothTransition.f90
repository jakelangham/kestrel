submodule (erosion_transition_types) erosion_smooth_transition

   use set_precision_module, only: wp
   use dict_module, only: Dict

   implicit none

contains

   module function SmoothErosionTransition_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(SmoothErosionTransition) :: this
      call this%initialize(params)
      this%name = varString('smooth transition')
      call params%get_or_default('smooth transition rate', this%rate, default=1e5_wp)
   end function SmoothErosionTransition_constructor

   module subroutine SmoothErosionTransition_validate(this)
      class(SmoothErosionTransition), intent(in) :: this
      call this%base_validate()
      if (this%rate <= 0) then
         call FatalErrorMessage('Smooth Transition Rate must be positive.  Received ', this%rate)
      end if
   end subroutine SmoothErosionTransition_validate

   module pure function SmoothErosionTransition_transition(this, Bt) result(factor)
      class(SmoothErosionTransition), intent(in) :: this
      real(kind=wp), intent(in) :: Bt
      real(kind=wp) :: factor

      factor = 0.5_wp * (1.0_wp + tanh(1e5_wp * (Bt + this%erosion_depth)))
   end function SmoothErosionTransition_transition

end submodule erosion_smooth_transition