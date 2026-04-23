submodule (erosion_types) mixed_erosion

! MixedErosion
! Solid concentration weighted combination of fluid and granular erosion.
! The weighting is determined using the switching function.

   use set_precision_module, only: wp
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use drag_interface, only: DragModel
   use closures_module, only: fSwitch, set_switcher

   implicit none

contains

   module function MixedErosion_constructor(params, fluid_erosion, granular_erosion) result(this)
      type(Dict), intent(in) :: params
      class(ErosionModel), intent(in) :: fluid_erosion
      class(ErosionModel), intent(in) :: granular_erosion
      type(MixedErosion) :: this

      this%name = varString("Mixed erosion")
      call this%initialize(params)

      this%fluid_erosion = fluid_erosion
      this%granular_erosion = granular_erosion

      if (params%has_key('voellmy switch rate')) then
         call params%get('voellmy switch rate', this%voellmy_switch_rate)
      end if
      if (params%has_key('voellmy switch value')) then
         call params%get('voellmy switch value', this%voellmy_switch_value)
      end if

      if (params%has_key('switch function')) then
         call params%get('switch function', this%switch_function)
      else
         call WarningMessage('switch function not given.  Using default value tanh')
         this%switch_function = varString('tanh')
      end if
      call set_switcher(this%switch_function%s, this%switcher)

      this%parameters = this%parameters + this%fluid_erosion%parameters
      this%parameters = this%parameters + this%granular_erosion%parameters

      call this%parameters%append('voellmy switch rate', this%voellmy_switch_rate)
      call this%parameters%append('voellmy switch value', this%voellmy_switch_value)
      call this%parameters%append('switch function', this%switch_function)
      
   end function MixedErosion_constructor

   module subroutine MixedErosion_validate(this)
      class(MixedErosion), intent(in) :: this
      call this%base_validate()
      call this%fluid_erosion%validate()
      call this%granular_erosion%validate()

      if (this%voellmy_switch_rate<0) then
         call FatalErrorMessage('Voellmy Switch Rate must be positive.  Received', this%voellmy_switch_rate)
      end if
      if (this%voellmy_switch_value<0) then
         call FatalErrorMessage('Voellmy Switch Value must be positive.  Received', this%voellmy_switch_value)
      end if
   end subroutine MixedErosion_validate

   module pure function MixedErosion_erosion(this, drag_model, RunParams, uvect) result(ero)
      class(MixedErosion), intent(in) :: this
      class(DragModel), intent(in) :: drag_model
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      real(kind=wp) :: psi
      real(kind=wp) :: fc
      real(kind=wp) :: fluid_ero
      real(kind=wp) :: granular_ero

      psi = uvect(RunParams%Vars%psi)

      ! Get the weighting.
      fc = this%switcher(this%voellmy_switch_rate, this%voellmy_switch_value, psi)

      fluid_ero = this%fluid_erosion%erosion(drag_model, RunParams, uvect)
      granular_ero = this%granular_erosion%erosion(drag_model, RunParams, uvect)

      ero = (1.0_wp - fc) * fluid_ero + fc * granular_ero

   end function MixedErosion_erosion

end submodule mixed_erosion