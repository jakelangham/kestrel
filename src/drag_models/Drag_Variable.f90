submodule (drag_types) drag_variable

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: FlowSquaredSpeedSlopeAligned, fSwitch, set_switcher
   use drag_interface, only: DragModel

   implicit none

   ! -- Drag closure: Variable drag --
   !
   !    Variable drag is a concentration-dependent
   !    variation between two DragModels: 
   !       fluid_drag for dilute flows
   !       granular_drag for dense flow.
   !    The transition is specified by a switching function fSwitch.
   !
   !    Parameters:
   !        fluid_drag
   !        granular_drag
   !
contains

   ! Constructor
   module function variable_constructor(params, fluid_drag, granular_drag) result(this)
      type(Dict), intent(in) :: params   
      class(DragModel), intent(in) :: fluid_drag
      class(DragModel), intent(in) :: granular_drag
      type(VariableDrag) :: this
      this%name = varString("Variable")
      if (params%has_key('voellmy switch rate')) then
         call params%get('voellmy switch rate', this%voellmy_switch_rate)
      else
         call WarningMessage('Voellmy Switch Rate not given.  Using default value ', this%voellmy_switch_rate)
      end if
      if (params%has_key('voellmy switch value')) then
         call params%get('voellmy switch value', this%voellmy_switch_value)
      else
         call WarningMessage('Voellmy Switch Value not given.  Using default value ', this%voellmy_switch_value)
      end if

      if (params%has_key('switch function')) then
         call params%get('switch function', this%switch_function)
      else
         call WarningMessage('switch function not given.  Using default value tanh')
         this%switch_function = varString('tanh')
      end if
      call set_switcher(this%switch_function%s, this%switcher)
      
      this%fluid = fluid_drag
      
      this%granular = granular_drag

      ! Parameters
      this%parameters = this%fluid%parameters + this%granular%parameters
      call this%parameters%append("switch function", this%switch_function)
      call this%parameters%append("voellmy switch rate", this%voellmy_switch_rate)
      call this%parameters%append("voellmy switch value", this%voellmy_switch_value)
   end function variable_constructor

   ! Validator
   module subroutine variable_validate(this)
      class(VariableDrag), intent(in) :: this

      real(kind=wp) :: f0
      type(varString) :: msg

      call this%fluid%validate()
      call this%granular%validate()

      if (this%voellmy_switch_rate<0) then
         call FatalErrorMessage('Voellmy Switch Rate must be positive.  Received', this%voellmy_switch_rate)
      end if
      if (this%voellmy_switch_value<0) then
         call FatalErrorMessage('Voellmy Switch Value must be positive.  Received', this%voellmy_switch_value)
      end if

      f0 = this%switcher(this%voellmy_switch_rate, this%voellmy_switch_value, 0.0_wp)
      if ( f0 > 0.0_wp) then
        msg = varString('Switch function has f(0) = ')
        msg = msg + varString(f0)
        msg = msg%trim() + ' > 0'
        call WarningMessage(msg%s)
      end if
   end subroutine variable_validate

   ! Friction
   module pure function variable_friction(this, RunParams, uvect) result(friction)
      class(VariableDrag), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: psi
      real(kind=wp) :: fluid_friction
      real(kind=wp) :: granular_friction
      real(kind=wp) :: fc

      psi = uvect(RunParams%Vars%psi)
      
      fluid_friction = this%fluid%friction(RunParams, uvect)
      granular_friction = this%granular%friction(RunParams, uvect)
      
      fc = this%switcher(this%voellmy_switch_rate, this%voellmy_switch_value, psi)

      friction = fluid_friction * (1.0_wp - fc) + granular_friction * fc
      
   end function variable_friction

end submodule drag_variable