submodule (drag_types) drag_chezy

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: FlowSquaredSpeedSlopeAligned

   implicit none

   ! -- Drag closure: Chezy drag --
   !
   !    Chezy drag has friction function F = Cd * |u|^2.
   !
   !    Parameters:
   !        chezy_drag > 0
   !
contains

   ! Constructor
   module function chezy_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(ChezyDrag) :: this
      this%name = varString("Chezy")
      
      if (params%has_key('chezy co')) then
         call params%get('chezy co', this%chezy_co)
      else
         call WarningMessage('Chezy Co not given.  Using default value ', this%chezy_co)
      end if

      ! Parameters
      call this%parameters%append("chezy co", this%chezy_co)
   end function chezy_constructor

   ! Validator
   module subroutine chezy_validate(this)
      class(ChezyDrag), intent(in) :: this
      if (this%chezy_co<0) then
         call FatalErrorMessage('Chezy Co must be positive.  Received', this%chezy_co)
      end if
   end subroutine chezy_validate

   ! Friction
   module pure function chezy_friction(this, RunParams, uvect) result(friction)
      class(ChezyDrag), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: modu2

      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)
      friction = this%chezy_co * modu2

   end function chezy_friction

end submodule drag_chezy