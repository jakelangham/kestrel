submodule (drag_types) drag_voellmy

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: FlowSquaredSpeedSlopeAligned, fSwitch
   use drag_types, only: ChezyDrag, CoulombDrag

   implicit none

   ! -- Drag closure: Voellmy drag --
   !
   !    Voellmy drag is the sum of Chezy and Coulomb drag.
   !
   !    Parameters:
   !        chezy_part
   !        coulomb_part
   !
   ! Voellmy drag is the sum of Chezy and Coulomb drag.
contains

   ! Constructor
   module function voellmy_constructor(chezy_part, coulomb_part) result(this)
      type(ChezyDrag), intent(in) :: chezy_part
      type(CoulombDrag), intent(in) :: coulomb_part
      type(VoellmyDrag) :: this
      this%name = varString("Voellmy")

      this%chezy_part = chezy_part
      
      this%coulomb_part = coulomb_part

      ! Parameters
      this%parameters = this%chezy_part%parameters + this%coulomb_part%parameters
      
   end function voellmy_constructor

   ! Validator
   module subroutine voellmy_validate(this)
      class(VoellmyDrag), intent(in) :: this

      call this%chezy_part%validate()
      call this%coulomb_part%validate()

   end subroutine voellmy_validate

   ! Friction
   module pure function voellmy_friction(this, RunParams, uvect) result(friction)
      class(VoellmyDrag), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      friction = this%chezy_part%friction(RunParams, uvect) + this%coulomb_part%friction(RunParams, uvect)            
      
   end function voellmy_friction

end submodule drag_voellmy