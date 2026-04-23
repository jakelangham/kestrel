submodule (drag_types) drag_coulomb

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: GeometricCorrectionFactor

   implicit none

   ! -- Drag closure: Coulomb drag --
   !
   !    Coulomb drag has friction function F = \mu * g * Hn.
   !
   !    Parameters:
   !        coulomb_co > 0
   !
   !    Global parameters (in RunSet):
   !        g
   !
contains

   ! Constructor
   module function coulomb_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(CoulombDrag) :: this
      this%name = varString("Coulomb")
      if (params%has_key('coulomb co')) then
         call params%get('coulomb co', this%coulomb_co)
      else
         call WarningMessage('Coulomb Co not given.  Using default value ', this%coulomb_co)
      end if

      ! Parameters
      call this%parameters%append("coulomb co", this%coulomb_co)
   end function coulomb_constructor

   ! Validator
   module subroutine coulomb_validate(this)
      class(CoulombDrag), intent(in) :: this
      if (this%coulomb_co<0) then
         call FatalErrorMessage('Coulomb Co must be positive.  Received', this%coulomb_co)
      end if
   end subroutine coulomb_validate

   ! Friction
   module pure function coulomb_friction(this, RunParams, uvect) result(friction)
      class(CoulombDrag), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: Hn
      real(kind=wp) :: gam
      real(kind=wp) :: g

      Hn = uvect(RunParams%Vars%Hn)
      gam = GeometricCorrectionFactor(RunParams, uvect)
      g = RunParams%g / gam

      friction = this%coulomb_co * g * Hn
      
   end function coulomb_friction

end submodule drag_coulomb