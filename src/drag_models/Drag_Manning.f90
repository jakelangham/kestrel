submodule (drag_types) drag_manning

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: GeometricCorrectionFactor, Reciprocal

   implicit none

   ! -- Drag closure: Manning drag --
   !
   !    Manning drag has friction function F = g * n^2 / (Hn^{1/3}) 
   !    where n is the Manning coefficient.
   !
   !    Parameters:
   !        n > 0
   !
   !    Global parameters (in RunSet):
   !        g

contains

   ! Constructor
   module function manning_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(ManningDrag) :: this
      this%name = varString("Manning")
      if (params%has_key('manning co')) then
         call params%get('manning co', this%manning_co)
      else
         call WarningMessage('Manning Co not given.  Using default value ', this%manning_co)
      end if

      ! Parameters
      call this%parameters%append("manning co", this%manning_co)
   end function manning_constructor

   ! Validator
   module subroutine manning_validate(this)
      class(ManningDrag), intent(in) :: this
      if (this%manning_co<0) then
         call FatalErrorMessage('Manning Co must be positive.  Received', this%manning_co)
      end if
   end subroutine manning_validate

   ! Friction
   module pure function manning_friction(this, RunParams, uvect) result(friction)
      class(ManningDrag), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: Hn, Hn_recip
      real(kind=wp) :: gam
      real(kind=wp) :: g
      real(kind=wp) :: ManningCo

        Hn = uvect(RunParams%Vars%Hn)
        gam = GeometricCorrectionFactor(RunParams, uvect)
        g = RunParams%g / gam

        ManningCo = this%manning_co
      
        Hn_recip = Reciprocal(Hn, RunParams%heightThreshold)

        friction = g * ManningCo*ManningCo * (Hn_recip**(1.0_wp/3.0_wp))
        
   end function manning_friction

end submodule drag_manning