submodule (drag_interface) drag_shared

   use set_precision_module, only: wp
   use messages_module, only: WarningMessage, FatalErrorMessage
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: FlowSquaredSpeedSlopeAligned
   use drag_interface

   implicit none

   ! -- Drag closures shared implementations --
   !
   !
   !    The 'validate' function does validation of model parameters
   !
   !    Because of the way the time stepper works (treating the drag
   !    implicitly), these are required to define a 'friction' function F such that
   !    the basal stress is \tau_x = -F \rho u/|(u,v)|, \tau_y = -F \rho v/|(u,v)|
   !
   !    The basal stress are given by the 'drag_xy' function.
   !    It has a shared implementation for all models
   
contains

   module pure subroutine drag_xy(this, RunParams, uvect, drag_x, drag_y)
      class(DragModel), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp), intent(out) :: drag_x, drag_y
      real(kind=wp) :: rho, u, v, modu
      
      modu = sqrt(FlowSquaredSpeedSlopeAligned(RunParams, uvect))
      rho = uvect(RunParams%Vars%rho)
      u = uvect(RunParams%Vars%u)
      drag_x = -rho * this%friction(RunParams, uvect) * u/modu
      if (RunParams%isOneD) then
         drag_y = 0.0_wp
      else
         v = uvect(RunParams%Vars%v)
         drag_y = -rho * this%friction(RunParams, uvect) * v/modu
      end if
      
   end subroutine drag_xy

end submodule drag_shared