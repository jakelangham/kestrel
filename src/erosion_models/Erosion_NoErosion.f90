submodule (erosion_types) no_erosion

! NoErosion
! Set erosion rate to zero.

   use set_precision_module, only: wp
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use drag_interface, only: DragModel

   implicit none

contains

   module function NoErosion_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(NoErosion) :: this
      this%name = varString("No erosion")
      call this%initialize(params)
   end function NoErosion_constructor

   module subroutine NoErosion_validate(this)
      class(NoErosion), intent(in) :: this
      call this%base_validate()
   end subroutine NoErosion_validate

   module pure function NoErosion_erosion(this, drag_model, RunParams, uvect) result(ero)
      class(NoErosion), intent(in) :: this
      class(DragModel), intent(in) :: drag_model
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      ero = 0.0_wp

   end function NoErosion_erosion

end submodule no_erosion