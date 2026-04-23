submodule (erosion_types) simple_erosion

! SimpleErosion
! A fluid (Shields number) based erosion *without* a critical Shields number.
! The erosion rate is given by E = up * eps * S where up is the particle
! speed scale (= sqrt(g' d)) eps is the user defined erosion rate S is the
! Shields number.

   use set_precision_module, only: wp
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use drag_interface, only: DragModel

   implicit none

contains

   module function SimpleErosion_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(SimpleErosion) :: this
      this%name = varString("Simple erosion")
      call this%initialize(params)
   end function SimpleErosion_constructor

   module subroutine SimpleErosion_validate(this)
      class(SimpleErosion), intent(in) :: this
      call this%base_validate()
   end subroutine SimpleErosion_validate

   module pure function SimpleErosion_erosion(this, drag_model, RunParams, uvect) result(ero)
      implicit none
      class(SimpleErosion), intent(in) :: this
      class(DragModel), intent(in) :: drag_model
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      real(kind=wp) :: shields, Bt

      shields = this%shields_number(drag_model, RunParams, uvect)

      ero = this%erosion_rate * shields
      ero = ero * this%particle_speed(RunParams, uvect)

      Bt = uvect(RunParams%Vars%bt)

      ero = ero * this%erosion_transition%transition(Bt)

   end function SimpleErosion_erosion

end submodule simple_erosion