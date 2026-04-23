submodule (erosion_types) fluid_erosion

! FluidErosion
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

   module function FluidErosion_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(FluidErosion) :: this
      
      real(kind=wp) :: R

      this%name = varString("Fluid erosion")
      call this%initialize(params)

      R = (this%gred/this%visc_w/this%visc_w)**(1.0_wp/3.0_wp)*this%solid_diameter ! Scaled grain size

      this%critical_shields = 0.3_wp/(1.0_wp+1.2_wp*R) + 0.055_wp*(1.0_wp - exp(-0.02_wp*R))

   end function FluidErosion_constructor

   module subroutine FluidErosion_validate(this)
      class(FluidErosion), intent(in) :: this
      call this%base_validate()
   end subroutine FluidErosion_validate

   module pure function FluidErosion_erosion(this, drag_model, RunParams, uvect) result(ero)
      class(FluidErosion), intent(in) :: this
      class(DragModel), intent(in) :: drag_model
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      real(kind=wp) :: shields

      shields = this%shields_number(drag_model, RunParams, uvect)

      if (shields > this%critical_shields) then
         ero = this%erosion_rate * (shields - this%critical_shields)
         ero = ero * this%particle_speed(RunParams, uvect) * this%erosion_transition%transition(uvect(RunParams%Vars%bt))
      else
         ero = 0.0_wp
      end if
   end function FluidErosion_erosion

end submodule fluid_erosion