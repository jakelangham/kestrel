submodule (deposition_interface) deposition_shared

   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage, WarningMessage
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: GeometricCorrectionFactor, FlowSquaredSpeedSlopeAligned
   use erosion_interface, only: ErosionModel
   use erosion_transition_types, only: SmoothErosionTransition, StepErosionTransition, NoErosionTransition

   implicit none

   ! -- Deposition closures shared implementations --
   !
contains

   module subroutine DepositionModel_constructor(this, params)
      ! initialize ErosionModel with shared parameters
      class(DepositionModel), intent(inout) :: this
      type(Dict), intent(in) :: params

      real(kind=wp) :: R
      real(kind=wp) :: ws0_d

      call params%get("g", this%g)
      call params%get("rhow", this%rhow)
      call params%get("rhos", this%rhos)
      call params%get("water viscosity", this%visc_w)
      call params%get("solid diameter", this%solid_diameter)
      call params%get("maxPack", this%maxPack)

      this%gred = (this%rhos/this%rhow-1.0_wp)*this%g

      R = (this%gred/this%visc_w/this%visc_w)**(1.0_wp/3.0_wp)*this%solid_diameter ! Scaled grain size

      ws0_d = this%visc_w/this%solid_diameter*(sqrt(10.36_wp*10.36_wp + 1.048_wp*R*R*R) - 10.36_wp) ! Clear water settling velocity

      call params%get_or_default('settling speed', this%ws0, ws0_d)

   end subroutine DepositionModel_constructor

   module subroutine DepositionModel_validate(this)
      class(DepositionModel), intent(in) :: this

      if (this%g <=0) call FatalErrorMessage('g must be positive.  Received ', this%g)
      if (this%rhow <=0) call FatalErrorMessage('rhow must be positive.  Received ', this%rhow)
      if (this%rhos <=0) call FatalErrorMessage('rhos must be positive.  Received ', this%rhos)
      if (this%visc_w <=0) call FatalErrorMessage('Water viscosity must be positive.  Received ', this%visc_w)
      if (this%ws0 <=0) call FatalErrorMessage('Speedling speed must be positive.  Received ', this%ws0)
      if (this%solid_diameter <=0) call FatalErrorMessage('Solid Diameter must be positive.  Received ', this%solid_diameter)
      if (this%maxPack <= 0) call FatalErrorMessage('maxPack must be positive.  Received ', this%maxPack)
      if (this%maxPack > 1) call FatalErrorMessage('maxPack must be less than 1.  Received ', this%maxPack)

      call this%validate()
   end subroutine DepositionModel_validate

   ! Given a solution vector uvect, compute the corresponding deposition rate,
   ! returning in dep.
   module pure function Deposition_rate(this, RunParams, uvect) result(rate)
      implicit none
      class(DepositionModel), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: rate

      real(kind=wp) :: psi, alpha

      psi = uvect(RunParams%Vars%psi)

      alpha = 0.0_wp
      if (psi >= this%maxPack) then
         alpha = 0.0_wp
      else if (psi > 0.0_wp) then
         alpha = deposition_model%hindering(psi)
      else
         alpha = 0.0_wp
      end if

      rate = this%ws0 * alpha
   end function Deposition_rate

end submodule deposition_shared