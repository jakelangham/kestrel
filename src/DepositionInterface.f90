module deposition_interface

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use messages_module, only: FatalErrorMessage, WarningMessage
   use dict_module, only: Dict
   use runsettings_module, only: RunSet

   implicit none

   private
   public :: DepositionModel, deposition_model

   ! -- Deposition closures --
   !
   !    Abstract type and abstract interface for deposition models.
   !    Specific closures are given in the /deposition_models/
   !
   !    The 'validate' function does validation of model parameters
   !
   !    The 'flux' function computes the deposition flux
   !
   type, abstract :: DepositionModel
      type(varString) :: name
      real(kind=wp) :: ws0
      real(kind=wp) :: maxPack
      real(kind=wp) :: solid_diameter
      real(kind=wp) :: g, rhow, rhos, gred
      real(kind=wp) :: visc_w = 1.2e-6_wp
      type(Dict) :: parameters
   contains
      procedure, non_overridable :: initialize => DepositionModel_constructor
      procedure, non_overridable :: base_validate => DepositionModel_validate
      procedure(hindering_func), deferred :: hindering
      procedure(deposition_validate), deferred :: validate
      procedure, non_overridable :: rate => Deposition_rate
   end type DepositionModel

   abstract interface
      subroutine deposition_validate(this)
         import :: DepositionModel
         class(DepositionModel), intent(in) :: this
      end subroutine deposition_validate
   end interface

   abstract interface
      pure function hindering_func(this, psi) result(deposition)
         import :: DepositionModel
         import :: wp
         class(DepositionModel), intent(in) :: this
         real(kind=wp), intent(in) :: psi
         real(kind=wp) :: deposition
      end function hindering_func
   end interface

   class(DepositionModel), allocatable :: deposition_model

contains

    subroutine DepositionModel_constructor(this, params, RunParams)
      ! initialize ErosionModel with shared parameters
      class(DepositionModel), intent(inout) :: this
      type(Dict), intent(in) :: params
      type(RunSet), intent(in) :: RunParams

      real(kind=wp) :: R
      real(kind=wp) :: ws0_d

      call params%get_or_default("g", this%g, RunParams%g)
      call params%get_or_default("rhow", this%rhow, RunParams%rhow)
      call params%get_or_default("rhos", this%rhos, RunParams%rhos)
      call params%get("water viscosity", this%visc_w)
      call params%get("solid diameter", this%solid_diameter)
      call params%get_or_default("maxPack", this%maxPack, RunParams%maxPack)

      this%gred = (this%rhos/this%rhow-1.0_wp)*this%g

      R = (this%gred/this%visc_w/this%visc_w)**(1.0_wp/3.0_wp)*this%solid_diameter ! Scaled grain size

      ws0_d = this%visc_w/this%solid_diameter*(sqrt(10.36_wp*10.36_wp + 1.048_wp*R*R*R) - 10.36_wp) ! Clear water settling velocity

      call params%get_or_default('settling speed', this%ws0, ws0_d)

      call this%parameters%append("g", this%g)
      call this%parameters%append("rhow", this%rhow)
      call this%parameters%append("rhos", this%rhos)
      call this%parameters%append("water viscosity", this%visc_w)
      call this%parameters%append("solid diameter", this%solid_diameter)
      call this%parameters%append("maxPack", this%maxPack)
      call this%parameters%append('settling speed', this%ws0)

   end subroutine DepositionModel_constructor

   subroutine DepositionModel_validate(this)
      class(DepositionModel), intent(in) :: this

      if (this%g <=0) call FatalErrorMessage('g must be positive.  Received ', this%g)
      if (this%rhow <=0) call FatalErrorMessage('rhow must be positive.  Received ', this%rhow)
      if (this%rhos <=0) call FatalErrorMessage('rhos must be positive.  Received ', this%rhos)
      if (this%visc_w <=0) call FatalErrorMessage('Water viscosity must be positive.  Received ', this%visc_w)
      if (this%ws0 <=0) call FatalErrorMessage('Speedling speed must be positive.  Received ', this%ws0)
      if (this%solid_diameter <=0) call FatalErrorMessage('Solid Diameter must be positive.  Received ', this%solid_diameter)
      if (this%maxPack <= 0) call FatalErrorMessage('maxPack must be positive.  Received ', this%maxPack)
      if (this%maxPack > 1) call FatalErrorMessage('maxPack must be less than 1.  Received ', this%maxPack)

   end subroutine DepositionModel_validate

   ! Given a solution vector uvect, compute the corresponding deposition rate,
   ! returning in dep.
   pure function Deposition_rate(this, RunParams, uvect) result(rate)
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

end module deposition_interface