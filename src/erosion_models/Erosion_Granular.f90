submodule (erosion_types) granular_erosion

! FluidErosion
! A fluid (Shields number) based erosion *without* a critical Shields number.
! The erosion rate is given by E = up * eps * S where up is the particle
! speed scale (= sqrt(g' d)) eps is the user defined erosion rate S is the
! Shields number.

   use set_precision_module, only: wp, pi
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: FlowSquaredSpeedSlopeAligned, GeometricCorrectionFactor
   use drag_interface, only: DragModel
   use drag_types, only: PouliquenDrag, Edwards2019Drag, VariableDrag

   implicit none

contains

   module function GranularErosion_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(GranularErosion) :: this

      this%name = varString("Granular erosion")
      call this%initialize(params)

      call params%get('granular erosion rate', this%granular_erosion_rate)
      call params%get('pouliquen min', this%pouliquen_min)
      call params%get('pouliquen max', this%pouliquen_max)
      call params%get('pouliquen beta', this%pouliquen_beta)

      this%t1 = tan(pi / 180.0_wp) ! t1 = tan(1 deg)

      this%pouliquen_static_slope = (this%pouliquen_min + this%t1) / (1.0_wp - this%pouliquen_min * this%t1)

      call this%parameters%append('granular erosion rate', this%granular_erosion_rate)
      call this%parameters%append('pouliquen min', this%pouliquen_min)
      call this%parameters%append('pouliquen max', this%pouliquen_max)
      call this%parameters%append('pouliquen beta', this%pouliquen_beta)
      
   end function GranularErosion_constructor

   module subroutine GranularErosion_validate(this)
      class(GranularErosion), intent(in) :: this
      call this%base_validate()
      if (this%granular_erosion_rate < 0) call FatalErrorMessage('Granular Erosion Rate must be non-negative.  Received ', this%granular_erosion_rate)
      if (this%pouliquen_min <= 0) call FatalErrorMessage('Pouliquen Min must be positive.  Received ', this%pouliquen_min)
      if (this%pouliquen_max <= 0) call FatalErrorMessage('Pouliquen Max must be positive.  Received ', this%pouliquen_max)
      if (this%pouliquen_beta <= 0) call FatalErrorMessage('Pouliquen Beta must be positive.  Received ', this%pouliquen_beta)
   end subroutine GranularErosion_validate

   module pure function GranularErosion_erosion(this, drag_model, RunParams, uvect) result(ero)
      implicit none
      class(GranularErosion), intent(in) :: this
      class(DragModel), intent(in) :: drag_model
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: ero

      real(kind=wp) :: Hn, modu2
      real(kind=wp) :: gcostheta, mu, muNeutral

      class(DragModel), allocatable :: granular_drag
      
      Hn = uvect(RunParams%Vars%Hn)
      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)

      ! granular erosion, eps_g*(mu - mu_n)
      gcostheta = this%g / GeometricCorrectionFactor(RunParams, uvect)
      select type(drag_model)
         type is (PouliquenDrag)
            mu = drag_model%friction_coefficient(RunParams%heightThreshold, gcostheta, Hn, sqrt(modu2))
         type is (Edwards2019Drag)
            mu = drag_model%friction_coefficient(RunParams%heightThreshold, gcostheta, Hn, sqrt(modu2))
         type is (VariableDrag)
            granular_drag = drag_model%granular
            select type (granular_drag)
               type is (PouliquenDrag)
                  mu = granular_drag%friction_coefficient(RunParams%heightThreshold, gcostheta, Hn, sqrt(modu2))
               type is (Edwards2019Drag)
                  mu = granular_drag%friction_coefficient(RunParams%heightThreshold, gcostheta, Hn, sqrt(modu2))
            end select
         class default
            error stop "GranularErosion requires either Pouliquen or Edwards2019 drag.  Received " // drag_model%name%s
      end select
      muNeutral = this%pouliquen_min + (this%pouliquen_static_slope - this%pouliquen_min) /  &
         (1.0_wp + (Hn / 25.0_wp / this%solid_diameter)**2.0_wp)
      if (mu > muNeutral) then
         ero = this%granular_erosion_rate * (mu - muNeutral)
         ero = ero* this%particle_speed(RunParams, uvect)

         ero = ero * this%erosion_transition%transition(uvect(RunParams%Vars%Bt))
      else
         ero = 0.0_wp
      end if

   end function GranularErosion_erosion

end submodule granular_erosion