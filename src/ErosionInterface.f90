module erosion_interface

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use messages_module, only: FatalErrorMessage, WarningMessage
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: GeometricCorrectionFactor, FlowSquaredSpeedSlopeAligned
   use drag_interface, only: DragModel
   use erosion_transition_interface, only: ErosionTransition
   use erosion_transition_types

   implicit none

   private
   public :: ErosionModel, erosion_model

   type, abstract :: ErosionModel
      type(varString) :: name
      real(kind=wp) :: erosion_rate = 0.001_wp ! Global default
      real(kind=wp) :: solid_diameter
      real(kind=wp) :: g, rhow, rhos, gred
      real(kind=wp) :: visc_w = 1.2e-6_wp
      type(Dict) :: parameters
      class(ErosionTransition), allocatable :: erosion_transition
   contains
      procedure :: initialize => ErosionModel_constructor
      procedure, non_overridable :: base_validate => ErosionModel_validate
      procedure, non_overridable :: shields_number => ErosionModel_FlowShieldsNumber
      procedure, non_overridable :: particle_speed => ErosionModel_FlowParticleSpeed
      procedure(erosion_validate), deferred :: validate
      procedure(erosion_func), deferred :: erosion
   end type ErosionModel

   abstract interface
      subroutine erosion_validate(this)
         import :: ErosionModel
         class(ErosionModel), intent(in) :: this
      end subroutine erosion_validate
   end interface

   abstract interface
      pure function erosion_func(this, drag_model, RunParams, uvect) result(erosion)
         import :: ErosionModel
         import :: DragModel
         import :: wp
         import :: RunSet
         class(ErosionModel), intent(in) :: this
         class(DragModel), intent(in) :: drag_model
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: erosion
      end function erosion_func
   end interface

   class(ErosionModel), allocatable :: erosion_model

contains

   ! initialize ErosionModel with shared parameters
   subroutine ErosionModel_constructor(this, params)
      class(ErosionModel), intent(inout) :: this
      type(Dict), intent(in) :: params

      type(varString) :: erosion_transition

      call params%get("g", this%g)
      call params%get("rhow", this%rhow)
      call params%get("rhos", this%rhos)
      call params%get("rhow", this%rhow)
      call params%get("erosion rate", this%erosion_rate)
      call params%get("solid diameter", this%solid_diameter)

      this%gred = (this%rhos/this%rhow-1.0_wp)*this%g

      if (params%has_key("erosion transition")) then
         call params%get("erosion transition", erosion_transition)

         select case (erosion_transition%s)
            case ('smooth')
            !    allocate(this%erosion_transition, source=SmoothErosionTransition(params))
               this%erosion_transition = SmoothErosionTransition(params)
            case ('step')
               this%erosion_transition = StepErosionTransition(params)
            case ('none')
               this%erosion_transition = NoErosionTransition(params)
            case default
               this%erosion_transition = SmoothErosionTransition(params) 
               call WarningMessage('Erosion Transition not recognized.  Using default Smooth')
         end select
      else
         this%erosion_transition = SmoothErosionTransition(params)
         call WarningMessage('Erosion Transition not given.  Using default Smooth')
      end if

      call this%parameters%append("g", this%g)
      call this%parameters%append("rhow", this%rhow)
      call this%parameters%append("rhos", this%rhos)
      call this%parameters%append("rhow", this%rhow)
      call this%parameters%append("erosion rate", this%erosion_rate)
      call this%parameters%append("solid diameter", this%solid_diameter)
      call this%parameters%append("erosion transition", this%erosion_transition%name)

   end subroutine ErosionModel_constructor

   subroutine ErosionModel_validate(this)
      class(ErosionModel), intent(in) :: this

      if (this%g <=0) call FatalErrorMessage('g must be positive.  Received ', this%g)
      if (this%rhow <=0) call FatalErrorMessage('rhow must be positive.  Received ', this%rhow)
      if (this%rhos <=0) call FatalErrorMessage('rhos must be positive.  Received ', this%rhos)
      if (this%erosion_rate <=0) call FatalErrorMessage('Erosion Rate must be positive.  Received ', this%erosion_rate)
      if (this%solid_diameter <=0) call FatalErrorMessage('Solid Diameter must be positive.  Received ', this%solid_diameter)

      if (.not. allocated(this%erosion_transition)) call FatalErrorMessage('Erosion Transition must be specified for ' // this%name%s)
      
   end subroutine ErosionModel_validate

   ! Shields number for a generic drag_model.
   pure function ErosionModel_FlowShieldsNumber(this, drag_model, RunParams, uvect) result(shields)
      class(ErosionModel), intent(in) :: this
      class(DragModel), intent(in) :: drag_model
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: shields

      real(kind=wp) :: gred, modu2, friction

      ! g' -> g'_\perp = g' cos(theta)
      gred = this%gred / GeometricCorrectionFactor(RunParams, uvect)

      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)

      friction = drag_model%friction(RunParams, uvect)

      shields = friction / (gred * erosion_model%solid_diameter)
   end function ErosionModel_FlowShieldsNumber

   pure function ErosionModel_FlowParticleSpeed(this, RunParams, uvect) result(u_p)
      class(ErosionModel), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: u_p

      real(kind=wp) :: gred

      ! g' -> g'_\perp = g' cos(theta)
      gred = this%gred / GeometricCorrectionFactor(RunParams, uvect)

      u_p = sqrt(gred * this%solid_diameter)

   end function ErosionModel_FlowParticleSpeed

end module erosion_interface