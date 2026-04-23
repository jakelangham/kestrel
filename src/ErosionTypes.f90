module erosion_types

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use messages_module, only: FatalErrorMessage, WarningMessage
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use drag_interface, only: DragModel
   use erosion_interface, only: ErosionModel
   use closures_module, only: fSwitch

   implicit none

   private
   public :: SimpleErosion, FluidErosion, GranularErosion, MixedErosion, NoErosion
   public :: create_erosion_model

!--- SimpleErosion
!--- A fluid (Shields number) based erosion *without* a critical Shields number.
!--- The erosion rate is given by E = up * eps * S where up is the particle
!--- speed scale (= sqrt(g' d)) eps is the user defined erosion rate S is the
!--- Shields number.
   type, extends(ErosionModel) :: SimpleErosion
      contains
         procedure :: validate => SimpleErosion_validate
         procedure :: erosion => SimpleErosion_erosion
   end type SimpleErosion

   interface SimpleErosion
      procedure SimpleErosion_constructor
   end interface SimpleErosion

   interface
      module function SimpleErosion_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(SimpleErosion) :: this
      end function SimpleErosion_constructor

      module subroutine SimpleErosion_validate(this)
         class(SimpleErosion), intent(in) :: this
      end subroutine SimpleErosion_validate

      module pure function SimpleErosion_erosion(this, drag_model, RunParams, uvect) result(erosion)
         class(SimpleErosion), intent(in) :: this
         class(DragModel), intent(in) :: drag_model
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: erosion
      end function SimpleErosion_erosion
   end interface

!--- Fluid erosion
!--- Shields stress erosion rate
   type, extends(ErosionModel) :: FluidErosion
      real(kind=wp) :: critical_shields
      contains
         procedure :: validate => FluidErosion_validate
         procedure :: erosion => FluidErosion_erosion
   end type FluidErosion

   interface FluidErosion
      procedure FluidErosion_constructor
   end interface FluidErosion

   interface
      module function FluidErosion_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(FluidErosion) :: this
      end function FluidErosion_constructor

      module subroutine FluidErosion_validate(this)
         class(FluidErosion), intent(in) :: this
      end subroutine FluidErosion_validate

      module pure function FluidErosion_erosion(this, drag_model, RunParams, uvect) result(erosion)
         class(FluidErosion), intent(in) :: this
         class(DragModel), intent(in) :: drag_model
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: erosion
      end function FluidErosion_erosion
   end interface

!--- Granular erosion
!--- Erosion rate determined by the excess of the granular friction over a neutral friction
   type, extends(ErosionModel) :: GranularErosion
      real(kind=wp) :: granular_erosion_rate = 4.0_wp ! default
      real(kind=wp) :: pouliquen_min = 0.1_wp
      real(kind=wp) :: pouliquen_max = 0.4_wp
      real(kind=wp) :: pouliquen_static_slope
      real(kind=wp) :: pouliquen_beta = 0.136_wp
      real(kind=wp) :: t1
      contains
         procedure :: validate => GranularErosion_validate
         procedure :: erosion => GranularErosion_erosion
   end type GranularErosion

   interface GranularErosion
      procedure GranularErosion_constructor
   end interface GranularErosion

   interface
      module function GranularErosion_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(GranularErosion) :: this
      end function GranularErosion_constructor

      module subroutine GranularErosion_validate(this)
         class(GranularErosion), intent(in) :: this
      end subroutine GranularErosion_validate

      module pure function GranularErosion_erosion(this, drag_model, RunParams, uvect) result(erosion)
         class(GranularErosion), intent(in) :: this
         class(DragModel), intent(in) :: drag_model
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: erosion
      end function GranularErosion_erosion
   end interface

!--- Mixed Erosion
!--- Combines fluid and granular erosion via a switching function
   type, extends(ErosionModel) :: MixedErosion
      class(ErosionModel), allocatable :: fluid_erosion
      class(ErosionModel), allocatable :: granular_erosion
      real(kind=wp) :: voellmy_switch_rate = 3.0_wp
      real(kind=wp) :: voellmy_switch_value = 0.2_wp
      type(varString) :: switch_function
      procedure(fSwitch), pointer, nopass :: switcher => null()
      contains
         procedure :: validate => MixedErosion_validate
         procedure :: erosion => MixedErosion_erosion
   end type MixedErosion

   interface MixedErosion
      procedure MixedErosion_constructor
   end interface MixedErosion

   interface
      module function MixedErosion_constructor(params, fluid_erosion, granular_erosion) result(this)
         type(Dict), intent(in) :: params
         class(ErosionModel), intent(in) :: fluid_erosion
         class(ErosionModel), intent(in) :: granular_erosion
         type(MixedErosion) :: this
      end function MixedErosion_constructor

      module subroutine MixedErosion_validate(this)
         class(MixedErosion), intent(in) :: this
      end subroutine MixedErosion_validate

      module pure function MixedErosion_erosion(this, drag_model, RunParams, uvect) result(erosion)
         class(MixedErosion), intent(in) :: this
         class(DragModel), intent(in) :: drag_model
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: erosion
      end function MixedErosion_erosion
   end interface

!--- NoErosion
!--- Set erosion rate to zero.
   type, extends(ErosionModel) :: NoErosion
      contains
         procedure :: validate => NoErosion_validate
         procedure :: erosion => NoErosion_erosion
   end type NoErosion

   interface NoErosion
      procedure NoErosion_constructor
   end interface NoErosion

   interface
      module function NoErosion_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(NoErosion) :: this
      end function NoErosion_constructor

      module subroutine NoErosion_validate(this)
         class(NoErosion), intent(in) :: this
      end subroutine NoErosion_validate

      module pure function NoErosion_erosion(this, drag_model, RunParams, uvect) result(erosion)
         class(NoErosion), intent(in) :: this
         class(DragModel), intent(in) :: drag_model
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: erosion
      end function NoErosion_erosion
   end interface

contains

   subroutine create_erosion_model(choice, ParamsDict, erosion_model)
      character(len=*), intent(in) :: choice
      type(Dict), intent(in) :: ParamsDict
      class(ErosionModel), allocatable, intent(out) :: erosion_model

      if (len_trim(choice) == 0) then
         call WarningMessage("In the 'Parameters' block 'erosion' is not given." // new_line('a') // &
                                      " Using default 'erosion = off'")
         allocate(erosion_model, source=NoErosion(ParamsDict))
         return
      end if

      select case (trim(choice))
         case ('mixed', 'on')
            allocate(erosion_model, source=MixedErosion(ParamsDict, FluidErosion(ParamsDict), GranularErosion(ParamsDict)))

         case ('simple')
            allocate(erosion_model, source=SimpleErosion(ParamsDict))

         case ('fluid')
            allocate(erosion_model, source=FluidErosion(ParamsDict))

         case ('granular')
            allocate(erosion_model, source=GranularErosion(ParamsDict))

         case ('off')
            allocate(erosion_model, source=NoErosion(ParamsDict))

         case default
            call WarningMessage("In the 'Parameters' block 'erosion = " // choice // "' is not recognized." // new_line('a') // &
                                " Using default 'erosion = off'")
            allocate(erosion_model, source=NoErosion(ParamsDict))
      end select
   end subroutine create_erosion_model

end module erosion_types