module drag_types

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use messages_module, only: FatalErrorMessage, WarningMessage
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use drag_interface, only: DragModel
   use closures_module, only: fswitch

   implicit none

   private
   public :: ChezyDrag, ManningDrag, PowerlawDrag
   public :: CoulombDrag, PouliquenDrag, Edwards2019Drag
   public :: VariableDrag, VoellmyDrag
   public :: create_drag_model

   ! -- Specific drag model interfaces --
   !
   !    N.B. Because of the way the time stepper works (treating the drag
   !    implicitly), these are required to define a 'friction' function F such that
   !    the basal stress is \tau_x = -F \rho u/|(u,v)|, \tau_y = -F \rho v/|(u,v)|
   !    The basal stress are given by the 'drag_xy' function.
   !
   !    The 'validate' function does validation of model parameters
   !
   !    Implementations of the constructor, validator and friction are in submodules.

   ! Chezy
   type, extends(DragModel) :: ChezyDrag
      real(kind=wp) :: chezy_co = 0.01_wp
   contains
      procedure :: validate => chezy_validate
      procedure :: friction => chezy_friction
   end type ChezyDrag

   interface ChezyDrag
      module procedure chezy_constructor
   end interface ChezyDrag
   interface

      module function chezy_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(ChezyDrag) :: this
      end function chezy_constructor

      module subroutine chezy_validate(this)
         class(ChezyDrag), intent(in) :: this
      end subroutine chezy_validate

      module pure function chezy_friction(this, RunParams, uvect) result(friction)
         class(ChezyDrag), intent(in) :: this
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: friction
      end function chezy_friction

   end interface

   ! Coulomb
   type, extends(DragModel) :: CoulombDrag
      real(kind=wp) :: coulomb_co = 0.1_wp
   contains
      procedure :: validate => coulomb_validate
      procedure :: friction => coulomb_friction
   end type CoulombDrag

   interface CoulombDrag
      module procedure coulomb_constructor
   end interface CoulombDrag

   interface
      module function coulomb_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(CoulombDrag) :: this
      end function coulomb_constructor

      module subroutine coulomb_validate(this)
         class(CoulombDrag), intent(in) :: this
      end subroutine coulomb_validate

      module pure function coulomb_friction(this, RunParams, uvect) result(friction)
         class(CoulombDrag), intent(in) :: this
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: friction
      end function coulomb_friction
   end interface

   ! Edwards 2019
   type, extends(DragModel) :: Edwards2019Drag
      real(kind=wp) :: pouliquen_min = 0.1_wp
      real(kind=wp) :: pouliquen_max = 0.4_wp
      real(kind=wp) :: pouliquen_intermediate = 0.2_wp
      real(kind=wp) :: pouliquen_beta = 0.136_wp
      real(kind=wp) :: edwards2019_betastar = 0.136_wp
      real(kind=wp) :: edwards2019_kappa = 1.0_wp
      real(kind=wp) :: edwards2019_gamma = 0.0_wp
      real(kind=wp) :: solid_diameter
   contains
      procedure :: validate => edwards2019_validate
      procedure :: friction => edwards2019_friction
      procedure :: friction_coefficient => edwards2019_friction_coefficient
   end type Edwards2019Drag

   interface Edwards2019Drag
      module procedure edwards2019_constructor
   end interface Edwards2019Drag

   interface
      module function edwards2019_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(Edwards2019Drag) :: this
      end function edwards2019_constructor

      module subroutine edwards2019_validate(this)
         class(Edwards2019Drag), intent(in) :: this
      end subroutine edwards2019_validate

      module pure function edwards2019_friction(this, RunParams, uvect) result(friction)
         class(Edwards2019Drag), intent(in) :: this
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: friction
      end function edwards2019_friction

      module pure function edwards2019_friction_coefficient(this, heightThreshold, gcostheta, Hn, modu) result(mu)
         class(Edwards2019Drag), intent(in) :: this
         real(kind=wp), intent(in) :: heightThreshold
         real(kind=wp), intent(in) :: gcostheta
         real(kind=wp), intent(in) :: Hn
         real(kind=wp), intent(in) :: modu
         real(kind=wp) :: mu
      end function edwards2019_friction_coefficient
   end interface

   ! Manning
   type, extends(DragModel) :: ManningDrag
      real(kind=wp) :: manning_co = 0.03_wp
   contains
      procedure :: validate => manning_validate
      procedure :: friction => manning_friction
   end type ManningDrag

   interface ManningDrag
      module procedure manning_constructor
   end interface ManningDrag

   interface
      module function manning_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(ManningDrag) :: this
      end function manning_constructor

      module subroutine manning_validate(this)
         class(ManningDrag), intent(in) :: this
      end subroutine manning_validate

      module pure function manning_friction(this, RunParams, uvect) result(friction)
         class(ManningDrag), intent(in) :: this
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: friction
      end function manning_friction
   end interface

   ! Pouliquen
   type, extends(DragModel) :: PouliquenDrag
      real(kind=wp) :: pouliquen_min = 0.1_wp
      real(kind=wp) :: pouliquen_max = 0.4_wp
      real(kind=wp) :: pouliquen_beta = 0.136_wp
      real(kind=wp) :: solid_diameter
   contains
      procedure :: validate => pouliquen_validate
      procedure :: friction => pouliquen_friction
      procedure :: friction_coefficient => pouliquen_friction_coefficient
   end type PouliquenDrag

   interface PouliquenDrag
      module procedure pouliquen_constructor
   end interface PouliquenDrag

   interface
      module function pouliquen_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(PouliquenDrag) :: this
      end function pouliquen_constructor

      module subroutine pouliquen_validate(this)
         class(PouliquenDrag), intent(in) :: this
      end subroutine pouliquen_validate

      module pure function pouliquen_friction(this, RunParams, uvect) result(friction)
         class(PouliquenDrag), intent(in) :: this
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: friction
      end function pouliquen_friction

      module pure function pouliquen_friction_coefficient(this, heightThreshold, gcostheta, Hn, modu) result(mu)
         class(PouliquenDrag), intent(in) :: this
         real(kind=wp), intent(in) :: heightThreshold
         real(kind=wp), intent(in) :: gcostheta
         real(kind=wp), intent(in) :: Hn
         real(kind=wp), intent(in) :: modu
         real(kind=wp) :: mu
      end function pouliquen_friction_coefficient
   end interface

   ! Powerlaw
   type, extends(DragModel) :: PowerlawDrag
      real(kind=wp) :: power_law_co = 0.1_wp
      real(kind=wp) :: power_law_power = 1.0_wp
      logical :: set_power_law_co = .FALSE.
      logical :: set_power_law_power = .FALSE.
   contains
      procedure :: validate => powerlaw_validate
      procedure :: friction => powerlaw_friction
   end type PowerlawDrag

   interface PowerlawDrag
      module procedure powerlaw_constructor
   end interface PowerlawDrag

   interface
      module function powerlaw_constructor(params) result(this)
         type(Dict), intent(in) :: params
         type(PowerlawDrag) :: this
      end function powerlaw_constructor

      module subroutine powerlaw_validate(this)
         class(PowerlawDrag), intent(in) :: this
      end subroutine powerlaw_validate

      module pure function powerlaw_friction(this, RunParams, uvect) result(friction)
         class(PowerlawDrag), intent(in) :: this
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: friction
      end function powerlaw_friction
   end interface

   ! Variable
   type, extends(DragModel) :: VariableDrag

      class(DragModel), allocatable :: fluid
      class(DragModel), allocatable :: granular
      real(kind=wp) :: voellmy_switch_rate = 3.0_wp
      real(kind=wp) :: voellmy_switch_value = 0.2_wp
      type(varString) :: switch_function
      procedure(fSwitch), pointer, nopass :: switcher => null()
   contains
      procedure :: validate => variable_validate
      procedure :: friction => variable_friction
   end type VariableDrag

   interface VariableDrag
      module procedure variable_constructor
   end interface VariableDrag

   interface
      module function variable_constructor(params, fluid_drag, granular_drag) result(this)
         type(Dict), intent(in) :: params   
         class(DragModel), intent(in) :: fluid_drag
         class(DragModel), intent(in) :: granular_drag
         type(VariableDrag) :: this
      end function variable_constructor

      module subroutine variable_validate(this)
         class(VariableDrag), intent(in) :: this
      end subroutine variable_validate

      module pure function variable_friction(this, RunParams, uvect) result(friction)
         class(VariableDrag), intent(in) :: this
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: friction
      end function variable_friction
   end interface

   ! Voellmy
   type, extends(DragModel) :: VoellmyDrag

      type(ChezyDrag) :: chezy_part = ChezyDrag()
      type(CoulombDrag) :: coulomb_part = CoulombDrag()
   contains
      procedure :: validate => voellmy_validate
      procedure :: friction => voellmy_friction
   end type VoellmyDrag

   interface VoellmyDrag
      module procedure voellmy_constructor
   end interface VoellmyDrag

   interface
      module function voellmy_constructor(chezy_part, coulomb_part) result(this)
         type(ChezyDrag), intent(in) :: chezy_part
         type(CoulombDrag), intent(in) :: coulomb_part
         type(VoellmyDrag) :: this
      end function voellmy_constructor

      module subroutine voellmy_validate(this)
         class(VoellmyDrag), intent(in) :: this
      end subroutine voellmy_validate

      module pure function voellmy_friction(this, RunParams, uvect) result(friction)
         class(VoellmyDrag), intent(in) :: this
         type(RunSet), intent(in) :: RunParams
         real(kind=wp), dimension(:), intent(in) :: uvect
         real(kind=wp) :: friction
      end function voellmy_friction
   end interface

contains

   subroutine create_drag_model(choice, ParamsDict, drag_model)
      character(len=*), intent(in) :: choice
      type(Dict), intent(in) :: ParamsDict
      class(DragModel), allocatable, intent(out) :: drag_model

      if (len_trim(choice) == 0) then
         call WarningMessage("In the 'Parameters' block 'drag' is not given." // new_line('a') // &
                             "Using default 'Chezy' drag.")
         allocate(drag_model, source=ChezyDrag(ParamsDict))
         return
      end if

      select case (trim(choice))
         case ('chezy')
            allocate(drag_model, source=ChezyDrag(ParamsDict))
         case ('coulomb')
            allocate(drag_model, source=CoulombDrag(ParamsDict))
         case ('manning')
            allocate(drag_model, source=ManningDrag(ParamsDict))
         case ('pouliquen')
            allocate(drag_model, source=PouliquenDrag(ParamsDict))
         case ('power law')
            allocate(drag_model, source=PowerlawDrag(ParamsDict))
         case ('edwards2019')
            allocate(drag_model, source=Edwards2019Drag(ParamsDict))
         case ('variable')
            allocate(drag_model, source=VariableDrag(ParamsDict, &
                     ChezyDrag(ParamsDict), PouliquenDrag(ParamsDict)))
         case ('voellmy')
            allocate(drag_model, source=VoellmyDrag(&
                     ChezyDrag(ParamsDict), CoulombDrag(ParamsDict)))
         case default
            call WarningMessage("In the 'Parameters' block 'drag = " // trim(choice) // "' is not recognized." // new_line('a') // &
                                "Using default 'Chezy' drag.")

            allocate(drag_model, source=ChezyDrag(ParamsDict))
      end select
   end subroutine create_drag_model

end module drag_types