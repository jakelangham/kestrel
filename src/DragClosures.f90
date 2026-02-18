module drag_closures_module

   use set_precision_module, only: wp
   use messages_module, only: WarningMessage
   use varstring_module, only: varString

   implicit none

   type, abstract :: DragModel
      type(varString) :: name
   contains
      procedure(validate_proc), deferred :: validate
   end type DragModel

   abstract interface
      subroutine validate_proc(this)
         import :: DragModel
         class(DragModel), intent(inout) :: this
      end subroutine validate_proc
   end interface

   ! Chezy
   type, extends(DragModel) :: ChezyDrag
      real(kind=wp) :: chezy_co = 0.01_wp
      logical :: set_chezy_co = .FALSE.
   contains
      procedure :: validate => chezy_validate
   end type ChezyDrag

   interface ChezyDrag
      procedure chezy_constructor
   end interface ChezyDrag

   ! Manning
   type, extends(DragModel) :: ManningDrag
      real(kind=wp) :: manning_co = 0.03_wp
      logical :: set_manning_co = .FALSE.
   contains
      procedure :: validate => manning_validate
   end type ManningDrag

   interface ManningDrag
      procedure manning_constructor
   end interface ManningDrag

   ! Powerlaw
   type, extends(DragModel) :: PowerlawDrag
      real(kind=wp) :: power_law_co = 0.1_wp
      real(kind=wp) :: power_law_power = 1.0_wp
      logical :: set_power_law_co = .FALSE.
      logical :: set_power_law_power = .FALSE.
   contains
      procedure :: validate => powerlaw_validate
   end type PowerlawDrag

   interface PowerlawDrag
      procedure powerlaw_constructor
   end interface PowerlawDrag

   ! Coulomb
   type, extends(DragModel) :: CoulombDrag
      real(kind=wp) :: coulomb_co = 0.1_wp
      logical :: set_coulomb_co = .FALSE.
   contains
      procedure :: validate => coulomb_validate
   end type CoulombDrag

   interface CoulombDrag
      procedure coulomb_constructor
   end interface CoulombDrag

   ! Pouliquen
   type, extends(DragModel) :: PouliquenDrag
      real(kind=wp) :: pouliquen_min = 0.1_wp
      real(kind=wp) :: pouliquen_max = 0.4_wp
      real(kind=wp) :: pouliquen_beta = 0.136_wp
      logical :: set_pouliquen_min = .FALSE.
      logical :: set_pouliquen_max = .FALSE.
      logical :: set_pouliquen_beta = .FALSE.
   contains
      procedure :: validate => pouliquen_validate
   end type PouliquenDrag

   interface PouliquenDrag
      procedure pouliquen_constructor
   end interface PouliquenDrag

   ! Edwards 2019
   type, extends(DragModel) :: Edwards2019Drag
      real(kind=wp) :: pouliquen_min = 0.1_wp
      real(kind=wp) :: pouliquen_max = 0.4_wp
      real(kind=wp) :: pouliquen_intermediate = 0.2_wp
      real(kind=wp) :: pouliquen_beta = 0.136_wp
      real(kind=wp) :: edwards2019_betastar = 0.136_wp
      real(kind=wp) :: edwards2019_kappa = 1.0_wp
      real(kind=wp) :: edwards2019_gamma = 0.0_wp
      logical :: set_pouliquen_min = .FALSE.
      logical :: set_pouliquen_max = .FALSE.
      logical :: set_pouliquen_intermediate = .FALSE.
      logical :: set_pouliquen_beta = .FALSE.
      logical :: set_edwards2019_betastar = .FALSE.
      logical :: set_edwards2019_kappa = .FALSE.
      logical :: set_edwards2019_gamma = .FALSE.
   contains
      procedure :: validate => edwards2019_validate
   end type Edwards2019Drag

   interface Edwards2019Drag
      procedure edwards2019_constructor
   end interface Edwards2019Drag

   ! Variable
   type, extends(DragModel) :: VariableDrag
      real(kind=wp) :: chezy_co = 0.01_wp
      real(kind=wp) :: pouliquen_min = 0.1_wp
      real(kind=wp) :: pouliquen_max = 0.4_wp
      real(kind=wp) :: pouliquen_beta = 0.136_wp
      real(kind=wp) :: voellmy_switch_rate = 3.0_wp
      real(kind=wp) :: voellmy_switch_value = 0.2_wp
      logical :: set_chezy_co = .FALSE.
      logical :: set_pouliquen_min = .FALSE.
      logical :: set_pouliquen_max = .FALSE.
      logical :: set_pouliquen_beta = .FALSE.
      logical :: set_voellmy_switch_rate = .FALSE.
      logical :: set_voellmy_switch_value = .FALSE.
   contains
      procedure :: validate => variable_validate
   end type VariableDrag

   interface VariableDrag
      procedure variable_constructor
   end interface VariableDrag

   ! Voellmy
   type, extends(DragModel) :: VoellmyDrag
      real(kind=wp) :: chezy_co = 0.01_wp
      real(kind=wp) :: coulomb_co = 0.1_wp
      logical :: set_chezy_co = .FALSE.
      logical :: set_coulomb_co = .FALSE.
   contains
      procedure :: validate => voellmy_validate
   end type VoellmyDrag

   interface VoellmyDrag
      procedure voellmy_constructor
   end interface VoellmyDrag


contains

   ! Chezy
   pure function chezy_constructor(chezy_co) result(this)
      real(kind=wp), intent(in), optional :: chezy_co
      type(ChezyDrag) :: this
      this%name = varString("Chezy")
      if (present(chezy_co)) then
         this%chezy_co = chezy_co
         this%set_chezy_co = .TRUE.
      end if
   end function chezy_constructor

   subroutine chezy_validate(this)
      class(ChezyDrag), intent(inout) :: this
      if (.not. this%set_chezy_co) then
         call Warning_DragDefaultValue(this%name%s, 'Chezy Co', this%chezy_co)
      end if
   end subroutine chezy_validate

   ! Manning
   pure function manning_constructor(manning_co) result(this)
      real(kind=wp), intent(in), optional :: manning_co
      type(ManningDrag) :: this
      this%name = varString("Manning")
      if (present(manning_co)) then
         this%manning_co = manning_co
         this%set_manning_co = .TRUE.
      end if
   end function manning_constructor

   subroutine manning_validate(this)
      class(ManningDrag), intent(inout) :: this
      if (.not. this%set_manning_co) then
         call Warning_DragDefaultValue(this%name%s, 'Manning Co', this%manning_co)
      end if
   end subroutine manning_validate

   ! Powerlaw
   pure function powerlaw_constructor(powerlaw_co, powerlaw_power) result(this)
      ! NOTE - expect named arguments are used in call to constructor
      real(kind=wp), intent(in), optional :: powerlaw_co
      real(kind=wp), intent(in), optional :: powerlaw_power
      type(PowerlawDrag) :: this
      this%name = varString("Powerlaw")
      if (present(powerlaw_co)) then
         this%power_law_co = powerlaw_co
         this%set_power_law_co = .TRUE.
      end if
      if (present(powerlaw_power)) then
         this%power_law_power = powerlaw_power
         this%set_power_law_power = .TRUE.
      end if
   end function powerlaw_constructor

   subroutine powerlaw_validate(this)
      class(PowerlawDrag), intent(inout) :: this
      if (.not. this%set_power_law_co) then
         call Warning_DragDefaultValue(this%name%s, 'Power law co', this%power_law_co)
      end if
      if (.not. this%set_power_law_power) then
         call Warning_DragDefaultValue(this%name%s, 'Power law power', this%power_law_power)
      end if
   end subroutine powerlaw_validate

   ! Coulomb
   pure function coulomb_constructor(coulomb_co) result(this)
      real(kind=wp), intent(in), optional :: coulomb_co
      type(CoulombDrag) :: this
      this%name = varString("Coulomb")
      if (present(coulomb_co)) then
         this%coulomb_co = coulomb_co
         this%set_coulomb_co = .TRUE.
      end if
   end function coulomb_constructor

   subroutine coulomb_validate(this)
      class(CoulombDrag), intent(inout) :: this
      if (.not. this%set_coulomb_co) then
         call Warning_DragDefaultValue(this%name%s, 'Coulomb Co', this%coulomb_co)
      end if
   end subroutine coulomb_validate

   ! Pouliquen
   pure function pouliquen_constructor(pouliquen_min, pouliquen_max, pouliquen_beta) result(this)
      real(kind=wp), intent(in), optional :: pouliquen_min
      real(kind=wp), intent(in), optional :: pouliquen_max
      real(kind=wp), intent(in), optional :: pouliquen_beta
      type(PouliquenDrag) :: this
      this%name = varString("Pouliquen")
      if (present(pouliquen_min)) then
         this%pouliquen_min = pouliquen_min
         this%set_pouliquen_min = .TRUE.
      end if
      if (present(pouliquen_max)) then
         this%pouliquen_max = pouliquen_max
         this%set_pouliquen_max = .TRUE.
      end if
      if (present(pouliquen_beta)) then
         this%pouliquen_beta = pouliquen_beta
         this%set_pouliquen_beta = .TRUE.
      end if
   end function pouliquen_constructor

   subroutine pouliquen_validate(this)
      class(PouliquenDrag), intent(inout) :: this
      if (.not. this%set_pouliquen_min) then
         call Warning_DragDefaultValue(this%name%s, 'Pouliquen min', this%pouliquen_min)
      end if
      if (.not. this%set_pouliquen_max) then
         call Warning_DragDefaultValue(this%name%s, 'Pouliquen max', this%pouliquen_max)
      end if
      if (.not. this%set_pouliquen_beta) then
         call Warning_DragDefaultValue(this%name%s, 'Pouliquen beta', this%pouliquen_beta)
      end if
   end subroutine pouliquen_validate

   ! Edwards 2019
   pure function edwards2019_constructor( &
      pouliquen_min, &
      pouliquen_max, &
      pouliquen_intermediate, &
      pouliquen_beta, &
      edwards2019_betastar, &
      edwards2019_kappa, &
      edwards2019_gamma &
   ) result(this)
      real(kind=wp), intent(in), optional :: pouliquen_min
      real(kind=wp), intent(in), optional :: pouliquen_max
      real(kind=wp), intent(in), optional :: pouliquen_intermediate
      real(kind=wp), intent(in), optional :: pouliquen_beta
      real(kind=wp), intent(in), optional :: edwards2019_betastar
      real(kind=wp), intent(in), optional :: edwards2019_kappa
      real(kind=wp), intent(in), optional :: edwards2019_gamma
      type(Edwards2019Drag) :: this
      this%name = varString("Edwards2019")
      if (present(pouliquen_min)) then
         this%pouliquen_min = pouliquen_min
         this%set_pouliquen_min = .TRUE.
      end if
      if (present(pouliquen_max)) then
         this%pouliquen_max = pouliquen_max
         this%set_pouliquen_max = .TRUE.
      end if
      if (present(pouliquen_intermediate)) then
         this%pouliquen_intermediate = pouliquen_intermediate
         this%set_pouliquen_intermediate = .TRUE.
      end if
      if (present(pouliquen_beta)) then
         this%pouliquen_beta = pouliquen_beta
         this%set_pouliquen_beta = .TRUE.
      end if
      if (present(edwards2019_betastar)) then
         this%edwards2019_betastar = edwards2019_betastar
         this%set_edwards2019_betastar = .TRUE.
      end if
      if (present(edwards2019_kappa)) then
         this%edwards2019_kappa = edwards2019_kappa
         this%set_edwards2019_kappa = .TRUE.
      end if
      if (present(edwards2019_gamma)) then
         this%edwards2019_gamma = edwards2019_gamma
         this%set_edwards2019_gamma = .TRUE.
      end if
   end function edwards2019_constructor

   subroutine edwards2019_validate(this)
      class(Edwards2019Drag), intent(inout) :: this
      if (.not. this%set_pouliquen_min) then
         call Warning_DragDefaultValue(this%name%s, 'Pouliquen min', this%pouliquen_min)
      end if
      if (.not. this%set_pouliquen_max) then
         call Warning_DragDefaultValue(this%name%s, 'Pouliquen max', this%pouliquen_max)
      end if
      if (.not. this%set_pouliquen_intermediate) then
         call Warning_DragDefaultValue(this%name%s, 'pouliquen intermediate', this%pouliquen_intermediate)
      end if
      if (.not. this%set_pouliquen_beta) then
         call Warning_DragDefaultValue(this%name%s, 'pouliquen beta', this%pouliquen_beta)
      end if
      if (.not. this%set_edwards2019_betastar) then
         call Warning_DragDefaultValue(this%name%s, 'edwards2019 betastar', this%edwards2019_betastar)
      end if
      if (.not. this%set_edwards2019_kappa) then
         call Warning_DragDefaultValue(this%name%s, 'edwards2019 kappa', this%edwards2019_kappa)
      end if
      if (.not. this%set_edwards2019_gamma) then
         call Warning_DragDefaultValue(this%name%s, 'edwards2019 gamma', this%edwards2019_gamma)
      end if
   end subroutine edwards2019_validate

   ! Variable
   pure function variable_constructor( &
      chezy_co, &
      pouliquen_min, &
      pouliquen_max, &
      pouliquen_beta, &
      voellmy_switch_rate, &
      voellmy_switch_value &
   ) result(this)
      real(kind=wp), intent(in), optional :: chezy_co
      real(kind=wp), intent(in), optional :: pouliquen_min
      real(kind=wp), intent(in), optional :: pouliquen_max
      real(kind=wp), intent(in), optional :: pouliquen_beta
      real(kind=wp), intent(in), optional :: voellmy_switch_rate
      real(kind=wp), intent(in), optional :: voellmy_switch_value
      type(VariableDrag) :: this
      this%name = varString("Variable")
      if (present(chezy_co)) then
         this%chezy_co = chezy_co
         this%set_chezy_co = .TRUE.
      end if
      if (present(pouliquen_min)) then
         this%pouliquen_min = pouliquen_min
         this%set_pouliquen_min = .TRUE.
      end if
      if (present(pouliquen_max)) then
         this%pouliquen_max = pouliquen_max
         this%set_pouliquen_max = .TRUE.
      end if
      if (present(pouliquen_beta)) then
         this%pouliquen_beta = pouliquen_beta
         this%set_pouliquen_beta = .TRUE.
      end if
      if (present(voellmy_switch_rate)) then
         this%voellmy_switch_rate = voellmy_switch_rate
         this%set_voellmy_switch_rate = .TRUE.
      end if
      if (present(voellmy_switch_value)) then
         this%voellmy_switch_value = voellmy_switch_value
         this%set_voellmy_switch_value = .TRUE.
      end if
   end function variable_constructor

   subroutine variable_validate(this)
      class(VariableDrag), intent(inout) :: this
      if (.not. this%set_chezy_co) then
         call Warning_DragDefaultValue(this%name%s, 'chezy co', this%chezy_co)
      end if
      if (.not. this%set_pouliquen_min) then
         call Warning_DragDefaultValue(this%name%s, 'pouliquen min', this%pouliquen_min)
      end if
      if (.not. this%set_pouliquen_max) then
         call Warning_DragDefaultValue(this%name%s, 'pouliquen max', this%pouliquen_max)
      end if
      if (.not. this%set_pouliquen_beta) then
         call Warning_DragDefaultValue(this%name%s, 'pouliquen beta', this%pouliquen_beta)
      end if
      if (.not. this%set_voellmy_switch_rate) then
         call Warning_DragDefaultValue(this%name%s, 'voellmy switch rate', this%voellmy_switch_rate)
      end if
      if (.not. this%set_voellmy_switch_value) then
         call Warning_DragDefaultValue(this%name%s, 'voellmy switch value', this%voellmy_switch_value)
      end if
   end subroutine variable_validate

   ! Voellmy
   pure function voellmy_constructor(chezy_co, coulomb_co) result(this)
      real(kind=wp), intent(in), optional :: chezy_co
      real(kind=wp), intent(in), optional :: coulomb_co
      type(VoellmyDrag) :: this
      this%name = varString("Voellmy")
      if (present(chezy_co)) then
         this%chezy_co = chezy_co
         this%set_chezy_co = .TRUE.
      end if
      if (present(coulomb_co)) then
         this%coulomb_co = coulomb_co
         this%set_coulomb_co = .TRUE.
      end if
   end function voellmy_constructor

   subroutine voellmy_validate(this)
      class(VoellmyDrag), intent(inout) :: this
      if (.not. this%set_chezy_co) then
         call Warning_DragDefaultValue(this%name%s, 'chezy co', this%chezy_co)
      end if
      if (.not. this%set_coulomb_co) then
         call Warning_DragDefaultValue(this%name%s, 'pouliquen min', this%coulomb_co)
      end if
   end subroutine voellmy_validate


   subroutine Warning_DragDefaultValue(drag,var,num)

      implicit none

      character(len=*), intent(in) :: drag
      character(len=*), intent(in) :: var
      real(kind=wp), intent(in) :: num

      call WarningMessage("In the 'Parameters' block " // trim(drag) // " drag is given without a corresponding " &
         // trim(var) // " value.  Using default " // trim(var) // " = ", num)

      return
   end subroutine Warning_DragDefaultValue

end module drag_closures_module