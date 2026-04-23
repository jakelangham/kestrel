! This file is part of the Kestrel software for simulations
! of sediment-laden Earth surface flows.
!
! Version v1.1.2
!
! Copyright 2023 Mark J. Woodhouse, Jake Langham, (University of Bristol).
!
! This program is free software: you can redistribute it and/or modify it 
! under the terms of the GNU General Public License as published by the Free 
! Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT 
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
! more details.
!
! You should have received a copy of the GNU General Public License along with 
! this program. If not, see <https://www.gnu.org/licenses/>. 


! This module assigns values to variables in RunSet associated with the
! simulation parameters, as read from the "Parameters" block of an input file.
! Some input validation is performed, producing Warning or FatalError messages.
!
! The parameters includes choice of model closures by pointing the generic
! closure functions to specific implementations.  See Closures.f90
module parameters_module

   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage, InputLabelUnrecognized, WarningMessage
   use varstring_module, only: varString
   use dict_module
   use runsettings_module, only: RunSet
   use closures_module
   use drag_interface, only: DragModel, drag_model
   use drag_types

   implicit none

   private
   public :: Params_Set

   character(len=4), parameter :: fswitch_d = 'tanh'
   procedure(fswitch), pointer :: fswitch_dfunc => tanhSwitch
   character(len=4), parameter :: morpho_damp_d = 'tanh'
   procedure(MorphoDamping), pointer :: morpho_damp_dfunc => tanhMorphoDamping
   logical, parameter :: geometric_factors_d = .true.
   real(kind=wp), parameter :: g_d = 9.81_wp
   real(kind=wp), parameter :: chezyco_d = 0.01_wp
   real(kind=wp), parameter :: manningco_d = 0.03_wp
   real(kind=wp), parameter :: coulombco_d = 0.1_wp
   real(kind=wp), parameter :: powerlawco_d = 0.1_wp
   real(kind=wp), parameter :: powerlawpower_d = 1.0_wp
   real(kind=wp), parameter :: pouliquenMinSlope_d = 0.1_wp
   real(kind=wp), parameter :: pouliquenMaxSlope_d = 0.4_wp
   real(kind=wp), parameter :: pouliquenIntermediateSlope_d = 0.2_wp
   real(kind=wp), parameter :: pouliquen_beta_d=0.136_wp
   real(kind=wp), parameter :: edwards2019_betastar_d=0.136_wp
   real(kind=wp), parameter :: edwards2019_kappa_d=1.0_wp
   real(kind=wp), parameter :: edwards2019_Gamma_d=0.0_wp
   real(kind=wp), parameter :: EroRate_d = 0.001_wp
   real(kind=wp), parameter :: EroRateGranular_d = 4.0_wp
   real(kind=wp), parameter :: EroDepth_d = 1.0_wp

   character(len=16), parameter :: DepositionClosure_d = 'Spearman Manning'
   procedure(DepositionClosure), pointer :: DepositionClosure_dfunc => SpearmanManningHinderedSettling
   character(len=5), parameter :: ErosionClosure_d = 'Mixed'
   procedure(ErosionClosure), pointer :: ErosionClosure_dfunc => MixedErosion
   character(len=6), parameter :: EroTransition_d = 'smooth'
   procedure(ErosionTransition), pointer :: EroTransition_dfunc => SmoothErosionTransition
   real(kind=wp), parameter :: EroCriticalHeight_d = 0.01_wp
   real(kind=wp), parameter :: VoellmySwitchRate_d = 3.0_wp
   real(kind=wp), parameter :: VoellmySwitchValue_d = 0.2_wp
   real(kind=wp), parameter :: maxPack_d = 0.65_wp
   real(kind=wp), parameter :: rhow_d = 1000.0_wp ! Density of water
   real(kind=wp), parameter :: rhos_d = 2000.0_wp ! Density of solids
   real(kind=wp), parameter :: BedPorosity_d = 0.35_wp ! Porosity of bed
   real(kind=wp), parameter :: SolidDiameter_d = 1e-3_wp

   real(kind=wp), parameter :: visc_w = 1.2e-6_wp ! Kinematic viscosity of water
   real(kind=wp), parameter :: EddyViscosity_d = 0.0_wp

contains


   subroutine Params_Set(ParamsDict, RunParams)
      ! Set parameters from input file.
      ! Inputs: ParamsDict [type: Dict] - Dictionary (key-value lists) of parameters
      ! Output: Params [type ParamSet] - The parameters information.

      implicit none

      type(Dict), intent(inout) :: ParamsDict
      type(RunSet), intent(inout) :: RunParams

      type(varString) :: DragChoice, DepositionChoice, ErosionChoice
      type(varString) :: geometric_factors, switcher, eroTransition
      type(varString) :: morpho_damp

      integer :: J

      ! Set defaults if not set
      call ParamsDict%append('rhow', rhow_d)
      call ParamsDict%append('rhos', rhos_d)
      call ParamsDict%append('bed porosity', BedPorosity_d)
      call ParamsDict%append('maxPack', maxPack_d)
      call ParamsDict%append('eddy viscosity', EddyViscosity_d)
      call ParamsDict%append('erosion', 'off')
      call ParamsDict%append('solid diameter', SolidDiameter_d)
      call ParamsDict%append('g', g_d)
      call ParamsDict%append('switch function', fswitch_d)

      ! set gravity
      call ParamsDict%get('g', RunParams%g)

      ! set densities
      call ParamsDict%get('rhow', RunParams%rhow)
      call ParamsDict%get('rhos', RunParams%rhos)

      ! set solid diameter
      call ParamsDict%get('solid diameter', RunParams%SolidDiameter)

      ! set bed porosity
      call ParamsDict%get('bed porosity', RunParams%BedPorosity)

      ! set maxPack
      call ParamsDict%get('maxPack', RunParams%maxPack)

      ! set eddy viscosity
      call ParamsDict%get('eddy viscosity', RunParams%EddyViscosity)

      ! set erosion on/off
      call ParamsDict%get('erosion', RunParams%ErosionChoice)
      if (RunParams%ErosionChoice /= 'off') then
         RunParams%MorphodynamicsOn = .TRUE.
      else
         RunParams%MorphodynamicsOn = .FALSE.
      end if

      ! set fswitch
      call ParamsDict%get('switch function', RunParams%fswitch)
         end select
      end do

      if (.not.set_fswitch) then
         RunParams%fswitch = varString(fswitch_d)
         fswitch => fswitch_dfunc
      end if

      if (.not.set_EroTransition) then
         RunParams%ErosionTransition = varString(EroTransition_d)
         ErosionTransition => EroTransition_dfunc
      end if

      if (.not.set_MorphoDamp) then
         RunParams%MorphoDamp = varString(morpho_damp_d)
         MorphoDamping => morpho_damp_dfunc
      end if

      if (.not.set_geometric_factors) then
         RunParams%geometric_factors = geometric_factors_d
         if (geometric_factors_d) then
            GeometricCorrectionFactor => IversonOuyangGeometricCorrectionFactor
            GeometricCorrectionFactor_gradin => IversonOuyangGeometricCorrectionFactor_gradin
         else
            GeometricCorrectionFactor => NoGeometricCorrectionFactor
            GeometricCorrectionFactor_gradin => NoGeometricCorrectionFactor_gradin
         end if
      end if
      if (.not.set_g) RunParams%g = g_d
      if (.not.set_drag) then
         RunParams%DragChoice = varString('Chezy')
         DragClosure => ChezyDrag
         call WarningMessage("In the 'Parameters' block 'Drag' is not given.  Using default 'Chezy' drag.")
      end if

      if ((RunParams%DragChoice%s=='Chezy').or.(RunParams%DragChoice%s=="Variable")) then
         if (.not.set_chezyco) then
            RunParams%ChezyCo = chezyco_d
            call Warning_DragDefaultValue(RunParams%DragChoice%s,"Chezy Co",RunParams%ChezyCo)
         end if
      end if

      if ((RunParams%DragChoice%s=="Coulomb").and.(.not.set_coulombco)) then
         RunParams%CoulombCo = coulombco_d
         call Warning_DragDefaultValue("Coulomb","Coulomb Co",RunParams%CoulombCo)
      end if

      ! Set drag model
      call ParamsDict%get('drag', DragChoice)
      DragChoice = DragChoice%to_lower()   
      call create_drag_model(DragChoice%s, ParamsDict, drag_model)
      call drag_model%validate()
      end if

      RunParams%nsettling = (4.7_wp+0.41_wp*RunParams%Rep**(0.75_wp))/(1.0_wp + 0.175_wp*RunParams%Rep**(0.75_wp)) ! Exponent in hindered settling form (Rowe 1987; A convenient empirical equation for estimation of the Richardson-Zaki exponent)

      RunParams%CriticalShields = 0.3_wp/(1.0_wp+1.2_wp*R) + 0.055_wp*(1.0_wp - exp(-0.02_wp*R))

      RunParams%diffusiveTimeScale = huge(RunParams%diffusiveTimeScale)
      if (RunParams%EddyViscosity > 0.0_wp) then
          RunParams%diffusiveTimeScale = min( &
              RunParams%deltaX * RunParams%deltaX / RunParams%EddyViscosity, &
              RunParams%deltaY * RunParams%deltaY / RunParams%EddyViscosity)
      end if

      return

   end subroutine Params_Set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine FatalError_Positive(InputFile,var)

      implicit none

      character(len=*), intent(in) :: InputFile
      character(len=*), intent(in) :: var

      call FatalErrorMessage("In the 'Parameters' block in the input file " &
         // trim(InputFile) // new_line('A') &
         // " the block variable '" // trim(var) // "' must be positive.")
      return
   end subroutine FatalError_Positive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine FatalError_NonNegative(InputFile,var)

      implicit none

      character(len=*), intent(in) :: InputFile
      character(len=*), intent(in) :: var

      call FatalErrorMessage("In the 'Parameters' block in the input file " &
         // trim(InputFile) // new_line('A') &
         // " the block variable '" // trim(var) // "' must be non-negative.")
      return
   end subroutine FatalError_NonNegative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine FatalError_MaxValue(InputFile,var,num)

      implicit none

      character(len=*), intent(in) :: InputFile
      character(len=*), intent(in) :: var
      integer, intent(in) :: num

      call FatalErrorMessage("In the 'Parameters' block in the input file " &
         // trim(InputFile) // new_line('A') &
         // " the block variable '" // trim(var) // "' must be less than ", num)
      return
   end subroutine FatalError_MaxValue

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine Warning_DragDefaultValue(drag,var,num)

      implicit none

      character(len=*), intent(in) :: drag
      character(len=*), intent(in) :: var
      real(kind=wp), intent(in) :: num

      call WarningMessage("In the 'Parameters' block " // trim(drag) // " drag is given without a corresponding " &
         // trim(var) // " value.  Using default " // trim(var) // " = ", num)

      return
   end subroutine Warning_DragDefaultValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine Warning_DefaultValue(var,num)

      implicit none

      character(len=*), intent(in) :: var
      real(kind=wp), intent(in) :: num

      call WarningMessage("In the 'Parameters' block " // trim(var) // " is not set. Using default value " // trim(var) // " = ", num)

      return
   end subroutine Warning_DefaultValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine Warning_ErosionDefaultValue(var,num)

      implicit none

      character(len=*), intent(in) :: var
      real(kind=wp), intent(in) :: num

      call WarningMessage("In the 'Parameters' block 'Erosion' is switched 'On' without a corresponding '" // trim(var) // "' value." // " Using default " // trim(var) // " = ", num)

      return
   end subroutine Warning_ErosionDefaultValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module parameters_module
