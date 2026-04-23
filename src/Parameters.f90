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
   use erosion_interface, only: ErosionModel, erosion_model
   use erosion_types
   use deposition_interface, only: DepositionModel, deposition_model
   use deposition_types

   implicit none

   private
   public :: Params_Set

   character(len=4), parameter :: fswitch_d = 'tanh'
   character(len=4), parameter :: morpho_damp_d = 'tanh'
   procedure(MorphoDamping), pointer :: morpho_damp_dfunc => tanhMorphoDamping

   logical, parameter :: geometric_factors_d = .false.
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
      type(varString) :: geometric_factors
      type(varString) :: morpho_damp

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

      ! set morphodynamic parameters, if needed
      if (RunParams%MorphodynamicsOn) then
         ! set Erosion critical height
         call ParamsDict%append('erosion critical height', EroCriticalHeight_d)
         call ParamsDict%get('erosion critical height', RunParams%EroCriticalHeight)

         ! set erosion depth
         call ParamsDict%append('erosion depth', EroDepth_d)
         call ParamsDict%get('erosion depth', RunParams%EroDepth)

         ! set morphodynamic damping
         if (ParamsDict%has_key('morphodynamic damping')) then
            
            call ParamsDict%get('morphodynamic damping', morpho_damp)
            morpho_damp = morpho_damp%to_lower()
            select case (morpho_damp%s)
               case ('none', 'off')
                  RunParams%MorphoDamp = varString('None')
                  MorphoDamping => NoMorphoDamping
               case ('tanh')
                  RunParams%MorphoDamp = varString('tanh')
                  MorphoDamping => tanhMorphoDamping
               case ('rat3')
                  RunParams%MorphoDamp = varString('rat3')
                  MorphoDamping => rat3MorphoDamping
               case default
                  RunParams%MorphoDamp = varString(morpho_damp_d)
                  MorphoDamping => morpho_damp_dfunc
            end select
         else
            RunParams%MorphoDamp = varString(morpho_damp_d)
            MorphoDamping => morpho_damp_dfunc
         end if
      end if

      ! set geometric factors
      if (ParamsDict%has_key('geometric factors')) then
         call ParamsDict%get('geometric factors', geometric_factors)
         geometric_factors = geometric_factors%to_lower()
         !0=Off, 1=On (perhaps will change in future for different formulations
         select case (geometric_factors%s)
            case ('off')
               RunParams%geometric_factors = .false.
               GeometricCorrectionFactor => NoGeometricCorrectionFactor
               GeometricCorrectionFactor_gradin => NoGeometricCorrectionFactor_gradin
            case ('on')
               RunParams%geometric_factors = .true.
               GeometricCorrectionFactor => IversonOuyangGeometricCorrectionFactor
               GeometricCorrectionFactor_gradin => IversonOuyangGeometricCorrectionFactor_gradin
            case default
               call WarningMessage("In the 'Parameters' block 'geometric factors = " // geometric_factors%s // "' is not recognized. " // new_line('a') // &
                                   " Using default 'geometric factors = off'")
               RunParams%geometric_factors = geometric_factors_d
               GeometricCorrectionFactor => NoGeometricCorrectionFactor
               GeometricCorrectionFactor_gradin => NoGeometricCorrectionFactor_gradin
         end select
      else
         call WarningMessage("In the 'Parameters' block 'geometric factors' is not given. " // new_line('a') // &
                             " Using default 'geometric factors = off'")
         RunParams%geometric_factors = geometric_factors_d
         GeometricCorrectionFactor => NoGeometricCorrectionFactor
         GeometricCorrectionFactor_gradin => NoGeometricCorrectionFactor_gradin
      end if

      ! Set drag model
      call ParamsDict%get('drag', DragChoice)
      DragChoice = DragChoice%to_lower()   
      call create_drag_model(DragChoice%s, ParamsDict, drag_model)
      call drag_model%validate()
      

      if (RunParams%MorphodynamicsOn) then
         ! set deposition closure
         call ParamsDict%get('deposition', DepositionChoice)
         DepositionChoice = DepositionChoice%to_lower()
         call create_deposition_model(DepositionChoice%s, ParamsDict, RunParams, deposition_model)
         call deposition_model%validate()

         ! set erosion closure
         call ParamsDict%get('erosion', ErosionChoice)
         ErosionChoice = ErosionChoice%to_lower()
         call create_erosion_model(ErosionChoice%s, ParamsDict, erosion_model)
         call erosion_model%validate()
      end if

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
