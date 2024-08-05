! This file is part of the Kestrel software for simulations
! of sediment-laden Earth surface flows.
!
! Version 1.0
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
   use runsettings_module, only: RunSet
   use closures_module

   implicit none

   private
   public :: Params_Set

   character(len=4), parameter :: fswitch_d = 'tanh'
   procedure(fswitch), pointer :: fswitch_dfunc => tanhSwitch
   character(len=4), parameter :: morpho_damp_d = 'tanh'
   procedure(MorphoDamping), pointer :: morpho_damp_dfunc => tanhMorphoDamping
   logical, parameter :: geometric_factors_d = .true.
   logical, parameter :: curvature_d = .true.
   real(kind=wp), parameter :: g_d = 9.81_wp
   real(kind=wp), parameter :: chezyco_d = 0.01_wp
   real(kind=wp), parameter :: manningco_d = 0.03_wp
   real(kind=wp), parameter :: coulombco_d = 0.1_wp
   real(kind=wp), parameter :: pouliquenMinSlope_d = 0.1_wp
   real(kind=wp), parameter :: pouliquenMaxSlope_d = 0.4_wp
   real(kind=wp), parameter :: pouliquen_beta_d=0.136_wp
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


   subroutine Params_Set(ParamLabels,ParamValues,RunParams)
      ! Set parameters from input file.
      ! Inputs: ParamLabels [type: varString] - Labels for Parameters;
      !         ParamValues [type: varString] - The values associated with the labels.
      ! Output: Params [type ParamSet] - The parameters information.

      implicit none

      type(varString), dimension(:), intent(in) :: ParamLabels
      type(varString), dimension(:), intent(in) :: ParamValues
      type(RunSet), intent(inout) :: RunParams

      type(varString) :: label
      type(varString) :: DragChoice, DepositionChoice, ErosionChoice
      type(varString) :: geometric_factors, curvature, switcher, eroTransition
      type(varString) :: morpho_damp

      integer :: J, N

      logical :: set_fswitch
      logical :: set_g
      logical :: set_drag
      logical :: set_chezyco
      logical :: set_manningco
      logical :: set_coulombco
      logical :: set_pouliquenMinSlope
      logical :: set_pouliquenMaxSlope
      logical :: set_pouliquen_beta
      logical :: set_pouliquen_L
      logical :: set_deposition
      logical :: set_erosion
      logical :: set_EroRate
      logical :: set_EroRateGranular
      logical :: set_EroDepth
      logical :: set_EroTransition
      logical :: set_EroCriticalHeight
      logical :: set_VoellmySwitchRate
      logical :: set_VoellmySwitchValue
      logical :: set_maxPack
      logical :: set_rhow
      logical :: set_rhos
      logical :: set_BedPorosity
      logical :: set_SolidDiameter
      logical :: set_EddyViscosity
      logical :: set_geometric_factors
      logical :: set_curvature
      logical :: set_SettlingSpeed
      logical :: set_MorphoDamp

      real(kind=wp) :: R

      N = size(ParamValues)

      set_fswitch=.FALSE.
      set_g=.FALSE.
      set_drag=.FALSE.
      set_chezyco=.FALSE.
      set_manningco=.FALSE.
      set_coulombco=.FALSE.
      set_pouliquenMinSlope=.FALSE.
      set_pouliquenMaxSlope=.FALSE.
      set_pouliquen_beta=.FALSE.
      set_pouliquen_L=.FALSE.
      set_deposition=.FALSE.
      set_erosion=.FALSE.
      set_EroRate=.FALSE.
      set_EroRateGranular=.FALSE.
      set_EroDepth=.FALSE.
      set_EroTransition=.FALSE.
      set_EroCriticalHeight=.FALSE.
      set_VoellmySwitchRate=.FALSE.
      set_VoellmySwitchValue=.FALSE.
      set_maxPack=.FALSE.
      set_rhow=.FALSE.
      set_rhos=.FALSE.
      set_BedPorosity=.FALSE.
      set_SolidDiameter=.FALSE.
      set_EddyViscosity=.FALSE.
      set_geometric_factors = .false.
      set_curvature = .false.
      set_SettlingSpeed = .false.
      set_MorphoDamp = .false.

      do J=1,N
         label = ParamLabels(J)%to_lower()
         select case (label%s)

            case ('switch function')
               set_fswitch=.TRUE.
               switcher = ParamValues(J)%to_lower()
               select case (switcher%s)
                  case ('tanh')
                     RunParams%fswitch = varString('tanh')
                     fswitch => tanhSwitch
                  case ('rat3')
                     RunParams%fswitch = varString('rat3')
                     fswitch => rat3Switch
                  case ('cos')
                     RunParams%fswitch = varString('cos')
                     fswitch => cosSwitch
                  case ('linear')
                     RunParams%fswitch = varString('linear')
                     fswitch => linearSwitch
                  case ('equal', '0.5')
                     RunParams%fswitch = varString('equal')
                     fswitch => equalSwitch
                  case ('off', '0', 'zero')
                     RunParams%fswitch = varString('zero')
                     fswitch => zeroSwitch
                  case ('1', 'one')
                     RunParams%fswitch = varString('one')
                     fswitch => oneSwitch
                  case ('step')
                     RunParams%fswitch = varString('step')
                     fswitch => stepSwitch
                  case default
                     RunParams%fswitch = varString(fswitch_d)
                     fswitch => fswitch_dfunc
               end select
         
            case ('deposition')
               set_deposition=.true.
               DepositionChoice = ParamValues(J)%to_lower()
               select case (DepositionChoice%s)
                  case ('none')
                     RunParams%DepositionChoice = varString('None')
                     DepositionClosure => NoDeposition
                  case ('simple')
                     RunParams%DepositionChoice = varString('Simple')
                     DepositionClosure => SimpleHinderedSettling
                  case ('spearman manning')
                     RunParams%DepositionChoice = varString('Spearman Manning')
                     DepositionClosure => SpearmanManningHinderedSettling
                  case default
                     RunParams%DepositionChoice = varString(DepositionClosure_d)
                     DepositionClosure => DepositionClosure_dfunc
               end select

            case ('morphodynamic damping')
               set_MorphoDamp=.true.
               morpho_damp = ParamValues(J)%to_lower()
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

            case ('iverson', 'geometric factors')
               set_geometric_factors=.true.
               geometric_factors = ParamValues(J)%to_lower()
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
                     set_geometric_factors = geometric_factors_d
                     GeometricCorrectionFactor => IversonOuyangGeometricCorrectionFactor
                     GeometricCorrectionFactor_gradin => IversonOuyangGeometricCorrectionFactor_gradin
               end select

            case ('curvature')
               curvature = ParamValues(J)%to_lower()
               select case(curvature%s)
                  case ('off')
                     set_curvature=.true.
                     RunParams%curvature=.false.
                  case ('on')
                     set_curvature=.true.
                     RunParams%curvature=.true.
                  case default
                     set_curvature=.false.
               end select
                
            case ('g')
               set_g=.TRUE.
               RunParams%g = ParamValues(J)%to_real()
               if (RunParams%g<=0) call FatalError_Positive(RunParams%InputFile%s,"g")
            
            case ('drag')
               set_drag=.TRUE.
               DragChoice = ParamValues(J)%to_lower()
               select case (DragChoice%s)
                  case ('chezy')
                     RunParams%DragChoice = varString('Chezy')
                     DragClosure => ChezyDrag
                  case ('coulomb')
                     RunParams%DragChoice = varString('Coulomb')
                     DragClosure => CoulombDrag
                  case ('voellmy')
                     RunParams%DragChoice = varString('Voellmy')
                     DragClosure => VoellmyDrag
                  case ('pouliquen')
                     RunParams%DragChoice = varString('Pouliquen')
                     DragClosure => PouliquenDrag
                  case ('variable')
                     RunParams%DragChoice = varString('Variable')
                     DragClosure => VariableDrag
                  case ('manning')
                     RunParams%DragChoice = varString('Manning')
                     DragClosure => ManningDrag
                  case default
                     call FatalErrorMessage("In the 'Parameters' block in the input file " // trim(RunParams%InputFile%s) // new_line('A') &
                        // "The block value for the variable 'Drag' " // DragChoice%s // " is not recognized." // new_line('A') &
                        // "Currently accepted 'Drag' values are 'Chezy', 'Coulomb', 'Voellmy', 'Pouliquen', 'Manning' and 'Variable'")
               end select

            case ('chezy co')
               set_chezyco=.TRUE.
               RunParams%ChezyCo = ParamValues(J)%to_real()
               if (RunParams%ChezyCo<=0) call FatalError_Positive(RunParams%InputFile%s,"Chezy Co")
          
            case ('manning co')
               set_manningco=.TRUE.
               RunParams%ManningCo = ParamValues(J)%to_real()
               if (RunParams%ManningCo<=0) call FatalError_Positive(RunParams%InputFile%s,"Manning Co")
            
            case ('coulomb co')
               set_coulombco=.TRUE.
               RunParams%CoulombCo = ParamValues(J)%to_real()
               if (RunParams%CoulombCo<=0) call FatalError_Positive(RunParams%InputFile%s,"Coulomb Co")
          
            case ('pouliquen min')
               set_pouliquenMinSlope=.TRUE.
               RunParams%PouliquenMinSlope = ParamValues(J)%to_real()
               if (RunParams%PouliquenMinSlope<=0) call FatalError_Positive(RunParams%InputFile%s,"Pouliquen Min")
          
            case ('pouliquen max')
               set_pouliquenMaxSlope=.TRUE.
               RunParams%PouliquenMaxSlope = ParamValues(J)%to_real()
               if (RunParams%PouliquenMaxSlope<=0) call FatalError_Positive(RunParams%InputFile%s,"Pouliquen Max")
          
            case ('pouliquen beta')
               set_pouliquen_beta=.TRUE.
               RunParams%PouliquenBeta = ParamValues(J)%to_real()
               if (RunParams%PouliquenBeta<=0) call FatalError_Positive(RunParams%InputFile%s,"Pouliquen beta")
          
            case ('voellmy switch rate')
               set_VoellmySwitchRate=.TRUE.
               RunParams%VoellmySwitchRate = ParamValues(J)%to_real()
               if (RunParams%VoellmySwitchRate<=0) call FatalError_Positive(RunParams%InputFile%s,"Voellmy Switch Rate")
          
            case ('voellmy switch value')
               set_VoellmySwitchValue=.TRUE.
               RunParams%VoellmySwitchValue = ParamValues(J)%to_real()
               if (RunParams%VoellmySwitchValue<0) call FatalError_Positive(RunParams%InputFile%s,"Voellmy Switch Value")
               if (RunParams%VoellmySwitchValue>1) call FatalError_MaxValue(RunParams%InputFile%s,"Voellmy Switch Value",1)
          
            case ('erosion')
               set_erosion=.TRUE.
               ErosionChoice = ParamValues(J)%to_lower()
               select case (ErosionChoice%s)
                  case ('mixed', 'on')
                     RunParams%ErosionChoice = varString('Mixed')
                     RunParams%MorphodynamicsOn = .TRUE.
                     ErosionClosure => MixedErosion
                  case ('simple')
                     RunParams%ErosionChoice = varString('Simple')
                     RunParams%MorphodynamicsOn = .TRUE.
                     ErosionClosure => SimpleErosion
                  case ('fluid')
                     RunParams%ErosionChoice = varString('Fluid')
                     RunParams%MorphodynamicsOn = .TRUE.
                     ErosionClosure => FluidErosion
                  case ('granular')
                     RunParams%ErosionChoice = varString('Granular')
                     RunParams%MorphodynamicsOn = .TRUE.
                     ErosionClosure => GranularErosion
                  case ('off')
                     RunParams%ErosionChoice = varString('Off')
                     RunParams%MorphodynamicsOn = .FALSE.
                     ErosionClosure => NoErosion
                  case default
                     call FatalErrorMessage("In the 'Parameters' block in the input file " // trim(RunParams%InputFile%s) // new_line('A') &
                        // "The block value for the variable 'Erosion' is not recognized." // new_line('A') &
                        // "Currently accepted 'Erosion' values are 'On, Off, Mixed, Simple, Fluid, Granular'.")
               end select
          
            case ('erosion rate')
               set_EroRate=.TRUE.
               RunParams%EroRate = ParamValues(J)%to_real()
               if (RunParams%EroRate<=0) call FatalError_Positive(RunParams%InputFile%s,"Erosion Rate")
          
            case ('granular erosion rate')
               set_EroRateGranular=.TRUE.
               RunParams%EroRateGranular = ParamValues(J)%to_real()
               if (RunParams%EroRateGranular<=0) call FatalError_Positive(RunParams%InputFile%s,"Granular Erosion Rate")
          
            case ('erosion depth')
               set_EroDepth=.TRUE.
               RunParams%EroDepth = ParamValues(J)%to_real()
               if (RunParams%EroDepth<0) call FatalError_NonNegative(RunParams%InputFile%s,"Erosion Depth")

            case ('erosion transition')
               set_EroTransition=.TRUE.
               eroTransition = ParamValues(J)%to_lower()
               select case (eroTransition%s)
                  case ('smooth')
                     RunParams%ErosionTransition = varString('smooth')
                     ErosionTransition => SmoothErosionTransition
                  case ('step')
                     RunParams%ErosionTransition = varString('step')
                     ErosionTransition => StepErosionTransition
                  case ('off')
                     RunParams%ErosionTransition = varString('off')
                     ErosionTransition => NoErosionTransition
                  case default
                     RunParams%ErosionTransition = varString(EroTransition_d)
                     ErosionTransition => EroTransition_dfunc
               end select
         
            case ('erosion critical height')
               set_EroCriticalHeight=.TRUE.
               RunParams%EroCriticalHeight = ParamValues(J)%to_real()
               if (RunParams%EroCriticalHeight<=0) call FatalError_Positive(RunParams%InputFile%s,"Erosion Critical Height")
          
            case ('bed porosity')
               set_BedPorosity=.TRUE.
               RunParams%BedPorosity = ParamValues(J)%to_real()
               if (RunParams%BedPorosity<0) call FatalError_NonNegative(RunParams%InputFile%s,"Bed Porosity")
               if (RunParams%BedPorosity>=1) call FatalError_MaxValue(RunParams%InputFile%s,"Bed Porosity",1)
          
            case ('rhow')
               set_rhow=.TRUE.
               RunParams%rhow = ParamValues(J)%to_real()
               if (RunParams%rhow<=0) call FatalError_Positive(RunParams%InputFile%s,"rhow")
          
            case ('rhos')
               set_rhos=.TRUE.
               RunParams%rhos = ParamValues(J)%to_real()
               if (RunParams%rhos<=0) call FatalError_Positive(RunParams%InputFile%s,"rhos")
          
            case ('maxpack', 'max pack')
               set_maxPack=.TRUE.
               RunParams%maxPack = ParamValues(J)%to_real()
               if (RunParams%maxPack<=0) call FatalError_Positive(RunParams%InputFile%s,"maxPack")
               if (RunParams%maxPack>1) call FatalError_MaxValue(RunParams%InputFile%s,"maxPack",1)
          
            case ('solid diameter')
               set_SolidDiameter=.TRUE.
               RunParams%SolidDiameter = ParamValues(J)%to_real()
               if (RunParams%SolidDiameter<=0) call FatalError_Positive(RunParams%InputFile%s,"Solid Diameter")
          
            case ('eddy viscosity')
               set_EddyViscosity=.TRUE.
               RunParams%EddyViscosity = ParamValues(J)%to_real()
               if (RunParams%EddyViscosity<0) call FatalError_NonNegative(RunParams%InputFile%s,"Eddy Viscosity")
          
            case ('settling speed')
               set_SettlingSpeed=.TRUE.
               RunParams%ws0 = ParamValues(J)%to_real()
          
            case default
               call InputLabelUnrecognized(ParamLabels(J)%s)

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

      if (.not.set_curvature) then
         RunParams%curvature=curvature_d
      end if

      if (.not.set_g) RunParams%g = g_d
      if (.not.set_drag) then
         RunParams%DragChoice = varString('Chezy')
         DragClosure => ChezyDrag
         call WarningMessage("In the 'Parameters' block 'Drag' is not given.  Using default 'Chezy' drag.")
      end if

      if ((RunParams%DragChoice%s=='Chezy').and.(.not.set_chezyco)) then
         RunParams%ChezyCo = chezyco_d
         call Warning_DragDefaultValue("Chezy","Chezy Co",RunParams%ChezyCo)
      end if

      if ((RunParams%DragChoice%s=="Coulomb").and.(.not.set_coulombco)) then
         RunParams%CoulombCo = coulombco_d
         call Warning_DragDefaultValue("Coulomb","Coulomb Co",RunParams%CoulombCo)
      end if

      if (RunParams%DragChoice%s=="Voellmy") then
         if (.not.set_chezyco) then
            RunParams%ChezyCo = chezyco_d
            call Warning_DragDefaultValue("Voellmy","Chezy Co",RunParams%ChezyCo)
         end if
         if (.not.set_coulombco) then
            RunParams%CoulombCo = coulombco_d
            call Warning_DragDefaultValue("Voellmy","Coulomb Co",RunParams%CoulombCo)
         end if
      end if

      if ((RunParams%DragChoice%s=="Pouliquen").or.(RunParams%DragChoice%s=="Variable")) then
         if (.not.set_pouliquenMinSlope) then
            RunParams%PouliquenMinSlope = pouliquenMinSlope_d
            call Warning_DragDefaultValue(RunParams%DragChoice%s,"Pouliquen Min",RunParams%pouliquenMinSlope)
         end if
         if (.not.set_pouliquenMaxSlope) then
            RunParams%PouliquenMaxSlope = pouliquenMaxSlope_d
            call Warning_DragDefaultValue(RunParams%DragChoice%s,"Pouliquen Max",RunParams%pouliquenMaxSlope)
         end if
         if (.not.set_pouliquen_beta) then
            RunParams%PouliquenBeta = pouliquen_beta_d
            call Warning_DragDefaultValue(RunParams%DragChoice%s,"Pouliquen Beta",RunParams%PouliquenBeta)
         end if
      end if

      if (RunParams%DragChoice%s=="Variable") then
         if (.not.set_VoellmySwitchRate) then
            RunParams%VoellmySwitchRate = VoellmySwitchRate_d
            call Warning_DragDefaultValue(RunParams%DragChoice%s,"Voellmy switch rate",RunParams%VoellmySwitchRate)
         end if
         if (.not.set_VoellmySwitchValue) then
            RunParams%VoellmySwitchValue = VoellmySwitchValue_d
            call Warning_DragDefaultValue(RunParams%DragChoice%s,"Voellmy switch value",RunParams%VoellmySwitchValue)
         end if
      end if

      if ((RunParams%DragChoice%s=="Manning").and.(.not.set_manningco)) then
         RunParams%ManningCo = manningco_d
         call Warning_DragDefaultValue("Manning","Manning Co",RunParams%ManningCo)
      end if

      if (.not.set_deposition) then
         RunParams%DepositionChoice = varString(DepositionClosure_d)
         DepositionClosure => DepositionClosure_dfunc
      end if

      if (.not.set_erosion) then
         RunParams%ErosionChoice = varString(ErosionClosure_d)
         RunParams%MorphodynamicsOn = .TRUE.
         ErosionClosure => ErosionClosure_dfunc
      end if
      
      if (((RunParams%ErosionChoice%s=='Fluid').or.(RunParams%ErosionChoice%s=='Mixed')).and.(.not.set_EroRate)) then
         RunParams%EroRate = EroRate_d
         call Warning_ErosionDefaultValue("Erosion Rate",RunParams%EroRate)
      end if
      if (((RunParams%ErosionChoice%s=='Granular').or.(RunParams%ErosionChoice%s=='Mixed')).and.(.not.set_EroRateGranular)) then
         RunParams%EroRateGranular = EroRateGranular_d
         call Warning_ErosionDefaultValue("Granular Erosion Rate",RunParams%EroRateGranular)
      end if
      if ((RunParams%ErosionChoice%s/="Off").and.(.not.set_EroDepth)) then
         RunParams%EroDepth = EroDepth_d
         call Warning_ErosionDefaultValue("Erosion Depth",RunParams%EroDepth)
      end if
      if ((RunParams%ErosionChoice%s/="Off").and.(.not.set_EroCriticalHeight)) then
         RunParams%EroCriticalHeight = EroCriticalHeight_d
         call Warning_ErosionDefaultValue("Erosion Critical Height",RunParams%EroCriticalHeight)
      end if

      if (.not.set_BedPorosity) then
         RunParams%BedPorosity = BedPorosity_d
         call Warning_DefaultValue("BedPorosity",RunParams%BedPorosity)
      end if

      if (.not.set_rhow) then
         RunParams%rhow = rhow_d
         call Warning_DefaultValue("rhow",RunParams%rhow)
      end if

      if (.not.set_rhos) then
         RunParams%rhos = rhos_d
         call Warning_DefaultValue("rhos",RunParams%rhos)
      end if

      if (.not.set_maxPack) then
         RunParams%maxPack = maxPack_d
         call Warning_DefaultValue("maxPack",RunParams%maxPack)
      end if

      if (.not.set_SolidDiameter) then
         RunParams%SolidDiameter = SolidDiameter_d
         call Warning_DefaultValue("Solid Diameter",RunParams%SolidDiameter)
      end if

      if (.not.set_EddyViscosity) then
         RunParams%EddyViscosity = EddyViscosity_d
         call Warning_DefaultValue("Eddy Viscosity",RunParams%EddyViscosity)
      end if

      RunParams%gred = (RunParams%rhos/RunParams%rhow-1.0_wp)*RunParams%g

      RunParams%Rep = sqrt(RunParams%g*RunParams%SolidDiameter)*RunParams%SolidDiameter/visc_w ! Particle Reynolds number

      R = (RunParams%gred/visc_w/visc_w)**(1.0_wp/3.0_wp)*RunParams%SolidDiameter ! Scaled grain size

      if (.not.set_SettlingSpeed) then
         RunParams%ws0 = visc_w/RunParams%SolidDiameter*(sqrt(10.36_wp*10.36_wp + 1.048_wp*R*R*R) - 10.36_wp) ! Clear water settling velocity
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
