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


! This module assigns values to variables in RunSet associated
! with the simulation outputs, as read from the "Output" block of an input file.
! Some input validation is performed, producing Warning or FatalError messages.
module outputsettings_module
   
   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage, InputLabelUnrecognized, WarningMessage
   use varstring_module, only: varString
   use utilities_module, only: CheckPath
   use runsettings_module, only: RunSet

   implicit none

   private
   public :: Output_Set

! default values
   character(len=7), parameter :: OutDir_d = "results"
   character(len=11), parameter :: InfoFilename_d = "RunInfo.txt"
   character(len=10), parameter :: MaxHeightFilename_d = "MaxHeights"
   character(len=9), parameter :: MaxSpeedFilename_d = "MaxSpeeds"
   character(len=10), parameter :: MaxErosionFilename_d = "MaxErosion"
   character(len=10), parameter :: MaxDepositFilename_d = "MaxDeposit"
   character(len=14), parameter :: InundationTimeFilename_d = "InundationTime"
   character(len=8), parameter :: MaximumsFilename_d = "Maximums"
   logical, parameter :: compressOutput_d = .false.
   integer, parameter :: Nout_d = 100 ! number of output files
   real(kind=wp), parameter :: kmlHeight_d = 0.1_wp ! Threshold for depths in kml files

contains

   subroutine Output_Set(OutputLabels,OutputValues,RunParams)
      ! Set output settings from input file.
      ! Inputs: OutputLabels [type:varString] - Labels for Output;
      !         OutputValues [type:varString] - The values associated with the labels.
      ! Output: RunParams [type RunSet] - The solver settings.

      implicit none

      type(varString), dimension(:), intent(in) :: OutputLabels
      type(varString), dimension(:), intent(in) :: OutputValues
      type(RunSet), intent(inout) :: RunParams

      type(varString) :: label
      type(varString) :: labelValue
      type(varString) :: OutDir
      type(varString) :: basePath
      type(varString) :: InfoFilename
      type(varString) :: MaxHeightFilename
      type(varString) :: MaxSpeedFilename
      type(varString) :: MaxErosionFilename
      type(varString) :: MaxDepositFilename
      type(varString) :: MaximumsFilename
      type(varString) :: InundationTimeFilename
      type(varString), dimension(:), allocatable :: outputFormat

      integer :: J, N, K

      logical :: set_OutDir
      logical :: set_basePath
      logical :: set_InfoFilename
      logical :: set_MaxHeightFilename
      logical :: set_MaxSpeedFilename
      logical :: set_MaxErosionFilename
      logical :: set_MaxDepositFilename
      logical :: set_InundationFilename
      logical :: set_InundationTimeFilename
      logical :: set_MaximumsFilename
      logical :: set_Nout
      logical :: set_kmlHeight
      logical :: set_compressOutput
      logical :: set_outputFormat

      integer :: N_fmts

      character(len=4096) :: cwd

      N = size(OutputValues)

      set_OutDir=.FALSE.
      set_basePath=.FALSE.
      set_InfoFilename=.FALSE.
      set_MaxHeightFilename=.FALSE.
      set_MaxSpeedFilename=.FALSE.
      set_MaxErosionFilename=.FALSE.
      set_MaxDepositFilename=.FALSE.
      set_InundationFilename=.FALSE.
      set_InundationTimeFilename=.FALSE.
      set_MaximumsFilename=.FALSE.
      set_Nout=.FALSE.
      set_kmlHeight=.FALSE.
      set_compressOutput=.FALSE.
      set_outputFormat=.FALSE.

      do J=1,N
        label = OutputLabels(J)%to_lower()
         select case (label%s)
          case ('directory')
            set_OutDir=.TRUE.
            OutDir = OutputValues(J)
          case ('base path')
            set_basePath=.TRUE.
            basePath = OutputValues(J)
          case ('info filename')
            set_InfoFilename=.TRUE.
            InfoFilename = OutputValues(J)
          case ('inundation time filename')
            set_InundationTimeFilename=.TRUE.
            InundationTimeFilename = OutputValues(J)
          case ('max height filename')
            set_MaxHeightFilename=.TRUE.
            MaxHeightFilename = OutputValues(J)
          case ('max speed filename')
            set_MaxSpeedFilename=.TRUE.
            MaxSpeedFilename = OutputValues(J)
          case ('max erosion filename')
            set_MaxErosionFilename=.TRUE.
            MaxErosionFilename = OutputValues(J)
          case ('max deposit filename')
            set_MaxDepositFilename=.TRUE.
            MaxDepositFilename = OutputValues(J)
          case ('maximums filename')
            set_MaximumsFilename=.TRUE.
            MaximumsFilename = OutputValues(J)
          case ('n out')
            set_Nout=.TRUE.
            RunParams%Nout = OutputValues(J)%to_int()
            if (RunParams%Nout<1) call FatalErrorMessage("In the 'Output' block in the input file "// trim(RunParams%InputFile%s) // new_line('A') &
               // " The block variable 'N out' must be positive.")
          case ('kml height')
            set_kmlHeight=.TRUE.
            RunParams%kmlHeight = OutputValues(J)%to_real()
            if (RunParams%kmlHeight.le.RunParams%heightThreshold) then
               call FatalErrorMessage("In the 'Output' block in the input file "// trim(RunParams%InputFile%s) // new_line('A') &
                  // " The block variable 'kml height' must be greater than ", RunParams%heightThreshold)
            end if

          case ('format')
            set_outputFormat=.TRUE.
            call OutputValues(J)%read_list(outputFormat, delimiter=',')
            N_fmts = OutputValues(J)%count_substring(',')+1
            do K=1,N_fmts
               select case (outputFormat(K)%s)
                case ('txt')
                  RunParams%out_txt = .TRUE.
                case ('kml')
                  RunParams%out_kml = .TRUE.
                case ('nc','netcdf')
#if HAVE_NETCDF4
                  RunParams%out_nc = .TRUE.
#else
                  call FatalErrorMessage('NetCDF not active')
#endif
                case default
                  call WarningMessage("In the 'Output' block the value of 'Format' is not recognised.")
               end select
            end do

          case ('compression', 'compress')
            set_compressOutput=.TRUE.
            labelValue = OutputValues(J)%to_lower()
            select case (labelValue%s)
             case ('on')
               RunParams%compressOutput = .true.
             case ('off')
               RunParams%compressOutput = .false.
             case default
               call WarningMessage("In the 'Output' block the value of 'Compression' is not recognised.")
               RunParams%compressOutput = .false.
            end select

          case default
            call InputLabelUnrecognized(OutputLabels(J)%s)
         end select
      end do

      if (.not.set_basePath) then
         call getcwd(cwd)
         basePath = varString(cwd, trim_str=.TRUE.)
         basePath = CheckPath(basePath)
         call WarningMessage("In the 'Output' block 'Base path' is not given.  Using current working directory  " // basePath%s)
      end if
      RunParams%basePath = basePath

      if (.not.set_OutDir) then
         OutDir = varString(OutDir_d)
         call WarningMessage("In the 'Output' block 'Directory' is not given.  Using default directory " // OutDir_d)
      end if
      RunParams%OutDir = RunParams%basePath + OutDir
      RunParams%OutDir = CheckPath(RunParams%OutDir)

      if (.not.set_InfoFilename) then
         InfoFilename = varString(InfoFilename_d)
      end if
      RunParams%InfoFilename = InfoFilename

      if (.not.set_Nout) then
         call FatalErrorMessage("In the 'Output' block in the input file "// trim(RunParams%InputFile%s) // new_line('A') &
            // " 'N out' is not given.")
      end if

      if (.not.set_outputFormat) then
         RunParams%out_txt = .TRUE.
         call WarningMessage("In the 'Output' block in the input file 'Format' is not given.  Using txt output format.")
      end if

      ! Set defaults for conditionally optional settings for txt output if not set.
      if (RunParams%out_txt) then

         if (.not. set_compressOutput) then
            RunParams%CompressOutput = compressOutput_d
            set_compressOutput = .TRUE.
         end if

         if (.not.set_MaxHeightFilename) then
            MaxHeightFilename = varString(MaxHeightFilename_d)
            set_MaxHeightFilename = .TRUE.
         end if
         RunParams%MaxHeightFilename = MaxHeightFilename
         
         if (.not.set_MaxSpeedFilename) then
            MaxSpeedFilename = varString(MaxSpeedFilename_d)
            set_MaxSpeedFilename = .TRUE.
         end if
         RunParams%MaxSpeedFilename = MaxSpeedFilename
         
         if (.not.set_MaxErosionFilename) then
            MaxErosionFilename = varString(MaxErosionFilename_d)
         end if
         RunParams%MaxErosionFilename = MaxErosionFilename

         if (.not.set_MaxDepositFilename) then
            MaxDepositFilename = varString(MaxDepositFilename_d)
         end if
         RunParams%MaxDepositFilename = MaxDepositFilename

         if (.not.set_InundationTimeFilename) then
            InundationTimeFilename = varString(InundationTimeFilename_d)
            set_InundationTimeFilename = .TRUE.
         end if
         RunParams%InundationTimeFilename = InundationTimeFilename
      end if

      ! Set defaults for conditionally optional settings for nc output if not set.
      if (RunParams%out_nc) then
         if (.not.set_MaximumsFilename) then
            MaximumsFilename = varString(MaximumsFilename_d)
         end if
         RunParams%MaximumsFilename = MaximumsFilename
      end if

      ! Set defaults for conditionally optional settings for txt output if not set.
      if (RunParams%out_kml) then
         if (.not.set_kmlHeight) then
            RunParams%kmlHeight = kmlHeight_d
         end if

         if (.not.set_MaxHeightFilename) then
            MaxHeightFilename = varString(MaxHeightFilename_d)
            set_MaxHeightFilename = .TRUE.
         end if
         RunParams%MaxHeightFilename = MaxHeightFilename
         
         if (.not.set_MaxSpeedFilename) then
            MaxSpeedFilename = varString(MaxSpeedFilename_d)
            set_MaxSpeedFilename = .TRUE.
         end if
         RunParams%MaxSpeedFilename = MaxSpeedFilename

         if (.not.set_InundationTimeFilename) then
            InundationTimeFilename = varString(InundationTimeFilename_d)
            set_InundationTimeFilename = .TRUE.
         end if
         RunParams%InundationTimeFilename = InundationTimeFilename
      end if

      RunParams%DeltaT = (RunParams%tend - RunParams%tstart) / real(RunParams%Nout,kind=wp)

   end subroutine Output_Set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module outputsettings_module
