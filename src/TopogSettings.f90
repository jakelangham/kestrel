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
! with the topography, as read from the "Topog" block of an input file.
! Some input validation is performed, producing Warning or FatalError messages.
module topog_settings_module

   use iso_c_binding, only: c_bool
   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage, InputLabelUnrecognized, WarningMessage
   use varstring_module, only: varString
   use utilities_module, only: PathTrail, Int2String
   use runsettings_module, only: RunSet
   use topog_funcs_module
   use utm_module

   implicit none

   private
   public :: Topog_Set

   ! default values
   logical(kind=c_bool), parameter :: EmbedRaster_d = .FALSE.
   logical, parameter :: RebuildDEM_d = .TRUE.

contains

   subroutine Topog_Set(TopogLabels,TopogValues,RunParams)
      ! Set topog settings from input file.
      ! Inputs: TopogLabels [character string array] - Labels for Topog;
      !         TopogValues [character string array] - The values associated with the labels.
      ! Output: RunParams [type RunSet] - The solver settings.

      implicit none

      type(varString), dimension(:), intent(in) :: TopogLabels
      type(varString), dimension(:), intent(in) :: TopogValues
      type(RunSet), intent(inout) :: RunParams

      type(varString) :: label
      type(varString) :: DemPath
      type(varString) :: RasterFile
      type(varString) :: Type
      type(varString) :: TopogFuncStr
      type(varString) :: SRTMPath
      type(varString) :: EmbedRaster
      type(varString) :: RebuildDEM

      integer :: J, N

      logical :: set_DemPath
      logical :: set_SRTMPath
      logical :: set_Type
      logical :: set_TopogFunc
      logical :: set_Topog_params
      logical :: set_RasterFile
      logical :: set_EmbedRaster
      logical :: set_RebuildDEM

      character(len=4096) :: cwd

      N = size(TopogValues)

      set_Type=.FALSE.
      set_DemPath=.FALSE.
      set_SRTMPath=.FALSE.
      set_RasterFile=.FALSE.
      set_TopogFunc=.FALSE.
      set_Topog_params=.FALSE.
      set_EmbedRaster=.FALSE.
      set_RebuildDEM=.FALSE.

      do J=1,N
         label = TopogLabels(J)%to_lower()
         select case (label%s)
            case ('type')
               set_Type=.TRUE.
               Type = TopogValues(J)%to_lower()
            case ('dem directory')
               set_DemPath=.TRUE.
               DemPath = TopogValues(J)
            case ('srtm directory')
               set_SRTMPath=.TRUE.
               SRTMPath = TopogValues(J)
            case ('raster file')
               set_RasterFile=.TRUE.
               RasterFile = TopogValues(J)
            case ('topog function')
               set_TopogFunc=.TRUE.
               TopogFuncStr = TopogValues(J)%to_lower()
            case ('topog params')
               set_Topog_params=.TRUE.
               call TopogValues(J)%read_set(RunParams%TopogFuncParams)
            case ('embed raster')
               set_EmbedRaster=.TRUE.
               EmbedRaster = TopogValues(J)%to_lower()
            case ('rebuild dem')
               set_RebuildDEM=.TRUE.
               RebuildDEM = TopogValues(J)%to_lower()
            case default
               call InputLabelUnrecognized(TopogLabels(J)%s)
         end select
      end do

      ! Check required variable Type is set
      if (.not.set_Type) then
         call FatalErrorMessage("In the 'Topog' block in the input file "// trim(RunParams%InputFile%s) // new_line('A') &
                  // "the required variable 'Type' is not given.")
      end if
      select case (Type%s)
         case ('dem', 'raster')
            RunParams%Topog = varString("raster")
         case ('srtm')
            RunParams%Topog = varString("srtm")
         case ('function')
            RunParams%Topog = varString("function")
         case default
            call FatalErrorMessage("In the 'Topog' block in the input file "// trim(RunParams%InputFile%s) // new_line('A') &
                  // "the required variable 'Type = "//Type%s//"' is not recognized.")
      end select

      ! Check conditionally required/optional settings for Type = DEM/Raster
      if (RunParams%Topog == "raster") then

         ! Get required RasterFile
         if (.not.set_RasterFile) then
            call FatalErrorMessage("In the 'Topog' block in the input file "// trim(RunParams%InputFile%s) // new_line('A') &
                  // " 'Raster File' is required when 'Type' is 'Raster' or 'DEM'.")
         end if
         RunParams%RasterFile = RasterFile
      
         ! Check optional, setting defaults if not given

         ! Get DemPath as CWD if not set
         if (.not. set_DemPath) then
            call getcwd(cwd)
            DemPath = varString(cwd, trim_str=.TRUE.)
            call WarningMessage("In the 'Topog' block 'dem directory' is not given.  Using current directory")
         end if   
         RunParams%DemPath = PathTrail(DemPath)

         ! Check EmbedRaster is set
         if (.not.set_EmbedRaster) then
            RunParams%EmbedRaster = EmbedRaster_d
         else
            select case (EmbedRaster%s)
               case ('on','yes')
                  RunParams%EmbedRaster = .TRUE.
               case ('off','no')
                  RunParams%EmbedRaster = .FALSE.
               case default
                  RunParams%EmbedRaster = EmbedRaster_d
            end select
         end if
         
         ! If Embed raster, then also need SRTM
         if (RunParams%EmbedRaster) then
            if (.not. set_SRTMPath) then
               call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " 'SRTM directory' is not given.")
            end if
            RunParams%SRTMPath = PathTrail(SRTMPath)
         end if

         RunParams%Georeference = .TRUE.
      end if
         
      ! Check conditionally required settings for Type = SRTM
      if (RunParams%Topog == "srtm") then
         
         ! Require SRTM path
         if (.not. set_SRTMPath) then
            call FatalErrorMessage("In the 'Topog' block in the input file " // &
               trim(RunParams%InputFile%s) // new_line('A') &
               // " 'SRTM directory' is not given.")
         end if
         RunParams%SRTMPath = PathTrail(SRTMPath)

         RunParams%Georeference = .TRUE.
      end if

      ! Check conditionally required/optional settings for Type = Function
      if (RunParams%Topog == "function") then

         ! Require TopogFunc
         if (.not. set_TopogFunc) call FatalErrorMessage("In the 'Topog' block in the input file "// trim(RunParams%InputFile%s) // new_line('A') &
            // " missing variable 'Topog function' that is required when 'Type' is 'Function'.")

         ! Check if the number of TopogFuncParams are appropriate for each TopogFunc
         select case (TopogFuncStr%s)

            case ('flat')
               TopogFunc => flat

            case ('xslope')
               if (size(RunParams%TopogFuncParams) < 1) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " one 'Topog params' is needed for 'Type = xslope'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xslope

            case ('yslope')
               if (size(RunParams%TopogFuncParams) < 1) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " one 'Topog params' is needed for 'Type = yslope'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => yslope
            
            case ('xyslope')
               if (size(RunParams%TopogFuncParams) < 2) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " two 'Topog params' are needed for 'Type = xyslope'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xyslope

            case ('x2slopes')
               if (size(RunParams%TopogFuncParams) < 3) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                     trim(RunParams%InputFile%s) // new_line('A') &
                     // " three 'Topog params' are needed for 'Type = x2slopes'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => x2Slopes

            case ('usgs')
               if (size(RunParams%TopogFuncParams) < 2) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                     trim(RunParams%InputFile%s) // new_line('A') &
                     // " two 'Topog params' are needed for 'Type = USGS'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => usgs
            
            case ('xbislope')
               if (size(RunParams%TopogFuncParams) < 3) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                     trim(RunParams%InputFile%s) // new_line('A') &
                     // " three 'Topog params' are needed for 'Type = xbislope'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xbislope            

            case ('xsinslope')
               if (size(RunParams%TopogFuncParams) < 1) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " one 'Topog params' is needed for 'Type = xsinslope'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xsinslope

            case ('xysinslope')
               if (size(RunParams%TopogFuncParams) < 1) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " one 'Topog params' is needed for 'Type = xysinslope'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xysinslope

            case ('xhump')
               if (size(RunParams%TopogFuncParams) < 2) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " two 'Topog params' are needed for 'Type = xhump'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xhump

            case ('xtanh')
               if (size(RunParams%TopogFuncParams) < 3) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " three 'Topog params' are needed for 'Type = xtanh'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xtanh
            
            case ('xparab')
               if (size(RunParams%TopogFuncParams) < 1) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " one 'Topog params' is needed for 'Type = xparab'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xparab
            
            case ('xyparab')
               if (size(RunParams%TopogFuncParams) < 2) then
                  call FatalErrorMessage("In the 'Topog' block in the input file " // &
                  trim(RunParams%InputFile%s) // new_line('A') &
                  // " two 'Topog params' are needed for 'Type = xyparab'; received " // Int2String(size(RunParams%TopogFuncParams)))
               end if
               TopogFunc => xyparab

            case default
               call FatalErrorMessage("In the 'Topog' block in the input file "// trim(RunParams%InputFile%s) // new_line('A') &
                  // " the 'Topog function = "// TopogFuncStr%s //"' is not recognized.")

         end select
         RunParams%Georeference = .FALSE.
      end if

      if (.not.set_RebuildDEM) then
         RunParams%RebuildDEM = RebuildDEM_d
      else
         label = RebuildDEM%to_lower()
         select case (label%s)
            case ('yes')
               RunParams%RebuildDEM = .TRUE.
            case ('no')
               RunParams%RebuildDEM = .FALSE.
            case default
               RunParams%RebuildDEM = RebuildDEM_d
         end select
      end if

      if ((RunParams%Topog%s=='raster').and.(.not.set_EmbedRaster)) then
         if (EmbedRaster_d) then
            call WarningMessage("In the 'Topog' block Type = Raster and 'Embed Raster' is not set." // new_line('A') &
               // " The raster file will be embedded into the background 30m SRTM DEM." // new_line('A') &
               // " Mismatches in the elevation of the raster and SRTM DEMs are not corrected.")
         else
            call WarningMessage("In the 'Topog' block Type = Raster and 'Embed Raster' is not set." // new_line('A') &
               // " The raster file will be used on its own." // new_line('A') &
               // " Flow off the raster is not modelled.")
         end if
      end if

      if (RunParams%Georeference) then
         RunParams%UTM_zone_number = LatLonToZoneNumber(RunParams%Lat, RunParams%Lon)
         if (RunParams%Lat<0.0_wp) then
            RunParams%UTM_zone_letter = "S"
         else
            RunParams%UTM_zone_letter = "N"
         end if
         RunParams%utmEPSG = LatLonToUtmEpsg(RunParams%Lat, RunParams%Lon)
         RunParams%projTransformer = proj_transformer(RunParams%utmEPSG)
         RunParams%centerUTM = RunParams%projTransformer%wgs84_to_utm(RunParams%Lat, RunParams%Lon)
      end if

   end subroutine Topog_Set

end module topog_settings_module
