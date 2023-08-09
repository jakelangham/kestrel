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


! This module contains routines that are used to output simulation
! results to file.
!
! Results can be output as delimiter text, NetCDF and KML files.
! The output choice is set in the input file via OutputSettings.f90.
!
! NetCDF is an optional feature, and requires netCDF libraries.
! If configuration does not specify HAVE_NETCDF then the NetCDF output
! is not allowed.
!
! The main public routines are:
! - OutputSolutionData -- outputs timestep solution data in required formats
! - OutputAggregateData -- outputs aggregated solution data (e.g. maximum over time of fields) in required formats
! - OutputInfo -- outputs information of simulation
! - CalculateVolume -- computes and outputs total volume and related information
!
module output_module

   use set_precision_module, only: wp
   use utilities_module, only: CheckFileExists, CheckPath, Int2String, KahanAdd, pair
   use grid_module, only: GridCoords, GridType, TileList, TileType
   use closures_module, only : FlowSquaredSpeedSlopeAligned, GeometricCorrectionFactor
   use runsettings_module, only : RunSet
   use varstring_module, only: varString
   use messages_module, only: InfoMessage, WarningMessage
   use morphodynamic_rhs_module, only: ComputeCellCentredTopographicData
#if HAVE_NETCDF4
   use netcdf_utils_module
   use netcdf
#endif

   implicit none

   private
   public :: GenerateOutputFilename
   public :: CreateOutDir
   public :: OutputSolutionData
   public :: OutputAggregateData
   public :: OutputInfo
   public :: CalculateVolume

   ! Generate a timestep output filename based on current output time step.
   ! Format is XXXXXX.txt
   ! Overloaded to allow for filename with or without basename
   interface GenerateOutputFilename
      module procedure GenerateOutputFilename_basename, GenerateOutputFilename_nobasename
   end interface GenerateOutputFilename

contains

   ! Generate a timestep output filename based on current output time step.
   ! Inputs: basename - string proceeding filename (possibly with path)
   !         n - integer label
   !         width - width of filename label (i.e. n is left zero padded to specified width)
   !         ext [optional, default='.txt'] - extension to append to filename (should include '.')
   ! Returns: filename - varString with filename in format basename_XXXXXX.ext
   function GenerateOutputFilename_basename(basename,n,width, ext) result(filename)
      type(varString), intent(in) :: basename
      integer, intent(in) :: n
      integer, intent(in) :: width
      character(len=*), optional, intent(in) :: ext
      type(varString) :: filename

      character(len=20) :: iwdth
      character(len=20) :: nstr

      integer :: basename_length

      type(varString) :: ext_in

      if (present(ext)) then
         ext_in = varString(ext)
      else
         ext_in = varString('.txt')
      end if

      write(iwdth,*) width
      write(nstr,"(I" // ADJUSTL(iwdth) // "." // ADJUSTL(iwdth) // ")") n

      basename_length = len_trim(basename%s)

      if (basename_length>0) then
         filename = basename
         if ((basename%s(basename_length-1:)=='/').or.(basename%s(basename_length-1:)=='\\')) then
            filename = filename + "_"
         end if
      end if

      filename = filename + trim(nstr) + ext_in

      return

   end function GenerateOutputFilename_basename


   ! Generate a timestep output filename based on current output time step.
   ! Inputs: n - integer label
   !         width - width of filename label (i.e. n is left zero padded to specified width)
   !         ext [optional, default='.txt'] - extension to append to filename (should include '.')
   ! Returns: filename - varString with filename in format XXXXXX.ext
   function GenerateOutputFilename_nobasename(n,width, ext) result(filename)
      integer, intent(in) :: n
      integer, intent(in) :: width
      character(len=*), optional, intent(in) :: ext
      type(varString) :: filename

      character(len=20) :: iwdth
      character(len=20) :: nstr

      type(varString) :: ext_in

      if (present(ext)) then
         ext_in = varString(ext)
      else
         ext_in = varString('.txt')
      end if

      write(iwdth,*) width
      write(nstr,"(I" // ADJUSTL(iwdth) // "." // ADJUSTL(iwdth) // ")") n

      filename = varString(trim(nstr)) + ext_in

      return

   end function GenerateOutputFilename_nobasename


   ! Create on output directory from user specified OutDir in RunParams
   subroutine CreateOutDir(RunParams)
      type(RunSet), intent(in) :: RunParams

      logical dir_exists

      inquire (file=RunParams%OutDir%s, exist=dir_exists)
      if (dir_exists) then
         call WarningMessage("Directory "//RunParams%OutDir%s//" exists and will be overwritten.")
      else
         ! workaround: it calls an extern program...
         call WarningMessage("Directory "//RunParams%OutDir%s//" will be created.")
         call execute_command_line('mkdir -p '//trim(adjustl(RunParams%OutDir%s)), wait=.true.)
      end if
   end subroutine CreateOutDir


   ! Write solution data at a timestep to text and/or NetCDF files
   ! Inputs: RunParams - run settings
   !         outnum - integer output label
   !         grid - numerical solution data
   subroutine OutputSolutionData(RunParams, outnum, grid)
      type(RunSet), intent(in) :: RunParams
      integer, intent(in) :: outnum
      type(GridType), target, intent(in) :: grid

      type(varString) :: filename

      if (RunParams%out_txt) then
         filename = GenerateOutputFilename(outnum, 6)
         call OutputSolutionData_txt(RunParams, filename, grid)
      end if

#if HAVE_NETCDF4
      if (RunParams%out_nc) then
         filename = GenerateOutputFilename(outnum, 6, ext = '.nc')
         call OutputSolutionData_NetCDF(RunParams, filename%s, grid)
      end if
#endif
      return
   end subroutine OutputSolutionData

   ! Write aggregated solution data to text, NetCDF files, and/or KML file
   ! Inputs: RunParams - run settings
   !         grid - numerical solution data
   subroutine OutputAggregateData(RunParams, grid)
      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(in) :: grid

      type(varString) :: MaximumFilenameTXT
      type(varString) :: MaxHeightFilenameKML
      type(varString) :: MaxSpeedFilenameKML
      type(varString) :: InundationTimeFilenameKML
#if HAVE_NETCDF4
      type(varString) :: MaximumFilenameNetCDF
#endif
      
      if (RunParams%out_txt) then
         MaximumFilenameTXT = RunParams%OutDir + RunParams%MaximumsFilename + '.txt'
         call OutputMaximums_txt(RunParams, MaximumFilenameTXT%s, grid)
      end if

      if (RunParams%out_kml) then
         MaxHeightFilenameKML = RunParams%OutDir + RunParams%MaxHeightFilename + '.kml'
         MaxSpeedFilenameKML = RunParams%OutDir + RunParams%MaxSpeedFilename + '.kml'
         InundationTimeFilenameKML = RunParams%OutDir + RunParams%InundationTimeFilename + '.kml'
         call OutputMaxHeights_kml(RunParams, MaxHeightFilenameKML%s, grid)
         call OutputMaxSpeeds_kml(RunParams, MaxSpeedFilenameKML%s, grid)
         call OutputInundationTime_kml(RunParams, InundationTimeFilenameKML%s, grid)
      end if

#if HAVE_NETCDF4
      if (RunParams%out_nc) then
         MaximumFilenameNetCDF = RunParams%OutDir + RunParams%MaximumsFilename + '.nc'
         call OutputMaximums_NetCDF(RunParams, MaximumFilenameNetCDF%s, grid)
      end if
#endif

   end subroutine OutputAggregateData


   ! Write simulation info file
   ! Inputs: RunParams - run settings
   subroutine OutputInfo(RunParams)
      type(RunSet), intent(in) :: RunParams
      
      type(varString) :: filename

      logical :: FileExists
      integer :: J, K

      type(pair) :: pt
      type(pair) :: distxy

      filename = RunParams%OutDir + RunParams%InfoFilename
      filename = filename%trim()

      inquire (file=filename%s, exist=FileExists)
      if (FileExists) then
         open (101, file=filename%s, status="replace", action="write")
      else
         open (101, file=filename%s, status="new", action="write")
      end if

      write (101, fmt="(a12,a4,a1,a2,a1,a2,a1,a4)") "Run start : ", RunParams%RunDate%s(1:4), "-", RunParams%RunDate%s(5:6), &
         "-", RunParams%RunDate%s(7:8), " ", RunParams%RunTime%s(1:4)
      write (101, fmt="(a18,a)") "Input file name : ", RunParams%InputFile%s
      write (101, *)

      write (101, fmt="(a)") "Domain settings:"

      if (RunParams%isOneD) then
         write (101, fmt="(a)") "Spatial dimensions = 1"
      else
         write (101, fmt="(a)") "Spatial dimensions = 2"
      end if

      if (RunParams%Georeference) then

         write (101, fmt="(a,2P G0)") "Latitude of domain centre = ", RunParams%Lat
         write (101, fmt="(a,3P G0)") "Longitude of domain centre = ", RunParams%Lon

         distxy%first = -0.5_wp*RunParams%nXPoints*RunParams%deltaX + RunParams%centerUTM%first
         distxy%second = -0.5_wp*RunParams%nYPoints*RunParams%deltaY + RunParams%centerUTM%second
         pt = RunParams%projTransformer%utm_to_wgs84(distxy%first, distxy%second)
         write (101, fmt="(a,2P G0)") "Latitude of domain lower left = ", pt%first
         write (101, fmt="(a,3P G0)") "Longitude of domain lower left = ", LongitudeOut(pt%second)

         distxy%first = -0.5_wp*RunParams%nXPoints*RunParams%deltaX + RunParams%centerUTM%first
         distxy%second = 0.5_wp*RunParams%nYPoints*RunParams%deltaY + RunParams%centerUTM%second
         pt = RunParams%projTransformer%utm_to_wgs84(distxy%first, distxy%second)
         write (101, fmt="(a,2P G0)") "Latitude of domain upper left = ", pt%first
         write (101, fmt="(a,3P G0)") "Longitude of domain upper left = ", LongitudeOut(pt%second)

         distxy%first = 0.5_wp*RunParams%nXPoints*RunParams%deltaX + RunParams%centerUTM%first
         distxy%second = -0.5_wp*RunParams%nYPoints*RunParams%deltaY + RunParams%centerUTM%second
         pt = RunParams%projTransformer%utm_to_wgs84(distxy%first, distxy%second)
         write (101, fmt="(a,2P G0)") "Latitude of domain lower right = ", pt%first
         write (101, fmt="(a,3P G0)") "Longitude of domain lower right = ", LongitudeOut(pt%second)

         distxy%first = 0.5_wp*RunParams%nXPoints*RunParams%deltaX + RunParams%centerUTM%first
         distxy%second = 0.5_wp*RunParams%nYPoints*RunParams%deltaY + RunParams%centerUTM%second
         pt = RunParams%projTransformer%utm_to_wgs84(distxy%first, distxy%second)
         write (101, fmt="(a,2P G0)") "Latitude of domain upper right = ", pt%first
         write (101, fmt="(a,3P G0)") "Longitude of domain upper right = ", LongitudeOut(pt%second)

         write (101, fmt="(a,i2,a1)") "UTM zone = ", RunParams%UTM_zone_number, RunParams%UTM_zone_letter
         write (101, fmt="(a,i0)") "EPSG code = ", RunParams%utmEPSG
         write (101, fmt="(a,2P,G0)") "Central Easting = ", RunParams%centerUTM%first
         write (101, fmt="(a,2P,G0)") "Central Northing = ", RunParams%centerUTM%second
      end if

      write (101, fmt="(a)") "Boundary conditions = " // RunParams%bcs%s
      if (RunParams%bcs%s=='dirichlet') then
         write (101, fmt="(a,2P,G0)") "bcsHnval = ", RunParams%bcsHnval
         write (101, fmt="(a,2P,G0)") "bcsuval = ", RunParams%bcsuval
         write (101, fmt="(a,2P,G0)") "bcsvval = ", RunParams%bcsvval
         write (101, fmt="(a,2P,G0)") "bcspsival = ", RunParams%bcspsival
      end if

      write (101, fmt="(a,i0)") "nTiles = ", RunParams%nTiles
      write (101, fmt="(a,i0)") "nXtiles = ", RunParams%nXtiles
      write (101, fmt="(a,i0)") "nYtiles = ", RunParams%nYtiles
      write (101, fmt="(a,i0)") "nXpertile = ", RunParams%nXpertile
      write (101, fmt="(a,i0)") "nYpertile = ", RunParams%nYpertile
      write (101, fmt="(a,i0)") "nXPoints = ", RunParams%nXPoints
      write (101, fmt="(a,i0)") "nYPoints = ", RunParams%nYPoints
      write (101, fmt="(a,G0)") "Xtilesize = ", RunParams%Xtilesize
      write (101, fmt="(a,G0)") "Ytilesize = ", RunParams%Ytilesize
      write (101, fmt="(a,G0)") "xSize = ", RunParams%xSize
      write (101, fmt="(a,G0)") "ySize = ", RunParams%ySize
      write (101, fmt="(a,G0)") "deltaX = ", RunParams%deltaX
      write (101, fmt="(a,G0)") "deltaY = ", RunParams%deltaY
      write (101, *)

      write (101, fmt="(a)") "Initial Conditions:"
      if (RunParams%set_Caps) then
         write (101, fmt="(a,i0,a)") "Number of cap sources = ", RunParams%nCaps
         do J = 1, RunParams%nCaps
            write (101, fmt="(a,i0,a)") "Cap ", J, " :"
            write (101, fmt="(a,G0)") "x position = ", RunParams%CapSources(J)%x
            write (101, fmt="(a,G0)") "y position = ", RunParams%CapSources(J)%y
            write (101, fmt="(a,2P G0)") "latitude = ", RunParams%CapSources(J)%Lat
            write (101, fmt="(a,3P G0)") "longitude = ", RunParams%CapSources(J)%Lon
            write (101, fmt="(a,G0)") "radius = ", RunParams%CapSources(J)%Radius
            write (101, fmt="(a,G0)") "volume = ", RunParams%CapSources(J)%Volume
            write (101, fmt="(a,G0)") "height = ", RunParams%CapSources(J)%Height
            write (101, fmt="(a,G0)") "solids fraction = ", RunParams%CapSources(J)%psi
            write (101, fmt="(a,a4)") "shape = ", RunParams%CapSources(J)%Shape
         end do
      end if
      if (RunParams%set_Cubes) then
         write (101, fmt="(a,i0,a)") "Number of cube sources = ", RunParams%nCubes
         do J = 1, RunParams%nCubes
            write (101, fmt="(a,i0,a)") "Cube ", J, " :"
            write (101, fmt="(a,G0)") "x position = ", RunParams%CubeSources(J)%x
            write (101, fmt="(a,G0)") "y position = ", RunParams%CubeSources(J)%y
            write (101, fmt="(a,2P G0)") "latitude = ", RunParams%CubeSources(J)%Lat
            write (101, fmt="(a,3P G0)") "longitude = ", RunParams%CubeSources(J)%Lon
            write (101, fmt="(a,G0)") "length = ", RunParams%CubeSources(J)%Length
            write (101, fmt="(a,G0)") "width = ", RunParams%CubeSources(J)%Width
            write (101, fmt="(a,G0)") "volume = ", RunParams%CubeSources(J)%Volume
            write (101, fmt="(a,G0)") "height = ", RunParams%CubeSources(J)%Height
            write (101, fmt="(a,G0)") "solids fraction = ", RunParams%CubeSources(J)%psi
         end do
      end if
      if (RunParams%set_Sources) then
         write (101, fmt="(a,i0,a)") "Number of flux sources = ", RunParams%nSources
         do J = 1, RunParams%nSources
            write (101, fmt="(a,i0,a)") "Source ", J, " :"
            write (101, fmt="(a,G0)") "x position = ", RunParams%FluxSources(J)%x
            write (101, fmt="(a,G0)") "y position = ", RunParams%FluxSources(J)%y
            write (101, fmt="(a,2P G0)") "latitude = ", RunParams%FluxSources(J)%Lat
            write (101, fmt="(a,3P G0)") "longitude = ", RunParams%FluxSources(J)%Lon
            write (101, fmt="(a,G0)") "radius = ", RunParams%FluxSources(J)%Radius
            write (101, fmt="(a8)", advance="no") "time = ("
            do K = 1, RunParams%FluxSources(J)%nFluxSeries
               write (101, fmt="(G0)", advance="no") RunParams%FluxSources(J)%time(K)
               if (K < RunParams%FluxSources(J)%nFluxSeries) write (101, fmt="(a2)", advance="no") ", "
            end do
            write (101, fmt="(a1)") ")"
            write (101, fmt="(a8)", advance="no") "flux = ("
            do K = 1, RunParams%FluxSources(J)%nFluxSeries
               write (101, fmt="(G0)", advance="no") RunParams%FluxSources(J)%flux(K)
               if (K < RunParams%FluxSources(J)%nFluxSeries) write (101, fmt="(a2)", advance="no") ", "
            end do
            write (101, fmt="(a1)") ")"
            write (101, fmt="(a8)", advance="no") "conc = ("
            do K = 1, RunParams%FluxSources(J)%nFluxSeries
               write (101, fmt="(G0)", advance="no") RunParams%FluxSources(J)%psi(K)
               if (K < RunParams%FluxSources(J)%nFluxSeries) write (101, fmt="(a2)", advance="no") ", "
            end do
            write (101, fmt="(a1)") ")"
         end do
      end if
      write (101, *)

      write (101, fmt="(a)") "Parameters:"
      write (101, fmt="(a,G0)") 'g = ', RunParams%g
      write (101, fmt="(a,G0)") 'rhow = ', RunParams%rhow
      write (101, fmt="(a,G0)") 'rhos = ', RunParams%rhos
      write (101, fmt="(a,G0)") 'reduced g = ', RunParams%gred
      select case (RunParams%DragChoice%s)
      case ("Chezy")
         write (101, fmt="(a)") 'Drag : Chezy'
         write (101, fmt="(a,G0)") 'Chezy coefficient = ', RunParams%ChezyCo
      case ("Coulomb")
         write (101, fmt="(a)") 'Drag : Coulomb'
         write (101, fmt="(a,G0)") 'Coulomb coefficient = ', RunParams%CoulombCo
      case ("Voellmy")
         write (101, fmt="(a)") 'Drag : Voellmy'
         write (101, fmt="(a,G0)") 'Chezy coefficient = ', RunParams%ChezyCo
         write (101, fmt="(a,G0)") 'Coulomb coefficient = ', RunParams%CoulombCo
         write (101, fmt="(a,G0)") 'Voellmy switch rate = ', RunParams%VoellmySwitchRate
         write (101, fmt="(a,G0)") 'Voellmy switch value = ', RunParams%VoellmySwitchValue
      case ("Pouliquen")
         write (101, fmt="(a)") 'Drag : Pouliquen'
         write (101, fmt="(a,G0)") 'Pouliquen Min Slope = ', RunParams%PouliquenMinSlope
         write (101, fmt="(a,G0)") 'Pouliquen Max Slope = ', RunParams%PouliquenMaxSlope
      case ("Variable")
         write (101, fmt="(a)") 'Drag : Variable'
         write (101, fmt="(a,G0)") 'Chezy coefficient = ', RunParams%ChezyCo
         write (101, fmt="(a,G0)") 'Pouliquen Min Slope = ', RunParams%PouliquenMinSlope
         write (101, fmt="(a,G0)") 'Pouliquen Max Slope = ', RunParams%PouliquenMaxSlope
         write (101, fmt="(a,G0)") 'Voellmy switch rate = ', RunParams%VoellmySwitchRate
         write (101, fmt="(a,G0)") 'Voellmy switch value = ', RunParams%VoellmySwitchValue
      end select
      write (101, fmt="(a)") 'Erosion : '//RunParams%ErosionChoice%s
      if (RunParams%MorphodynamicsOn) then
         write (101, fmt="(a,G0)") 'Erosion rate = ', RunParams%EroRate
         write (101, fmt="(a,G0)") 'Granular erosion rate = ', RunParams%EroRateGranular
         write (101, fmt="(a,G0)") 'Critical shields = ', RunParams%CriticalShields
         write (101, fmt="(a,G0)") 'Erosion depth = ', RunParams%EroDepth
         write (101, fmt="(a,G0)") 'Erosion critical height = ', RunParams%EroCriticalHeight
         write (101, fmt="(a,G0)") 'Bed porosity = ', RunParams%BedPorosity
      end if
      write (101, fmt="(a,G0)") 'Maximum packing fraction = ', RunParams%maxPack
      write (101, fmt="(a,G0)") 'Solid diameter = ', RunParams%SolidDiameter
      write (101, fmt="(a,G0)") 'Particle Reynolds number = ', RunParams%Rep
      if (RunParams%MorphodynamicsOn) then
         write (101, fmt="(a,G0)") 'Deposition closure = ', RunParams%DepositionChoice%s
         write (101, fmt="(a,G0)") 'Particle settling speed in clear water = ', RunParams%ws0
         write (101, fmt="(a,G0)") 'Exponent in hindered settling term = ', RunParams%nsettling
         write (101, fmt="(a,G0)") 'Morphodynamic damping function = ', RunParams%MorphoDamp%s
      end if
      write (101, *)

      write (101, fmt="(a)") "Solver settings:"
      write (101, fmt="(a10,a7)") 'limiter : ', RunParams%limiter%s
      write (101, fmt="(a,G0)") 'height threshold = ', RunParams%heightThreshold
      write (101, fmt="(a,i0)") 'tile buffer = ', RunParams%TileBuffer
      write (101, fmt="(a,G0)") 'cfl number = ', RunParams%cfl
      write (101, fmt="(a,G0)") 'max time step = ', RunParams%maxdt
      write (101, fmt="(a,G0)") 'start time = ', RunParams%tstart
      write (101, fmt="(a,G0)") 'end time = ', RunParams%tend
      write (101, *)

      write (101, fmt="(a)") "Output settings:"
      write (101, fmt="(a,i0)") 'Number of output files = ', RunParams%Nout
      write (101, fmt="(a,G0)") 'Time step between outputs = ',  &
         (RunParams%tend - RunParams%tstart) / RunParams%Nout
      write (101, fmt="(a,i0)") 'Last output file = ', RunParams%CurrentOut

      write (101, *)
      write (101, fmt="(a)") "Topographic settings:"
      write (101, fmt="(a)") "Topography type = "//RunParams%Topog%s
      select case (RunParams%Topog%s)
         case ('srtm')
            write (101, fmt="(a)") "Topography path = "//RunParams%SRTMPath%s
            write (101, fmt="(a,a)") "SRTM virtual file = ", RunParams%OutDir%s//"SRTM.vrt"
         case ('raster', 'dem')
            write (101, fmt="(a)") "Topography path = "//RunParams%demPath%s
            write (101, fmt="(a)") "Raster file = "//RunParams%RasterFile%s
            if (RunParams%EmbedRaster) then
               write (101, fmt="(a)") "Embedded raster = on"
               write (101, fmt="(a,a)") "SRTM path = ", RunParams%SRTMPath%s
               write (101, fmt="(a,a)") "SRTM virtual file = ", RunParams%OutDir%s//"SRTM.vrt"
         else
            write (101, fmt="(a)") "Embedded raster = off"
         end if
        case default
            if (allocated(RunParams%TopogFuncParams)) then
               write (101, fmt="(a25)", advance="no") "Topography parameters = ("
               do J=1,size(RunParams%TopogFuncParams)
                  write (101, fmt="(G0,a2)", advance="no") RunParams%TopogFuncParams(J), ", "
               end do
               write (101, fmt="(a1)") ")"
            end if
      end select

      close (101)

   end subroutine OutputInfo

   
   ! Calculate volume and related quantities from solution data
   ! at a timestep and write to file.
   ! Inputs: RunParams - run settings
   !         t - current time
   !         grid - numerical solution data
   !         create [optional, default=.FALSE.] - should new file be created?
   subroutine CalculateVolume(RunParams, t, grid, create)
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: t
      type(GridType), target, intent(in) :: grid
      logical, intent(in), optional :: create

      type(varString) :: filename

      logical :: create_new

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      logical :: FileExists
      real(kind=wp) :: vol, mass, bed, solids
      real(kind=wp) :: c_vol, c_mass, c_bed, c_solids ! Compensated summation accumulators
      real(kind=wp) :: dv
      integer :: i, j
      integer :: nTiles, tt, ttk

      integer :: iHn, irho, iHnpsi, ibt

      real(kind=wp) :: rhow, rhos, rhob
      real(kind=wp) :: Hn, rho, gam

      filename = RunParams%OutDir + 'Volume.txt'

      create_new = .false.

      if (present(create)) create_new = create

      iHn = RunParams%Vars%Hn
      irho = RunParams%Vars%rho
      iHnpsi = RunParams%Vars%Hnpsi
      ibt = RunParams%Vars%bt

      rhow = RunParams%rhow
      rhos = RunParams%rhos
      rhob = rhow*RunParams%BedPorosity + rhos*(1.0_wp - RunParams%BedPorosity)

      nTiles = RunParams%nTiles

      inquire (file=filename%s, exist=FileExists)
      if ((.not.FileExists).or.(create_new)) then
         open (101, file=filename%s, status='replace')
         write(101,fmt="(a12,a2,6(a24,a2))") '        time', ', ', &
           '                  volume', ', ', '        total bed volume', ', ', &
           '              total mass', ', ', '                bed mass', ', ', &
           '       total solids mass', ', ', '         bed solids mass', ', '
      else
         open (101, file=filename%s, status="old", position="append")
      end if

      vol = 0.0_wp
      mass = 0.0_wp
      bed = 0.0_wp
      solids = 0.0_wp

      c_vol = 0.0_wp
      c_mass = 0.0_wp
      c_bed = 0.0_wp
      c_solids = 0.0_wp

      ! Pointers to simplify notation
      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)

         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile
               gam = GeometricCorrectionFactor(RunParams, tileContainer(ttk)%u(:, i, j))
               Hn = tileContainer(ttk)%u(iHn, i, j)
               rho = tileContainer(ttk)%u(irho, i, j)
               dv = (tileContainer(ttk)%u(RunParams%Vars%w, i, j) - tileContainer(ttk)%u(RunParams%Vars%b0, i, j) - &
                  tileContainer(ttk)%u(RunParams%Vars%bt, i, j))*gam*gam
               ! vol = vol + dv
               call KahanAdd(dv, vol, c_vol)
               ! mass = mass + rho*Hn*gam
               call KahanAdd(rho*Hn*gam, mass, c_mass)
               ! bed = bed + tileContainer(ttk)%u(ibt, i, j)
               call KahanAdd(tileContainer(ttk)%u(ibt, i, j), bed, c_bed)
               ! solids = solids + tileContainer(ttk)%u(iHnpsi, i, j)*gam
               call KahanAdd(tileContainer(ttk)%u(iHnpsi, i, j)*gam, solids, c_solids)
            end do
         end do
      end do

      if (RunParams%isOneD) then
         vol = vol*grid%deltaX
         mass = mass*grid%deltaX
         bed = bed*grid%deltaX
         solids = solids*grid%deltaX
      else
         vol = vol*grid%deltaX*grid%deltaY
         mass = mass*grid%deltaX*grid%deltaY
         bed = bed*grid%deltaX*grid%deltaY
         solids = solids*grid%deltaX*grid%deltaY
      end if

      ! Write out
      ! t, volume, total mass, total bed mass, total solid mass, total bed solid mass
      write(101,fmt="(f12.2,a2,6(es24.15,a2))") t, ', ', vol, ', ', bed, ', ', mass, ', ', bed * rhob, ', ', solids * rhos, ', ', bed * rhos * (1.0_wp - RunParams%BedPorosity), ', '

      close (101)

   end subroutine CalculateVolume

   !!!! The following routines are private

   pure function LongitudeOut(long_in) result(long_out)

      real(kind=wp), intent(in) :: long_in
      real(kind=wp) :: long_out

      if (long_in > 180_wp) then
         long_out = long_in - 360_wp
      else
         long_out = long_in
      end if

      return

   end function LongitudeOut

   subroutine OutputSolutionData_txt(RunParams, filename, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(varString), intent(in) :: filename
      type(GridType), target, intent(in) :: grid

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      logical :: FileExists
      !integer :: d
      integer :: i, j
      integer :: nTiles, tt, ttk, nYvertices
      type(varString) :: filename_full
      type(varString) :: topo_filename_full

      real(kind=wp) :: Hn, spd
      type(pair) :: latlon

      type(varString) :: compress_cmd

      nTiles = RunParams%nTiles

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      filename_full = CheckPath(RunParams%OutDir) + filename
      topo_filename_full = CheckPath(RunParams%OutDir) + filename + '_topo'

      inquire (file=filename_full%s, exist=FileExists)
      if (FileExists) then
         open (101, file=filename_full%s, status='replace')
      else
         open (101, file=filename_full%s, status='new')
      end if
      inquire (file=topo_filename_full%s, exist=FileExists)
      if (FileExists) then
         open (102, file=topo_filename_full%s, status='replace')
      else
         open (102, file=topo_filename_full%s, status='new')
      end if

      if (RunParams%isOneD) then
         write (101, fmt="(a8,a2,11(a18,a2),a18)", advance='yes') &
            '    tile', ', ', & ! column 1 -- tile number
            '                 x', ', ', & ! column 2 -- x coordinate from origin at domain centre
            '        flow_depth', ', ', & ! column 3 -- flow depth
            '        flow_speed', ', ', & ! column 4 -- flow speed
            '        x_velocity', ', ', & ! column 5 -- x velocity, u
            '           density', ', ', & ! column 6 -- density, rho
            '   solids_fraction', ', ', & ! column 7 -- solids fraction, psi
            '            x_flux', ', ', & ! column 8 -- x flux, rhoHnu
            '             Hnpsi', ', ', & ! column 9 -- solids volume per unit area, Hnpsi
            '                 w', ', ', & ! column 10 -- w
            '    base_elevation', ', ', & ! column 11 -- base topographic elevation, b0
            '  elevation_change', ', ', & ! column 12 -- change in topographic elevation, bt
            '           x_slope'          ! column 13 -- x-slope, db/dx
      else
         write (101, fmt="(a8,a2,18(a18,a2),a18)", advance='yes') &
            '    tile', ', ', &           ! column 1 -- tile number
            '                 x', ', ', & ! column 2 -- x coordinate from origin at domain centre
            '                 y', ', ', & ! column 3 -- y coordinate from origin at domain centre
            '          latitude', ', ', & ! column 4 -- latitude
            '         longitude', ', ', & ! column 5 -- longitude
            '        flow_depth', ', ', & ! column 6 -- flow depth
            '        flow_speed', ', ', & ! column 7 -- flow speed
            '        x_velocity', ', ', & ! column 8 -- x velocity, u
            '        y_velocity', ', ', & ! column 9 -- y velocity, v
            '           density', ', ', & ! column 10 -- density, rho
            '   solids_fraction', ', ', & ! column 11 -- solids fraction, psi
            '            x_flux', ', ', & ! column 12 -- x flux, rhoHnu
            '            y_flux', ', ', & ! column 13 -- y flux, rhoHnv
            '             Hnpsi', ', ', & ! column 14 -- solids volume per unit area, Hnpsi
            '                 w', ', ', & ! column 15 -- w
            '    base_elevation', ', ', & ! column 16 -- base topographic elevation, b0
            '  elevation_change', ', ', & ! column 17 -- change in topographic elevation, bt
            '           x_slope', ', ', & ! column 18 -- x-slope, db/dx
            '           y_slope'          ! column 19 -- y-slope
      end if
      write (102, fmt="(a8, a2, a18, a2, a18)", advance='yes') &
         '    tile', ', ', &           ! column 1 -- tile number
         '         B0_vertex', ', ', & ! column 2 -- base topographic elevation
         '         Bt_vertex'          ! column 3 -- change in topographic elevation

      if (.not. RunParams%Georeference) then
         latlon%first = 0.0_wp
         latlon%second = 0.0_wp
      end if

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile
               Hn = tileContainer(ttk)%u(RunParams%Vars%Hn, i, j)
               
               if (Hn > RunParams%heightThreshold) then
                  spd = sqrt(FlowSquaredSpeedSlopeAligned(Runparams, tileContainer(ttk)%u(:, i, j)))
               else
                  spd = 0.0_wp
               end if

               if (RunParams%isOneD) then
                  write (101, fmt="(i8, a2, 11(ES18.10E3, a2), ES18.10E3)", advance='yes') &
                     ttk, ', ', &
                     tileContainer(ttk)%x(i), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%Hn, i, j), ', ', &
                     spd, ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%u, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%rho, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%psi, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%rhoHnu, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%Hnpsi, i ,j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%w, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%b0, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%bt, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%dbdx, i, j)
               else
                  if (RunParams%Georeference) latlon = RunParams%projTransformer%utm_to_wgs84(tileContainer(ttk)%x(i), tileContainer(ttk)%y(j))
                  write (101, fmt="(i8, a2, 18(ES18.10E3, a2), ES18.10E3)", advance='yes') &
                     ttk, ', ', &
                     tileContainer(ttk)%x(i), ', ', &
                     tileContainer(ttk)%y(j), ', ', &
                     latlon%first, ', ', &
                     latlon%second, ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%Hn, i, j), ', ', &
                     spd, ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%u, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%v, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%rho, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%psi, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%rhoHnu, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%rhoHnv, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%Hnpsi, i ,j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%w, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%b0, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%bt, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%dbdx, i, j), ', ', &
                     tileContainer(ttk)%u(RunParams%Vars%dbdy, i, j)
               end if
            end do
            write (101, *)
         end do
      end do

      close (101)

      ! output vertex-centred topographic data
      nYvertices = 1
      if (.not. RunParams%isOneD) nYvertices = nYvertices + RunParams%nYpertile
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         do j = 1, nYvertices
            do i = 1, RunParams%nXpertile + 1
               write (102, fmt="(i8, a2, ES18.10E3, a2, ES18.10E3)", advance='yes') &
                  ttk, ', ', &
                  tileContainer(ttk)%b0(i, j), ', ', &
                  tileContainer(ttk)%bt(i, j)
            end do
            write (102, *)
         end do
      end do
      close (102)

      if (RunParams%compressOutput) then
         call InfoMessage("compressing... ")
         compress_cmd = varString('tar czf ')
         compress_cmd = compress_cmd + RunParams%OutDir + filename + '.tar.gz ' &
            + RunParams%OutDir + filename &
            + ' && rm ' + RunParams%OutDir + filename_full
         call InfoMessage(compress_cmd%s)
         call execute_command_line(compress_cmd%s)
      end if

   end subroutine OutputSolutionData_txt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if HAVE_NETCDF4
   subroutine OutputSolutionData_NetCDF(RunParams, filename, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      character(len=*), intent(in) :: filename
      type(GridType), target, intent(in) :: grid

      type(varString) :: filename_full

      integer :: tt, ttk, i, j

      integer :: ncid
      integer :: x_dim_id, x_vertex_dim_id
      integer :: y_dim_id, y_vertex_dim_id
      integer, dimension(:) :: dimids(2), vertex_dim_ids(2)
      integer :: x_id, y_id
      integer :: crs_id
      integer :: w_id, Hn_id, u_id, v_id, spd_id, psi_id
      integer :: rho_id, rhoHnu_id, rhoHnv_id, Hnpsi_id
      integer :: b0_id, bt_id, dbdx_id, dbdy_id

      integer :: x_vertex_id, y_vertex_id
      integer :: B0_vertex_id, Bt_vertex_id

      integer :: nc_status

      integer :: nXpertile, nYpertile, nTiles
      integer :: nXtiles, nYtiles
      integer :: nX, nX_vertex
      integer :: nY, nY_vertex

      integer :: tile_i, tile_j

      integer :: x_start, y_start
      integer :: x_vertex_start, y_vertex_start
      integer :: nX_vertex_pertile, nY_vertex_pertile
      integer, dimension(:) :: xy_start(2)
      integer, dimension(:) :: nXYpertile(2)
      integer, dimension(:) :: nXY_vertex_pertile(2)

      character(len=1) :: hemisphere

      integer :: tile_left, tile_bottom

      real(kind=wp), dimension(:, :), allocatable :: spd

      if (RunParams%Lat < 0) then
         hemisphere = 'S'
      else
         hemisphere = 'N'
      end if

      tile_left = maxval(grid%ActiveTiles%List)
      tile_bottom = tile_left

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile
      nTiles = grid%ActiveTiles%Size

      nXYpertile = [nXpertile, nYpertile]
      nX_vertex_pertile = nXpertile + 1
      nY_vertex_pertile = nYpertile + 1
      nXY_vertex_pertile = [nX_vertex_pertile, nY_vertex_pertile]

      nXtiles = int((grid%xmax - grid%xmin + RunParams%deltaX)/(RunParams%Xtilesize))
      nYtiles = int((grid%ymax - grid%ymin + RunParams%deltaY)/(RunParams%Ytilesize))

      nX = nXtiles*nXpertile
      nY = nYtiles*nYpertile

      nX_vertex = nXtiles*(nXpertile + 1)
      nY_vertex = nYtiles*(nYpertile + 1)

      do tt = 1, nTiles

         ttk = grid%ActiveTiles%List(tt)

         call GridCoords(ttk, grid, tile_i, tile_j)

         if (tile_i < tile_left) tile_left = tile_i
         if (tile_j < tile_bottom) tile_bottom = tile_j
      end do

      filename_full = RunParams%OutDir + trim(filename)

      nc_status = nf90_create(path=filename_full%s, cmode=NF90_NETCDF4, ncid=ncid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_create '//filename_full%s)

      call put_nc_att(ncid, 'title', 'Kestrel: '//trim(filename))
      call put_nc_att(ncid, 'institution', 'University of Bristol')
      call put_nc_att(ncid, 'source', 'Kestrel')

      call put_nc_att(ncid, 'Conventions', 'CF-1.11-draft')
      call put_nc_att(ncid, '_FillValue', -9999.9_wp)

      call netcdf_put_params(ncid,RunParams)

      call put_nc_att(ncid, 'tiles', grid%ActiveTiles%List)

      nc_status = nf90_def_dim(ncid, "x", nX, x_dim_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_dim x_dim_id')

      nc_status = nf90_def_dim(ncid, "y", nY, y_dim_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_dim y_dim_id')

      nc_status = nf90_def_dim(ncid, "x_vertex", nX_vertex, x_vertex_dim_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_dim x_vertex_dim_id')

      nc_status = nf90_def_dim(ncid, "y_vertex", nY_vertex, y_vertex_dim_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_dim y_vertex_dim_id')

      call define_nc_var_real(ncid, "x", [x_dim_id], x_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='x-coordinate', &
                              standard_name='projection_x_coordinate', &
                              units='m', &
                              axis='X')

      call define_nc_var_real(ncid, "y", [y_dim_id], y_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='y-coordinate', &
                              standard_name='projection_y_coordinate', &
                              units='m', &
                              axis='Y')

      call define_nc_var_real(ncid, "x_vertex", [x_vertex_dim_id], x_vertex_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='x-coordinate of cell vertices', &
                              standard_name='x_cell_vertices', &
                              units='m', &
                              axis='X')

      call define_nc_var_real(ncid, "y_vertex", [y_vertex_dim_id], y_vertex_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='y-coordinate of cell vertices', &
                              standard_name='y_cell_vertices', &
                              units='m', &
                              axis='Y')

      call set_nc_crs(ncid, RunParams%UTM_zone_number, hemisphere, RunParams%utmEPSG, crs_id, include_wkt=.FALSE.)

      dimids = [x_dim_id, y_dim_id]
      vertex_dim_ids = [x_vertex_dim_id, y_vertex_dim_id]

      allocate (spd(nXpertile, nYpertile))

      

      call define_nc_var_real(ncid, "flow_depth", dimids, Hn_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='flow depth', &
                              standard_name='flow_depth', &
                              units='m', &
                              coordinates='x y', &
                              grid_mapping='crs')
    
      call define_nc_var_real(ncid, "flow_speed", dimids, spd_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='flow speed', &
                              standard_name='flow_speed', &
                              units='m/s', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "x_velocity", dimids, u_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='flow velocity in the Easting direction', &
                              standard_name='x_velocity', &
                              units='m/s', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "y_velocity", dimids, v_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='flow velocity in the Northing direction', &
                              standard_name='y_velocity', &
                              units='m/s', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "density", dimids, rho_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='bulk density', &
                              standard_name='density', &
                              units='kg/m^3', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "solids_fraction", dimids, psi_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='solids fraction', &
                              standard_name='solids_fraction', &
                              units='1', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "x_flux", dimids, rhoHnu_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='momentum flux per unit area in the Easting direction', &
                              standard_name='x_flux', &
                              units='kg/m/s', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "y_flux", dimids, rhoHnv_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='momentum flux per unit area in the Northing direction', &
                              standard_name='y_flux', &
                              units='kg/m/s', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "Hnpsi", dimids, Hnpsi_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='volume of solids per unit area', &
                              standard_name='Hnpsi', &
                              units='m', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "w", dimids, w_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='w', &
                              standard_name='w', &
                              units='m', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "base_elevation", dimids, b0_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='initial topographic elevation', &
                              standard_name='base_elevation', &
                              units='m', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "elevation_change", dimids, bt_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='change in topographic elevation', &
                              standard_name='elevation_change', &
                              units='m', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "x_slope", dimids, dbdx_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='topographic slope in Easting direction', &
                              standard_name='x_slope', &
                              units='1', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "y_slope", dimids, dbdy_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='topographic slope in Northing direction', &
                              standard_name='y_slope', &
                              units='1', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "B0_vertex", vertex_dim_ids, B0_vertex_id, &
                              fill_value=-9999.9_wp, deflate=1, &
                              units='m', &
                              long_name="initial topographic elevation at cell vertices", &
                              standard_name="B0_vertex", &
                              coordinates='x_vertex y_vertex')

      call define_nc_var_real(ncid, "Bt_vertex", vertex_dim_ids, Bt_vertex_id, &
                              fill_value=-9999.9_wp, deflate=1, &
                              units='m', &
                              long_name="time varying topographic elevation at cell vertices", &
                              standard_name="Bt_vertex", &
                              coordinates='x_vertex y_vertex')

      nc_status = nf90_enddef(ncid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_enddef')

      if (RunParams%Georeference) then
         call put_nc_var(ncid, x_id, [(grid%xmin + RunParams%deltaX*i + RunParams%centerUTM%first, i=0, nX - 1)], start=[1], count=[nX])
         call put_nc_var(ncid, y_id, [(grid%ymin + RunParams%deltaY*i + RunParams%centerUTM%second, i=0, nY - 1)], start=[1], count=[nY])
      else
         call put_nc_var(ncid, x_id, [(grid%xmin + RunParams%deltaX*i, i=0, nX - 1)], start=[1], count=[nX])
         call put_nc_var(ncid, y_id, [(grid%ymin + RunParams%deltaY*i, i=0, nY - 1)], start=[1], count=[nY])
      end if

      do tt = 1, nTiles

         ttk = grid%ActiveTiles%List(tt)

         call GridCoords(ttk, grid, tile_i, tile_j)

         x_start = (tile_i - tile_left)*nXpertile + 1
         y_start = (tile_j - tile_bottom)*nYpertile + 1
         xy_start = [x_start, y_start]

         x_vertex_start = (tile_i - tile_left)*(nXpertile + 1) + 1
         y_vertex_start = (tile_j - tile_bottom)*(nYpertile + 1) + 1

         call put_nc_var(ncid, w_id, grid%tileContainer(ttk)%u(RunParams%Vars%w, :, :), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, Hn_id, grid%tileContainer(ttk)%u(RunParams%Vars%Hn, :, :), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, u_id,  grid%tileContainer(ttk)%u(RunParams%Vars%u, :, :),  start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, v_id,  grid%tileContainer(ttk)%u(RunParams%Vars%v, :, :),  start=xy_start, count=nXYpertile)

         do i = 1, RunParams%nXpertile
            do j = 1, RunParams%nYpertile
               spd(i, j) = sqrt(FlowSquaredSpeedSlopeAligned(Runparams, grid%tileContainer(ttk)%u(:, i, j)))
            end do
         end do

         call put_nc_var(ncid, spd_id, spd, start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, psi_id, grid%tileContainer(ttk)%u(RunParams%Vars%psi, :, :), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, rho_id, grid%tileContainer(ttk)%u(RunParams%Vars%rho, :, :), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, rhoHnu_id, grid%tileContainer(ttk)%u(RunParams%Vars%rhoHnu, :, :), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, rhoHnv_id, grid%tileContainer(ttk)%u(RunParams%Vars%rhoHnv, :, :), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, Hnpsi_id, grid%tileContainer(ttk)%u(RunParams%Vars%Hnpsi, : ,:), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, b0_id, grid%tileContainer(ttk)%u(RunParams%Vars%b0, :, :), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, bt_id, grid%tileContainer(ttk)%u(RunParams%Vars%bt, :, :), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, dbdx_id, grid%tileContainer(ttk)%u(RunParams%Vars%dbdx, :, :), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, dbdy_id, grid%tileContainer(ttk)%u(RunParams%Vars%dbdy, :, :), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, x_vertex_id, grid%tileContainer(ttk)%x_vertex(:), start=[x_vertex_start], count=[nX_vertex_pertile])
         call put_nc_var(ncid, y_vertex_id, grid%tileContainer(ttk)%y_vertex(:), start=[y_vertex_start], count=[nY_vertex_pertile])
         call put_nc_var(ncid, B0_vertex_id, grid%tileContainer(ttk)%B0(:, :), start=[x_vertex_start, y_vertex_start], count=nXY_vertex_pertile)
         call put_nc_var(ncid, Bt_vertex_id, grid%tileContainer(ttk)%Bt(:, :), start=[x_vertex_start, y_vertex_start], count=nXY_vertex_pertile)

      end do

      nc_status = nf90_close(ncid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_close')

      deallocate (spd)

   end subroutine OutputSolutionData_NetCDF
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine OutputMaximums_txt(RunParams, filename, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      character(len=*), intent(in) :: filename
      type(GridType), target, intent(in) :: grid

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      type(pair) :: latlon

      logical :: FileExists
      integer :: i, j
      integer :: nTiles, tt, ttk

      nTiles = RunParams%nTiles

      ! These pointers and variables simplify notation.
      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      inquire (file=filename, exist=FileExists)
      if (FileExists) then
         open (102, file=filename, status='replace')
      else
         open (102, file=filename, status='new')
      end if

      write (102, fmt="(a8,a2,14(a18,a2), a18)", advance='yes') &
         '    tile', ', ', &           ! column 1 -- tile number
         '        x_distance', ', ', & ! column 2 -- x coordinate from origin at domain centre
         '        y_distance', ', ', & ! column 3 -- y coordinate from origin at domain centre
         '          latitude', ', ', & ! column 4 -- latitude
         '         longitude', ', ', & ! column 5 -- longitude
         '         max_depth', ', ', & ! column 6 -- maximum flow depth
         '       t_max_depth', ', ', & ! column 7 -- time of maximum flow depth
         '         max_speed', ', ', & ! column 8 -- maximum flow speed
         '       t_max_speed', ', ', & ! column 9 -- time of maximum flow speed
         '       max_erosion', ', ', & ! column 10 -- maximum eroded depth
         '     t_max_erosion', ', ', & ! column 11 -- time of maximum eroded depth
         '       max_deposit', ', ', & ! column 12 -- maximum deposited depth
         '     t_max_deposit', ', ', & ! column 13 -- time of maximum deposited depth
         '   max_solids_frac', ', ', & ! column 14 -- maximum solids fraction
         ' t_max_solids_frac', ', ', & ! column 15 -- time of maximum maximum solids fraction
         '   inundation_time'          ! column 16 -- time of first inundation

      if (.not. RunParams%Georeference) then
         latlon%first = 0.0_wp
         latlon%second = 0.0_wp
      end if

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile
               if (RunParams%Georeference) then
                  latlon = RunParams%projTransformer%utm_to_wgs84(tileContainer(ttk)%x(i), tileContainer(ttk)%y(j))
               end if
               write (102, fmt="(i8,a2,14(ES18.10E3, a2), ES18.10E3)", advance='yes') &
                  ttk, ', ', &
                  tileContainer(ttk)%x(i), ', ', &
                  tileContainer(ttk)%y(j), ', ', &
                  latlon%first, ', ', &
                  latlon%second, ', ', &
                  tileContainer(ttk)%Hnmax(i, j, 1), ', ', &
                  tileContainer(ttk)%Hnmax(i, j, 2), ', ', &
                  tileContainer(ttk)%umax(i, j, 1), ', ', &
                  tileContainer(ttk)%umax(i, j, 2), ', ', &
                  tileContainer(ttk)%emax(i, j, 1), ', ', &
                  tileContainer(ttk)%emax(i, j, 2), ', ', &
                  tileContainer(ttk)%dmax(i, j, 1), ', ', &
                  tileContainer(ttk)%dmax(i, j, 2), ', ', &
                  tileContainer(ttk)%psimax(i, j, 1), ', ', &
                  tileContainer(ttk)%psimax(i, j, 2), ', ', &
                  tileContainer(ttk)%tfirst(i, j)
            end do
            write (102, *)
         end do
         write (102, *)
      end do

      close (102)

   end subroutine OutputMaximums_txt

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if HAVE_NETCDF4

   subroutine OutputMaximums_NetCDF(RunParams, filename, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      character(len=*), intent(in) :: filename
      type(GridType), target, intent(in) :: grid

      integer :: tt, ttk

      integer :: ncid
      integer :: x_dim_id
      integer :: y_dim_id

      integer :: x_id, y_id
      integer :: crs_id
      integer :: Hn_id, tHn_id
      integer :: spd_id, tspd_id
      integer :: e_id, te_id
      integer :: d_id, td_id
      integer :: psi_id, tpsi_id
      integer :: it_id
      integer :: nc_status

      integer :: nXpertile, nYpertile, nTiles
      integer :: nXtiles, nYtiles
      integer :: nX, nY

      integer :: tile_i, tile_j, i

      integer :: x_start, y_start
      integer, dimension(:) :: xy_start(2)
      integer, dimension(:) :: nXYpertile(2)

      character(len=1) :: hemisphere

      integer :: tile_left, tile_bottom

      if (RunParams%Lat < 0) then
         hemisphere = 'S'
      else
         hemisphere = 'N'
      end if

      tile_left = maxval(grid%ActiveTiles%List)
      tile_bottom = tile_left

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile
      nTiles = grid%ActiveTiles%Size

      nXYpertile = [nXpertile, nYpertile]

      nXtiles = int((grid%xmax - grid%xmin + RunParams%deltaX)/(RunParams%Xtilesize))
      nYtiles = int((grid%ymax - grid%ymin + RunParams%deltaY)/(RunParams%Ytilesize))

      nX = nXtiles*nXpertile
      nY = nYtiles*nYpertile

      do tt = 1, nTiles

         ttk = grid%ActiveTiles%List(tt)

         call GridCoords(ttk, grid, tile_i, tile_j)

         if (tile_i < tile_left) tile_left = tile_i
         if (tile_j < tile_bottom) tile_bottom = tile_j
      end do

      nc_status = nf90_create(path=filename, cmode=NF90_NETCDF4, ncid=ncid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_create '//filename)

      call put_nc_att(ncid, 'title', 'Kestrel: Maximums')
      call put_nc_att(ncid, 'institution', 'University of Bristol')
      call put_nc_att(ncid, 'source', 'Kestrel')

      call put_nc_att(ncid, 'Conventions', 'CF-1.11-draft')
      call put_nc_att(ncid, '_FillValue', -9999.9_wp)

      call netcdf_put_params(ncid,RunParams)

      call put_nc_att(ncid, 'tiles', grid%ActiveTiles%List)

      nc_status = nf90_def_dim(ncid, "x", nX, x_dim_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_dim x_dim_id')

      nc_status = nf90_def_dim(ncid, "y", nY, y_dim_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_dim y_dim_id')
      
      call define_nc_var_real(ncid, "x", [x_dim_id], x_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='x-coordinate', &
                              standard_name='projection_x_coordinate', &
                              units='m', &
                              axis='X')

      call define_nc_var_real(ncid, "y", [y_dim_id], y_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='y-coordinate', &
                              standard_name='projection_y_coordinate', &
                              units='m', &
                              axis='Y')

      call set_nc_crs(ncid, RunParams%UTM_zone_number, hemisphere, RunParams%utmEPSG, crs_id, include_wkt=.FALSE.)

      call define_nc_var_real(ncid, "max_depth", [x_dim_id, y_dim_id], Hn_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='maximum flow depth', &
                              units='m', &
                              standard_name='max_depth', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "t_max_depth", [x_dim_id, y_dim_id], tHn_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='time of maximum flow depth', &
                              units='s', &
                              standard_name='t_max_speed', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "max_speed", [x_dim_id, y_dim_id], spd_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='maximum flow speed', &
                              units='m/s', &
                              standard_name='max_speed', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "t_max_speed", [x_dim_id, y_dim_id], tspd_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='time of maximum flow speed', &
                              units='s', &
                              standard_name='t_max_speed', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "max_erosion", [x_dim_id, y_dim_id], e_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='maximum depth of erosion', &
                              units='m', &
                              standard_name='max_erosion', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "t_max_erosion", [x_dim_id, y_dim_id], te_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='time of maximum depth of erosion', &
                              units='s', &
                              standard_name='t_max_erosion', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "max_deposit", [x_dim_id, y_dim_id], d_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='maximum depth of deposit', &
                              units='m', &
                              standard_name='max_deposit', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "t_max_deposit", [x_dim_id, y_dim_id], td_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='time of maximum depth of deposit', &
                              units='s', &
                              standard_name='t_max_deposit', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "max_solids_frac", [x_dim_id, y_dim_id], psi_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='maximum solids fraction', &
                              units='1', &
                              standard_name='max_solids_frac', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "t_max_solids_frac", [x_dim_id, y_dim_id], tpsi_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='time of maximum solids fraction', &
                              units='s', &
                              standard_name='t_max_solids_frac', &
                              coordinates='x y', &
                              grid_mapping='crs')

      call define_nc_var_real(ncid, "inundation_time", [x_dim_id, y_dim_id], it_id, fill_value=-9999.9_wp, deflate=1, &
                              long_name='time of first inundation', &
                              units='s', &
                              standard_name='inundation_time', &
                              coordinates='x y', &
                              grid_mapping='crs')

      nc_status = nf90_enddef(ncid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_enddef')

      call put_nc_var(ncid, x_id, [(grid%xmin + RunParams%deltaX*i + RunParams%centerUTM%first, i=0, nX - 1)], start=[1], count=[nX])
      call put_nc_var(ncid, y_id, [(grid%ymin + RunParams%deltaY*i + RunParams%centerUTM%second, i=0, nY - 1)], start=[1], count=[nY])

      call put_nc_var(ncid, x_id, [(grid%xmin + RunParams%deltaX*i + RunParams%centerUTM%first, i=0, nX - 1)], start=[1], count=[nX])
      call put_nc_var(ncid, y_id, [(grid%ymin + RunParams%deltaY*i + RunParams%centerUTM%second, i=0, nY - 1)], start=[1], count=[nY])

      do tt = 1, nTiles

         ttk = grid%ActiveTiles%List(tt)

         call GridCoords(ttk, grid, tile_i, tile_j)

         x_start = (tile_i - tile_left)*nXpertile + 1
         y_start = (tile_j - tile_bottom)*nYpertile + 1

         xy_start = [x_start, y_start]

         call put_nc_var(ncid, Hn_id, grid%tileContainer(ttk)%Hnmax(:, :, 1), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, tHn_id, grid%tileContainer(ttk)%Hnmax(:, :, 2), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, spd_id, grid%tileContainer(ttk)%umax(:, :, 1), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, tspd_id, grid%tileContainer(ttk)%umax(:, :, 2), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, e_id, grid%tileContainer(ttk)%emax(:, :, 1), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, te_id, grid%tileContainer(ttk)%emax(:, :, 2), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, d_id, grid%tileContainer(ttk)%dmax(:, :, 1), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, td_id, grid%tileContainer(ttk)%dmax(:, :, 2), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, psi_id, grid%tileContainer(ttk)%psimax(:, :, 1), start=xy_start, count=nXYpertile)
         call put_nc_var(ncid, tpsi_id, grid%tileContainer(ttk)%psimax(:, :, 2), start=xy_start, count=nXYpertile)

         call put_nc_var(ncid, it_id, grid%tileContainer(ttk)%tfirst(:, :), start=xy_start, count=nXYpertile)

      end do

      nc_status = nf90_close(ncid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_close(ncid)')

   end subroutine OutputMaximums_NetCDF
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine OutputMaxHeights_KML(RunParams, filename, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      character(len=*), intent(in) :: filename
      type(GridType), target, intent(in) :: grid

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      logical :: FileExists
      integer :: i, j, k
      integer :: nTiles, tt, ttk

      type(pair), dimension(4) :: CellLatLon
      character(len=8) :: cstr

      nTiles = RunParams%nTiles

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      inquire (file=filename, exist=FileExists)
      if (FileExists) then
         open (101, file=filename, status='replace')
      else
         open (101, file=filename, status='new')
      end if

      write (101, fmt="(a38)") '<?xml version="1.0" encoding="UTF-8"?>'
      write (101, fmt="(a44)") '<kml xmlns="http://www.opengis.net/kml/2.2">'
      write (101, fmt="(2x,a10)") '<Document>'
      !write(101,fmt="(2x,a8)") '<Folder>'
      !write(101,fmt="(2x,a69)") '<ListStyle><listItemType>checkHideChildren</listItemType></ListStyle>'
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile

               if (tileContainer(ttk)%Hnmax(i, j, 1) > RunParams%kmlHeight) then
                  CellLatLon = LatLongCell(RunParams,tileContainer(ttk)%x(i),tileContainer(ttk)%y(j),RunParams%deltaX,RunParams%deltaY)

                  cstr = HEX_heights(tileContainer(ttk)%Hnmax(i, j, 1), 10.0_wp)

                  write (101, fmt="(4x,a11)") '<Placemark>'
                  !write(101,fmt="(6x,a34)") '<styleUrl>#transRedPoly</styleUrl>'
                  write (101, fmt="(6x,a29)") '<Style id="examplePolyStyle">'
                  write (101, fmt="(8x,a11)") '<PolyStyle>'
                  write (101, fmt="(10x,a7)") '<color>'
                  write (101, fmt="(12x,a8)") cstr
                  write (101, fmt="(10x,a8)") '</color>'
                  write (101, fmt="(10x,a29)") '<colorMode>normal</colorMode>'
                  write (101, fmt="(10x,a20)") '<outline>0</outline>'
                  write (101, fmt="(8x,a12)") '</PolyStyle>'
                  write (101, fmt="(6x,a8)") '</Style>'
                  write (101, fmt="(6x,a9)") '<Polygon>'
                  write (101, fmt="(8x,a20)") '<extrude>1</extrude>'
                  write (101, fmt="(8x,a45)") '<altitudeMode>relativeToGround</altitudeMode>'
                  write (101, fmt="(8x,a17)") '<outerBoundaryIs>'
                  write (101, fmt="(10x,a12)") '<LinearRing>'
                  write (101, fmt="(12x,a13)", advance='no') '<coordinates>'
                  do k = 1, 4
                     if (k == 1) then
                        write (101, fmt="(f0.16,a1,f0.16,a1,f0.16)", advance='yes') &
                           CellLatLon(k)%second, ',', &
                           CellLatLon(k)%first, ',', &
                           tileContainer(ttk)%Hnmax(i, j, 1)
                     else
                        write (101, fmt="(14x,f0.16,a1,f0.16,a1,f0.16)", advance='yes') &
                           CellLatLon(k)%second, ',', &
                           CellLatLon(k)%first, ',', &
                           tileContainer(ttk)%Hnmax(i, j, 1)
                     end if
                  end do
                  write (101, fmt="(14x,f0.16,a1,f0.16,a1,f0.16)", advance='no') &
                     CellLatLon(1)%second, ',', &
                     CellLatLon(1)%first, ',', &
                     tileContainer(ttk)%Hnmax(i, j, 1)
                  write (101, fmt="(a14)") '</coordinates>'
                  write (101, fmt="(10x,a13)") '</LinearRing>'
                  write (101, fmt="(8x,a18)") '</outerBoundaryIs>'
                  write (101, fmt="(6x,a10)") '</Polygon>'
                  write (101, fmt="(4x,a12)") '</Placemark>'
               end if
            end do
         end do
      end do
      !write(101,fmt="(2x,a9)") '</Folder>'
      write (101, fmt="(2x,a11)") '</Document>'
      write (101, fmt="(a6)") '</kml>'

      close (101)

      return

   end subroutine OutputMaxHeights_KML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine OutputInundationTime_KML(RunParams, filename, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      character(len=*), intent(in) :: filename
      type(GridType), target, intent(in) :: grid

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      logical :: FileExists
      integer :: i, j, k
      integer :: nTiles, tt, ttk

      type(pair), dimension(4) :: CellLatLon
      character(len=8) :: cstr

      nTiles = RunParams%nTiles

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      inquire (file=filename, exist=FileExists)
      if (FileExists) then
         open (101, file=filename, status='replace')
      else
         open (101, file=filename, status='new')
      end if

      write (101, fmt="(a38)") '<?xml version="1.0" encoding="UTF-8"?>'
      write (101, fmt="(a44)") '<kml xmlns="http://www.opengis.net/kml/2.2">'
      write (101, fmt="(2x,a10)") '<Document>'
      !write(101,fmt="(2x,a8)") '<Folder>'
      !write(101,fmt="(2x,a69)") '<ListStyle><listItemType>checkHideChildren</listItemType></ListStyle>'
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile

               if (tileContainer(ttk)%tfirst(i, j) .ne. -1.0_wp) then

                  CellLatLon = LatLongCell(RunParams,tileContainer(ttk)%x(i),tileContainer(ttk)%y(j),RunParams%deltaX,RunParams%deltaY)

                  cstr = HEX2KML(HEX_times(tileContainer(ttk)%tfirst(i, j), RunParams%Tend))

                  write (101, fmt="(4x,a11)") '<Placemark>'
                  !write(101,fmt="(6x,a34)") '<styleUrl>#transRedPoly</styleUrl>'
                  write (101, fmt="(6x,a29)") '<Style id="examplePolyStyle">'
                  write (101, fmt="(8x,a11)") '<PolyStyle>'
                  write (101, fmt="(10x,a7)") '<color>'
                  write (101, fmt="(12x,a8)") cstr
                  write (101, fmt="(10x,a8)") '</color>'
                  write (101, fmt="(10x,a29)") '<colorMode>normal</colorMode>'
                  write (101, fmt="(10x,a20)") '<outline>0</outline>'
                  write (101, fmt="(8x,a12)") '</PolyStyle>'
                  write (101, fmt="(6x,a8)") '</Style>'
                  write (101, fmt="(6x,a9)") '<Polygon>'
                  write (101, fmt="(8x,a20)") '<extrude>1</extrude>'
                  write (101, fmt="(8x,a45)") '<altitudeMode>relativeToGround</altitudeMode>'
                  write (101, fmt="(8x,a17)") '<outerBoundaryIs>'
                  write (101, fmt="(10x,a12)") '<LinearRing>'
                  write (101, fmt="(12x,a13)", advance='no') '<coordinates>'
                  do k = 1, 4
                     if (k == 1) then
                        write (101, fmt="(f0.16,a1,f0.16,a1,f0.16)", advance='yes') &
                           CellLatLon(k)%second, ',', &
                           CellLatLon(k)%first, ',', &
                           tileContainer(ttk)%Hnmax(i, j, 1)
                     else
                        write (101, fmt="(14x,f0.16,a1,f0.16,a1,f0.16)", advance='yes') &
                           CellLatLon(k)%second, ',', &
                           CellLatLon(k)%first, ',', &
                           tileContainer(ttk)%Hnmax(i, j, 1)
                     end if
                  end do
                  write (101, fmt="(14x,f0.16,a1,f0.16,a1,f0.16)", advance='no') &
                     CellLatLon(1)%second, ',', &
                     CellLatLon(1)%first, ',', &
                     tileContainer(ttk)%Hnmax(i, j, 1)
                  write (101, fmt="(a14)") '</coordinates>'
                  write (101, fmt="(10x,a13)") '</LinearRing>'
                  write (101, fmt="(8x,a18)") '</outerBoundaryIs>'
                  write (101, fmt="(6x,a10)") '</Polygon>'
                  write (101, fmt="(4x,a12)") '</Placemark>'
               end if
            end do
         end do
      end do
      !write(101,fmt="(2x,a9)") '</Folder>'
      write (101, fmt="(2x,a11)") '</Document>'
      write (101, fmt="(a6)") '</kml>'

      close (101)

      return

   end subroutine OutputInundationTime_KML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine OutputMaxSpeeds_KML(RunParams, filename, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      character(len=*), intent(in) :: filename
      type(GridType), target, intent(in) :: grid

      type(tileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      logical :: FileExists
      integer :: i, j, k
      integer :: nTiles, tt, ttk

      type(pair), dimension(4) :: CellLatLon
      character(len=8) :: cstr

      nTiles = RunParams%nTiles

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      inquire (file=filename, exist=FileExists)
      if (FileExists) then
         open (101, file=filename, status='replace')
      else
         open (101, file=filename, status='new')
      end if

      write (101, fmt="(a38)") '<?xml version="1.0" encoding="UTF-8"?>'
      write (101, fmt="(a44)") '<kml xmlns="http://www.opengis.net/kml/2.2">'
      write (101, fmt="(2x,a10)") '<Document>'
      !write(101,fmt="(2x,a8)") '<Folder>'
      !write(101,fmt="(2x,a69)") '<ListStyle><listItemType>checkHideChildren</listItemType></ListStyle>'
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile

               if (tileContainer(ttk)%Hnmax(i, j, 1) > RunParams%kmlHeight) then
                  CellLatLon = LatLongCell(RunParams,tileContainer(ttk)%x(i),tileContainer(ttk)%y(j),RunParams%deltaX,RunParams%deltaY)

                  cstr = HEX_heights(tileContainer(ttk)%umax(i, j, 1), 10.0_wp)

                  write (101, fmt="(4x,a11)") '<Placemark>'
                  !write(101,fmt="(6x,a34)") '<styleUrl>#transRedPoly</styleUrl>'
                  write (101, fmt="(6x,a29)") '<Style id="examplePolyStyle">'
                  write (101, fmt="(8x,a11)") '<PolyStyle>'
                  write (101, fmt="(10x,a7)") '<color>'
                  write (101, fmt="(12x,a8)") cstr
                  write (101, fmt="(10x,a8)") '</color>'
                  write (101, fmt="(10x,a29)") '<colorMode>normal</colorMode>'
                  write (101, fmt="(10x,a20)") '<outline>0</outline>'
                  write (101, fmt="(8x,a12)") '</PolyStyle>'
                  write (101, fmt="(6x,a8)") '</Style>'
                  write (101, fmt="(6x,a9)") '<Polygon>'
                  write (101, fmt="(8x,a20)") '<extrude>1</extrude>'
                  write (101, fmt="(8x,a45)") '<altitudeMode>relativeToGround</altitudeMode>'
                  write (101, fmt="(8x,a17)") '<outerBoundaryIs>'
                  write (101, fmt="(10x,a12)") '<LinearRing>'
                  write (101, fmt="(12x,a13)", advance='no') '<coordinates>'
                  do k = 1, 4
                     if (k == 1) then
                        write(101,fmt="(f0.16,a1,f0.16,a1,f0.16)",advance='yes') CellLatLon(k)%second, ',', CellLatLon(k)%first, ',', tileContainer(ttk)%umax(i,j,1)
                     else
                        write(101,fmt="(14x,f0.16,a1,f0.16,a1,f0.16)",advance='yes') CellLatLon(k)%second, ',', CellLatLon(k)%first, ',', tileContainer(ttk)%umax(i,j,1)
                     end if
                  end do
                  write(101,fmt="(14x,f0.16,a1,f0.16,a1,f0.16)",advance='no') CellLatLon(1)%second, ',', CellLatLon(1)%first, ',', tileContainer(ttk)%umax(i,j,1)
                  write (101, fmt="(a14)") '</coordinates>'
                  write (101, fmt="(10x,a13)") '</LinearRing>'
                  write (101, fmt="(8x,a18)") '</outerBoundaryIs>'
                  write (101, fmt="(6x,a10)") '</Polygon>'
                  write (101, fmt="(4x,a12)") '</Placemark>'
               end if
            end do
         end do
      end do
      !write(101,fmt="(2x,a9)") '</Folder>'
      write (101, fmt="(2x,a11)") '</Document>'
      write (101, fmt="(a6)") '</kml>'

      close (101)

      return

   end subroutine OutputMaxSpeeds_KML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function HEX_heights(h, hmax) result(cstr)

      implicit none

      real(kind=wp), intent(in) :: h
      real(kind=wp), intent(in) :: hmax
      character(len=8) :: cstr

      integer :: j

      do j = 1, 9
         if (h < j*hmax/9.0_wp) exit
      end do

      if (j == 1) then
         cstr = 'fffff5f0'
      elseif (j == 2) then
         cstr = 'fffee0d2'
      elseif (j == 3) then
         cstr = 'fffcbba1'
      elseif (j == 4) then
         cstr = 'fffc9272'
      elseif (j == 5) then
         cstr = 'fffb6a4a'
      elseif (j == 6) then
         cstr = 'ffef3b2c'
      elseif (j == 7) then
         cstr = 'ffcb181d'
      elseif (j == 8) then
         cstr = 'ffa50f15'
      else
         cstr = 'ff67000d'
      end if

      return

   end function HEX_heights

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function HEX_times(t, tmax) result(cstr)

      implicit none

      real(kind=wp), intent(in) :: t
      real(kind=wp), intent(in) :: tmax
      character(len=8) :: cstr

      integer :: j

      do j = 1, 9
         if (t < j*tmax/9.0_wp) exit
      end do

      if (j == 1) then
         cstr = '440154'
      elseif (j == 2) then
         cstr = '482878'
      elseif (j == 3) then
         cstr = '3e4989'
      elseif (j == 4) then
         cstr = '31688e'
      elseif (j == 5) then
         cstr = '26828e'
      elseif (j == 6) then
         cstr = '1f9e89'
      elseif (j == 7) then
         cstr = '35b779'
      elseif (j == 8) then
         cstr = '6ece58'
      else
         cstr = 'b5de2b'
      end if

      return

   end function HEX_times

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function HEX2KML(hexcode, alpha) result(kmlcode)
      ! Converts the format from RRGGBB to AABBGGRR.
      implicit none

      character(len=6), intent(in) :: hexcode
      character(len=2), intent(in), optional :: alpha
      character(len=8) :: kmlcode

      character(len=2) :: a

      if (present(alpha)) then
         a = alpha
      else
         a = "ff"
      end if

      kmlcode = a//hexcode(5:6)//hexcode(3:4)//hexcode(1:2)

      return

   end function HEX2KML
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if HAVE_NETCDF4

   subroutine netcdf_put_params(ncid,RunParams)

      integer, intent(in) :: ncid
      type(RunSet), intent(in) :: RunParams

      call put_nc_att(ncid, "time", RunParams%tstart + RunParams%CurrentOut * RunParams%DeltaT)

      call put_nc_att(ncid, "nXpertile", RunParams%nXpertile)
      call put_nc_att(ncid, "nYpertile", RunParams%nYpertile)

      call put_nc_att(ncid, "nXtiles", RunParams%nXtiles)
      call put_nc_att(ncid, "nYtiles", RunParams%nYtiles)
      call put_nc_att(ncid, "nXPoints", RunParams%nXPoints)
      call put_nc_att(ncid, "nYPoints", RunParams%nYPoints)
      call put_nc_att(ncid, "Xtilesize", RunParams%Xtilesize)
      call put_nc_att(ncid, "Ytilesize", RunParams%Ytilesize)
      call put_nc_att(ncid, "xSize", RunParams%xSize)
      call put_nc_att(ncid, "ySize", RunParams%ySize)
      call put_nc_att(ncid, "deltaX", RunParams%deltaX)
      call put_nc_att(ncid, "deltaY", RunParams%deltaY)
      call put_nc_att(ncid, "centre_latitude", RunParams%Lat)
      call put_nc_att(ncid, "centre_longitude", RunParams%Lon)
      call put_nc_att(ncid, "limiter", RunParams%limiter)
      call put_nc_att(ncid, "heightThreshold", RunParams%heightThreshold)
      if (RunParams%SpongeLayer) then
         call put_nc_att(ncid, "SpongeLayer", 'On')
      else
         call put_nc_att(ncid, "SpongeLayer", 'Off')
      end if
      call put_nc_att(ncid, "SpongeStrength", RunParams%SpongeStrength)
      call put_nc_att(ncid, "TileBuffer", RunParams%TileBuffer)
      call put_nc_att(ncid, "cfl", RunParams%cfl)
      call put_nc_att(ncid, "maxdt", RunParams%maxdt)
      call put_nc_att(ncid, "tstart", RunParams%tstart)
      call put_nc_att(ncid, "tend", RunParams%tend)
      call put_nc_att(ncid, "DeltaT", RunParams%DeltaT)
      if (RunParams%Restart) then
         call put_nc_att(ncid, "Restart", 'On')
      else
         call put_nc_att(ncid, "Restart", 'Off')
      end if
      call put_nc_att(ncid, "InitialCondition", RunParams%InitialCondition)
      call put_nc_att(ncid, "nSubSteps", RunParams%nSubSteps)
      if (RunParams%Georeference) then
         call put_nc_att(ncid, "central_easting", RunParams%centerUTM%first)
         call put_nc_att(ncid, "central_northing", RunParams%centerUTM%second)
         call put_nc_att(ncid, "crs_epsg", 'epsg:'//Int2String(RunParams%utmEPSG))
      end if
      if (RunParams%geometric_factors) then
         call put_nc_att(ncid, "Geometric factors", 'On')
      else
         call put_nc_att(ncid, "Geometric factors", 'Off')
      end if
      call put_nc_att(ncid, "g", RunParams%g)
      call put_nc_att(ncid, "water density", RunParams%rhow) ! Density of water
      call put_nc_att(ncid, "solids density", RunParams%rhos) ! Density of solids
      call put_nc_att(ncid, "reduced gravity", RunParams%gred) ! reduced gravity g' = (rhos-rhow)*g/rhow
      call put_nc_att(ncid, "bed porosity", RunParams%BedPorosity)
      call put_nc_att(ncid, "maximum packing", RunParams%maxPack)
      call put_nc_att(ncid, "solid diameter", RunParams%SolidDiameter)
      call put_nc_att(ncid, "eddy viscosity", RunParams%EddyViscosity)
      call put_nc_att(ncid, "particle Reynolds number", RunParams%Rep)
      call put_nc_att(ncid, "drag choice", RunParams%DragChoice%s)
      select case (RunParams%DragChoice%s)
         case ("Chezy")
            call put_nc_att(ncid, "Chezy coefficient", RunParams%ChezyCo)
         case ("Coulomb")
            call put_nc_att(ncid, "Coulomb coefficient", RunParams%CoulombCo)
         case ("Voellmy")
            call put_nc_att(ncid, "Chezy coefficient", RunParams%ChezyCo)
            call put_nc_att(ncid, "Coulomb coefficient", RunParams%CoulombCo)
         case ("Pouliquen")
            call put_nc_att(ncid, "Pouliquen minimum slope", RunParams%PouliquenMinSlope)
            call put_nc_att(ncid, "Pouliquen maximum slope", RunParams%PouliquenMaxSlope)
            call put_nc_att(ncid, "Pouliquen beta", RunParams%PouliquenBeta)
         case ("Variable")
            call put_nc_att(ncid, "Chezy coefficient", RunParams%ChezyCo)
            call put_nc_att(ncid, "Pouliquen minimum slope", RunParams%PouliquenMinSlope)
            call put_nc_att(ncid, "Pouliquen maximum slope", RunParams%PouliquenMaxSlope)
            call put_nc_att(ncid, "Pouliquen beta", RunParams%PouliquenBeta)
            call put_nc_att(ncid, "Voellmy switch function", RunParams%fswitch%s)
            call put_nc_att(ncid, "Voellmy switch rate", RunParams%VoellmySwitchRate)
            call put_nc_att(ncid, "Voellmy switch value", RunParams%VoellmySwitchValue)
         case ("Manning")
            call put_nc_att(ncid, "Manning coefficient", RunParams%ManningCo)
      end select
      if (RunParams%MorphodynamicsOn) then
         call put_nc_att(ncid, "morphodynamics time stepping", "on")
         call put_nc_att(ncid, "erosion choice", RunParams%ErosionChoice%s)
         select case (RunParams%ErosionChoice%s)
            case ('Fluid')
               call put_nc_att(ncid, "fluid erosion rate", RunParams%EroRate)
               call put_nc_att(ncid, "critical Shields number", RunParams%CriticalShields)
            case ('Granular')
               call put_nc_att(ncid, "granular erosion rate", RunParams%EroRateGranular)
            case ('Mixed')
               call put_nc_att(ncid, "fluid erosion rate", RunParams%EroRate)
               call put_nc_att(ncid, "granular erosion rate", RunParams%EroRateGranular)
               call put_nc_att(ncid, "critical Shields number", RunParams%CriticalShields)
         end select
         select case (RunParams%ErosionChoice%s)
            case ('Fluid', 'Granular', 'Mixed')
               call put_nc_att(ncid, "erosion depth", RunParams%EroDepth)
               call put_nc_att(ncid, "erosion critical height", RunParams%EroCriticalHeight)
               call put_nc_att(ncid, "erosion transition function", RunParams%ErosionTransition)
         end select
         call put_nc_att(ncid, "deposition closure", RunParams%DepositionChoice%s)
         call put_nc_att(ncid, "clear water settling speed", RunParams%ws0)
         call put_nc_att(ncid, "hindered settling exponent", RunParams%nsettling)
         call put_nc_att(ncid, "morphodynamic damping function", RunParams%MorphoDamp)
      else
         call put_nc_att(ncid, "morphodynamics time stepping", "off")
      end if
      
   end subroutine netcdf_put_params
#endif

   function LatLongCell(RunParams,Cellx,Celly,dx,dy) result(EdgeLatLon)
      !!! Gets the lat and long of the edges of a cell (with size dx*dy) whose centre point Cellx, Celly east/north of clatlong (approximate due to orthographic projection)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: Cellx, Celly
      real(kind=wp), intent(in) :: dx, dy
      type(pair), dimension(4) :: EdgeLatLon

      type(pair) :: Ptxy

      Ptxy%first = Cellx-0.5_wp*dx + RunParams%centerUTM%first
      Ptxy%second = Celly-0.5_wp*dy + RunParams%centerUTM%second
      EdgeLatLon(1) = RunParams%projTransformer%utm_to_wgs84(Ptxy%first, Ptxy%second)
      if (EdgeLatLon(1)%second>180_wp) EdgeLatLon(1)%second = EdgeLatLon(1)%second-360.0_wp

      Ptxy%first = Cellx+0.5_wp*dx + RunParams%centerUTM%first
      Ptxy%second = Celly-0.5_wp*dy + RunParams%centerUTM%second
      EdgeLatLon(2) = RunParams%projTransformer%utm_to_wgs84(Ptxy%first, Ptxy%second)
      if (EdgeLatLon(2)%second>180_wp) EdgeLatLon(2)%second = EdgeLatLon(2)%second-360.0_wp

      Ptxy%first = Cellx+0.5_wp*dx + RunParams%centerUTM%first
      Ptxy%second = Celly+0.5_wp*dy + RunParams%centerUTM%second
      EdgeLatLon(3) = RunParams%projTransformer%utm_to_wgs84(Ptxy%first, Ptxy%second)
      if (EdgeLatLon(3)%second>180_wp) EdgeLatLon(3)%second = EdgeLatLon(3)%second-360.0_wp

      Ptxy%first = Cellx-0.5_wp*dx + RunParams%centerUTM%first
      Ptxy%second = Celly+0.5_wp*dy + RunParams%centerUTM%second
      EdgeLatLon(4) = RunParams%projTransformer%utm_to_wgs84(Ptxy%first, Ptxy%second)
      if (EdgeLatLon(4)%second>180_wp) EdgeLatLon(4)%second = EdgeLatLon(4)%second-360.0_wp

      return

   end function LatLongCell

end module output_module
