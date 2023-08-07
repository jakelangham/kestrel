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


! The main purpose of this module is to define a structure (RunSet) for storing
! parameters and settings for the run, most of which are user-settable and read
! in at runtime. One RunSet structure is loaded when the program initialises and
! this is passed around to most routines.
module runsettings_module

   use, intrinsic :: iso_c_binding
   use set_precision_module, only: wp
   use varstring_module, only: varString
   use utilities_module, only: pair
   use utm_module, only: proj_transformer

   implicit none

   private
   public :: RunSet
   public :: CapStr, CubeStr, SourceStr
   public :: RasterData, FortranRasterData
   public :: SetImmutableRunSettings

   ! First we define a bunch of useful types that RunSet needs to know about.

   ! Indicies for the fields of the solution vector.
   type VarIndex
      integer :: w, rhoHnu, rhoHnv, Hnpsi ! Conserved quantities
      integer :: Hn, hp, u, v, psi, rho ! Derived variables
      integer :: b0, bt, dbdx, dbdy ! Specified elevation and slopes
   end type VarIndex

   ! Cap initial condition type. Defines an initial release of a volume of
   ! material with circular base.
   type Caps
      real(kind=wp) :: x, y ! x and y location of centre point.
      real(kind=wp) :: Lat, Lon ! latitude and longitude of centre point
      real(kind=wp) :: radius, volume, height, psi
      real(kind=wp) :: u, v ! Cap x and y velocities
      ! Shape sets the free surface and can be flat (const flow depth),
      ! parabolic, or 'level' (const elevation).
      character(len=6) :: shape
   end type Caps

   ! Associated strings for caps used in initialisation.
   type CapStr
      type(varString), dimension(:), allocatable :: x, y
      type(varString), dimension(:), allocatable :: radius
      type(varString), dimension(:), allocatable :: volume
      type(varString), dimension(:), allocatable :: height
      type(varString), dimension(:), allocatable :: u, v
      type(varString), dimension(:), allocatable :: psi
      type(varString), dimension(:), allocatable :: lat, lon
      type(varString), dimension(:), allocatable :: shape
   end type CapStr

   ! Cube initial condition type. Defines an initial release of material with a
   ! rectangular base.
   type Cubes
      real(kind=wp) :: x, y
      real(kind=wp) :: lat, lon
      real(kind=wp) :: length, width, volume, height, psi
      real(kind=wp) :: u, v
      ! Shape sets the free surface and can be flat (const flow depth) or
      ! 'level' (const elevation).
      character(len=6) :: shape 
   end type Cubes

   ! Associated strings for cubes used in initialisation.
   type CubeStr
      type(varString), dimension(:), allocatable :: x, y
      type(varString), dimension(:), allocatable :: width
      type(varString), dimension(:), allocatable :: length
      type(varString), dimension(:), allocatable :: height
      type(varString), dimension(:), allocatable :: u, v
      type(varString), dimension(:), allocatable :: psi
      type(varString), dimension(:), allocatable :: lat, lon
      type(varString), dimension(:), allocatable :: shape
   end type CubeStr

   ! Type for 'sources' of material, defined as fluxes across an area with given
   ! centre (x,y) and radius.
   type Sources
      real(kind=wp) :: x, y
      real(kind=wp) :: radius
      real(kind=wp) :: lat, lon
      integer :: nFluxSeries ! Number of elements in flux time series.
      ! Times series data: flux, time, and solids fraction.
      real(kind=wp), dimension(:), allocatable :: flux, time, psi
      integer :: NumCellsInSrc ! Number of cells within the source region.
   end type Sources

   ! Associated strings for sources used in initialisation.
   type SourceStr
      type(varString), dimension(:), allocatable :: x, y
      type(varString), dimension(:), allocatable :: radius
      type(varString), dimension(:), allocatable :: lat, lon
      type(varString), dimension(:), allocatable :: flux, time, psi
   end type SourceStr

   type, bind(C) :: RasterData
      integer(c_int) :: read_success
      integer(c_int32_t) :: x_size
      integer(c_int32_t) :: y_size
      integer(c_int32_t) :: size
      real(c_double) :: pixel_width
      real(c_double) :: pixel_height
      real(c_double) :: origin_x
      real(c_double) :: origin_y
      real(c_double) :: row_rotation
      real(c_double) :: column_rotation
      integer(c_int) :: EPSG_code
      logical(c_bool) :: latlon
      type(c_ptr) :: NW, NE, SE, SW
      type(c_ptr) :: values
      real(c_double) :: no_data
   end type

  ! RasterData structure (fortran version). This is linked with c_f_pointer to
  ! the version above.
   type :: FortranRasterData
      type(varString) :: RasterName
      integer(c_int) :: ReadSuccess
      integer(c_int32_t) :: RasterXSize
      integer(c_int32_t) :: RasterYSize
      integer(c_int32_t) :: RasterSize
      real(kind=wp) :: pixelWidth
      real(kind=wp) :: pixelHeight
      real(kind=wp) :: originX
      real(kind=wp) :: originY
      real(kind=wp) :: rowRotation
      real(kind=wp) :: columnRotation
      integer(c_int) :: EPSG_code
      logical(c_bool) :: latlon
      type(pair) :: NW, NE, SE, SW
      real(kind=wp), dimension(:), pointer :: values
      real(kind=wp) :: nodata
   end type

   ! Main structure. RunSet contains essentially all settings and parameters for
   ! the simulation. It is subdivided into sections which each have their own
   ! module for initialisation in the files *Settings.f90 and Parameters.f90. 
   type RunSet

   ! -- Domain settings. --

      ! Number grid points per tile in x and y.
      integer :: nXpertile, nYpertile
      ! Number of tiles in the domain in x and y.
      integer :: nXtiles, nYtiles
      ! nTiles = nXtiles * nYtiles.
      integer :: nTiles
      ! Total number of cells in the domain in x and y, i.e.
      ! nXPoints = nXpertile * nXtiles.
      integer :: nXPoints, nYPoints
      ! Total number of variables.
      integer :: nDimensions 
      ! Number of 'flux' variables, i.e. primary fields.
      integer :: nFlux 
      ! Dimensions of tiles in x and y.
      real(kind=wp) :: Xtilesize, Ytilesize
      ! Size of domain in x and y.
      real(kind=wp) :: xSize, ySize
      ! Size of finite volume cells in x and y and their reciprocals.
      real(kind=wp) :: deltaX, deltaY, deltaXRecip, deltaYRecip
      logical :: isOneD ! True for 1D simulations.
      ! Latitude / longitude of domain centre and associated georeferencing
      ! data.
      real(kind=wp) :: lat, lon
      type(pair) :: LatLon
      real(kind=wp) :: easting, northing
      type(pair) :: centerUTM
      type(proj_transformer) :: projTransformer
      integer :: UTM_zone_number
      character(len=1) :: UTM_zone_letter
      integer :: utmEPSG
      ! Boundary conditions at domain edges + variables for proscribed values,
      ! if needed.
      type(varString) :: bcs
      real(kind=wp) :: bcsHnval, bcsuval, bcsvval, bcspsival

   ! -- Initial Conditions. --
      logical :: set_Caps
      integer :: nCaps
      type(Caps), dimension(:), allocatable :: CapSources

      logical :: set_Cubes
      integer :: nCubes
      type(Cubes), dimension(:), allocatable :: CubeSources

      logical :: set_Sources
      integer :: nSources
      type(Sources), dimension(:), allocatable :: FluxSources

   ! -- Parameters. --
      ! Use geometrically corrected coordinates?
      logical :: geometric_factors
      real(kind=wp) :: g ! gravity
      real(kind=wp) :: rhow ! Density of water
      real(kind=wp) :: rhos ! Density of solids
      real(kind=wp) :: gred ! reduced gravity g' = (rhos-rhow)*g/rhow
      ! Drag parameters.
      real(kind=wp) :: ChezyCo ! Chezy coefficient
      real(kind=wp) :: ManningCo ! Manning coefficient
      real(kind=wp) :: CoulombCo ! Coulomb coefficient
      real(kind=wp) :: PouliquenMinSlope ! Minimum slope angle for steady flow with Pouliquen friction
      real(kind=wp) :: PouliquenMaxSlope ! Maximum slope angle for steady flow with Pouliquen friction
      real(kind=wp) :: PouliquenBeta ! beta parameter in Pouliquen friction
      real(kind=wp) :: VoellmySwitchRate ! Steepness of transition in Voellmy drag
      real(kind=wp) :: VoellmySwitchValue ! Centre value of transition in Voellmy drag
      ! Morphodynamic parameters.
      logical :: MorphodynamicsOn ! Whether or not to time step the morphodynamics.
      real(kind=wp) :: EroRate ! Erosion rate (dimensional mass flux)
      real(kind=wp) :: EroRateGranular ! Granular erosion rate (dimensional mass flux)
      real(kind=wp) :: CriticalShields ! Threshold on Shield number above which erosion occurs
      real(kind=wp) :: EroDepth ! Maximum depth of material that can be eroded
      real(kind=wp) :: EroCriticalHeight ! Flow depth that must be exceeded for erosion to occur
      real(kind=wp) :: BedPorosity ! Porosity of the bed
      real(kind=wp) :: maxPack ! maximum packing volume fraction
      real(kind=wp) :: SolidDiameter ! Diameter of solids
      real(kind=wp) :: Rep ! Particle Reynolds number
      real(kind=wp) :: ws0 ! Particle settling speed in clear water
      real(kind=wp) :: nsettling ! Exponent in hindered settling term
      ! Sets a simple turbulent diffusivity \nabla (\nu \rho h \nabla u)
      real(kind=wp) :: EddyViscosity 

      ! Closure functions.

      ! Choice of switching function. Can be: tanh [default], rat3, cos, linear, equal, zero/0/off, one/1
      type(varString) :: fswitch 
      ! Choice of function for erosion transition. Can be: smooth [default], step, off
      type(varString) :: ErosionTransition 
      ! Choice of drag function. Can be: Chezy, Coulomb, Voellmy, Pouliquen, Variable [default], Manning
      type(varString) :: DragChoice 
      ! Choice of deposition function.  Can be: None, Simple, Spearman Manning [default]
      type(varString) :: DepositionChoice
      ! Choice of erosion function. Can be: Fluid, Granular, Mixed, Off
      type(varString) :: ErosionChoice 
      ! Choice of morphodynamic damping function.  Can be: tanh [default], rat3, none/off
      type(varString) :: MorphoDamp 

   ! -- Solver settings. --
      ! Limiting function for spatial derivatives. 
      type(varString) :: limiter
      ! Small height below which desingularisation for derived variables (e.g.
      ! u & v) is applied in the numerical scheme. Also used elsewhere when we
      ! occasionally need to define cutoffs based on flow depth.
      real(kind=wp) :: heightThreshold
      ! Sets a buffer around flowing area where tiles are activated.
      integer :: TileBuffer
      real(kind=wp) :: cfl ! Courant-Friedrichs-Lewy time stepping constant
      real(kind=wp) :: diffusiveTimeScale ! min{deltaX^2 / nu, deltaY^2 / nu}
      real(kind=wp) :: maxdt ! User set maximum time step
      real(kind=wp) :: tstart ! In seconds, start simulation at t = tstart
      real(kind=wp) :: tend   ! In seconds, end simulation at t = tend
      logical :: Restart ! Is this a restart of a previous run?
      type(varString) :: InitialCondition ! File containing initial condition
      integer :: nSubSteps ! number of substeps for multi-step method
      ! Settings for 'sponge layer' to damp out flow near domain edge.
      logical :: SpongeLayer
      real(kind=wp) :: SpongeStrength

   ! -- Output settings. --
      ! Various output file names.
      type(varString) :: basePath
      type(varString) :: OutDir
      type(varString) :: TopogFilename
      type(varString) :: InfoFilename
      type(varString) :: MaxHeightFilename
      type(varString) :: MaxSpeedFilename
      type(varString) :: MaxErosionFilename
      type(varString) :: MaxDepositFilename
      type(varString) :: InundationTimeFilename
      type(varString) :: MaximumsFilename
      ! Depth threshold above which data should beFilename included in KML
      ! output.
      real(kind=wp) :: kmlHeight 
      ! Which output types to use?: txt, NetCDF, KML.
      logical :: out_txt = .false.
      logical :: out_nc = .false.
      logical :: out_kml = .false.
      logical :: CompressOutput ! In the case of .txt files, compress with tar.

      integer :: Nout ! Number of output files
      real(kind=wp) :: DeltaT ! Time interval between outputs
      ! These variables track the iteration over the output intervals
      integer :: FirstOut = 0 ! Gets overwritten if restaring a simulation
      integer :: CurrentOut

   ! -- Topography settings. --
      type(varString) :: Topog
      type(varString) :: demPath
      type(varString) :: RasterFile
      type(varString) :: SRTMPath
      character(len=11), dimension(:), allocatable :: SRTMfiles
      real(kind=wp), dimension(:), allocatable :: TopogFuncParams
      type(FortranRasterData) :: SRTMtiles
      logical(kind=c_bool) :: EmbedRaster
      logical :: RebuildDEM
      logical :: Georeference ! Georeference only if using a DEM or SRTM

   ! -- Misc other data. --
      ! Variable indices.
      type(VarIndex) :: Vars
      integer, dimension(:), allocatable :: iFlux
      integer, dimension(:), allocatable :: iDesing
      ! Date / time of run
      type(varString) :: RunDate, RunTime
      ! Name of input file.
      type(varString) :: InputFile
      ! Indexes which equations are stepped semi-implicitly.
      logical, dimension(:), allocatable :: ImplicitStep

   end type RunSet

contains

   ! Initialise some mandatory run settings. Must always call this before
   ! starting a run.
   subroutine SetImmutableRunSettings(RunParams)
      implicit none

      type(RunSet), intent(inout) :: RunParams

      ! Indexing arrays that specify the primary and derived flow variables.
      allocate (RunParams%iFlux(4))
      RunParams%iFlux(1) = RunParams%Vars%w
      RunParams%iFlux(2) = RunParams%Vars%rhoHnu
      RunParams%iFlux(3) = RunParams%Vars%rhoHnv
      RunParams%iFlux(4) = RunParams%Vars%Hnpsi
      allocate (RunParams%iDesing(5))
      RunParams%iDesing(1) = RunParams%Vars%Hn
      RunParams%iDesing(2) = RunParams%Vars%u
      RunParams%iDesing(3) = RunParams%Vars%v
      RunParams%iDesing(4) = RunParams%Vars%psi
      RunParams%iDesing(5) = RunParams%Vars%rho

      RunParams%nFlux = size(RunParams%iFlux)

      ! Time stepper settings. Perhaps one day different solvers could be set
      ! from input file, for now just use default.

      ! Records which variables are to be time stepped implicitly.
      allocate (RunParams%ImplicitStep(RunParams%nDimensions))
      RunParams%ImplicitStep(:) = .false.
      RunParams%ImplicitStep(RunParams%Vars%rhoHnu) = .true.
      RunParams%ImplicitStep(RunParams%Vars%rhoHnv) = .true.

      RunParams%nSubsteps = 4

   end subroutine SetImmutableRunSettings

end module runsettings_module
