! This file is part of the Kestrel software for simulations
! of sediment-laden Earth surface flows.
!
! Version v1.1.1
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


! This modules handles checkpointing for the simulations. Flows can either be
! integrated forward from a given (user-set) initial condition, or there is the
! option to restart an existing simulation - useful if the simulation aborted
! for some reason. Data for restarts is recorded in RunParams%InfoFilename
! (default: RunInfo.txt).
module restart_module

   use iso_c_binding
   use set_precision_module, only: wp
   use runsettings_module, only: RunSet
   use utilities_module, only: CheckFileExists
   use messages_module, only: InfoMessage, WarningMessage, FatalErrorMessage
   use varstring_module, only: ReadFileLine, varString
   use output_module, only: GenerateOutputFilename
   use closures_module, only : GeometricCorrectionFactor
   use grid_module, only: GridCoords, GridType, TileList, TileType
   use update_tiles_module, only: ActivateTile, AddTile, AddToActiveTiles, AllocateTile
   use morphodynamic_rhs_module, only: ComputeInterfacialTopographicData
#if HAVE_NETCDF4
   use netcdf
   use netcdf_utils_module
#endif

   implicit none

   private
   public :: LoadInitialCondition

contains

   ! Load the initial condition for the simulation. Our preference is that if
   ! the user specified an initial condition in the input file, we load that,
   ! otherwise we check if the user wants to restart and try to do that.
   subroutine LoadInitialCondition(RunParams, grid)
      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(GridType), target, intent(inout) :: grid

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: activeTiles

      type(varString) :: InitFile, TopoInitFile
      type(varString) :: ext

      character(len=4096) :: cwd_str
      type(varString) :: cwd
      type(varString) :: decompress_cmd

      integer :: nXpertile, nYpertile, tt, ttk

      ! If restarting, need to read in settings from the InfoFile - in case of
      ! conflict we choose the InfoFile settings.
      if (RunParams%Restart .and. RunParams%InitialCondition%len() == 0) then
         call ReadRunInfoToRunParams(RunParams)
      end if
      
      grid%t = RunParams%tstart + RunParams%DeltaT * RunParams%FirstOut

      tileContainer => grid%tileContainer
      activeTiles => grid%activeTiles

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      if (RunParams%nSources >= 1) then
         RunParams%FluxSources(:)%NumCellsInSrc = 0
      end if

      ! Get the initial condition file name.
      if (RunParams%InitialCondition%len()>0) then
         InitFile = RunParams%InitialCondition
      else if (RunParams%Restart) then
         call GetLastResultFile(RunParams, InitFile, ext='')
         InitFile = RunParams%basePath + InitFile
      else
         return ! No initial conditions to load
      end if

      ! Load initial condition file, priority order: netCDF, txt, tar.gz
      if (CheckFileExists(InitFile + '.nc')) then
#if HAVE_NETCDF4
         TopoInitFile = InitFile + '.topo.nc'
         InitFile = InitFile + '.nc'
#else
         call FatalErrorMessage('NetCDF not active')
#endif
      else if (CheckFileExists(InitFile + '.txt')) then
         InitFile = InitFile + '.txt'
         TopoInitFile = InitFile + '_topo'
      else if (CheckFileExists(InitFile + 'txt.tar.gz')) then
         InitFile = InitFile + '.txt'
         TopoInitFile = InitFile + '_topo'
         call InfoMessage("uncompressing " // InitFile%s // '.tar.gz')
         cwd = varString()
         call getcwd(cwd_str)
         cwd = varString(cwd_str, trim_str=.TRUE.)
         decompress_cmd = varString('cd ') + RunParams%out_path &
                                   + ' && tar xzf ' + InitFile &
                                   + '.tar.gz ' + ' && cd ' + cwd%s
         call execute_command_line(decompress_cmd%s)
      else
         call FatalErrorMessage('Unable to locate initial condition file ' // &
                                InitFile%s // ' with valid extension.')
      end if

      ext = InitFile%get_extension()
      select case (ext%s)
      case ('txt')
         call LoadInitialCondition_txt(RunParams, grid, InitFile%s, TopoInitFile%s)
#if HAVE_NETCDF4
      case ('nc')
         call LoadInitialCondition_nc(RunParams, grid, InitFile%s, &
                                      RunParams%out_path%s // "Maximums.nc")
#endif
      end select

      ! This data is needed for time step computation etc, but isn't stored 
      ! explicitly in save states so we need to compute it here.
      do tt = 1, grid%ActiveTiles%Size
         ttk = grid%ActiveTiles%List(tt)
         call ComputeInterfacialTopographicData(RunParams, grid, grid%tileContainer, ttk)
      end do

   end subroutine LoadInitialCondition

   ! Load initial condition from .txt data in the files resultfile (containing
   ! the cell-centred solution data) and topofile (containing the topographic
   ! data at cell vertices).
   subroutine LoadInitialCondition_txt(RunParams, grid, resultfile, topofile)
      implicit none

      type(RunSet), intent(inout) :: RunParams
      character(len=*), intent(in) :: resultfile
      character(len=*), intent(in) :: topofile
      type(GridType), target, intent(inout) :: grid

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: activeTiles

      type(varString) :: MaximumsFilename
      character(len=2) :: com

      logical :: FileExists
      integer :: iostatus ! Status of input/output

      real(kind=wp) :: x, y, lat, lon, Hn, spd, u, v, rho, psi

      integer :: i, j, n
      integer :: kk, ttk
      integer :: nTiles, nYvertices

      real(kind=wp) :: srcX, srcY, srcR
      real(kind=wp) :: R2, b0

      open (52, file=resultfile, status='OLD', action='read') ! Open input file

      inquire (file=topofile, exist=FileExists)
      if (.not. FileExists) then
         call FatalErrorMessage('Unable to locate topography file ' // topofile // '.')
      else
         call InfoMessage('Reading file ' // topofile)
      end if
      open (153, file=topofile, status='OLD', action='read') ! Open input file

      tileContainer => grid%tileContainer
      activeTiles => grid%activeTiles

      ! Go through result file, scan list of tiles and add all of them before
      ! doing anything else.
      read (52, *) ! Move past header line
      nTiles = 0
      do
         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile
               read (52, fmt="(i8)", iostat=iostatus) ttk
               if (iostatus /= 0) goto 10
            end do
            if (.not. RunParams%isOneD .and. j /= RunParams%nYpertile) read (52, *)
         end do
         call AddToActiveTiles(grid, ttk)
         call AllocateTile(RunParams, grid, ttk)
         nTiles = nTiles + 1
         read (52, *)
      end do
10    rewind (52)

      ! Load in the data from resultfile.
      read (52, *) ! header line
      do n = 1, nTiles
         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile
               if (RunParams%isOneD) then
                  read (52, fmt="(i8, a2, 11(ES18.10E3, a2), ES18.10E3)") ttk, com, &
                     tileContainer(ttk)%x(i), com, &
                     Hn, com, spd, com, u, com, rho, com, psi, com, &
                     tileContainer(ttk)%u(RunParams%Vars%rhoHnu, i, j), com, &
                     tileContainer(ttk)%u(RunParams%Vars%Hnpsi, i, j), com, &
                     tileContainer(ttk)%u(RunParams%Vars%w, i, j), com, &
                     b0, com, &
                     tileContainer(ttk)%u(RunParams%Vars%bt, i, j), com, &
                     tileContainer(ttk)%u(RunParams%Vars%dbdx, i, j)
               else
                  read (52, fmt="(i8, a2, 18(ES18.10E3, a2), ES18.10E3)") ttk, com, &
                     tileContainer(ttk)%x(i), com, &
                     tileContainer(ttk)%y(j), com, &
                     lat, com, &
                     lon, com, &
                     Hn, com, spd, com, u, com, v, com, rho, com, psi, com, &
                     tileContainer(ttk)%u(RunParams%Vars%rhoHnu, i, j), com, &
                     tileContainer(ttk)%u(RunParams%Vars%rhoHnv, i, j), com, &
                     tileContainer(ttk)%u(RunParams%Vars%Hnpsi, i, j), com, &
                     tileContainer(ttk)%u(RunParams%Vars%w, i, j), com, &
                     b0, com, &
                     tileContainer(ttk)%u(RunParams%Vars%bt, i, j), com, &
                     tileContainer(ttk)%u(RunParams%Vars%dbdx, i, j), com, &
                     tileContainer(ttk)%u(RunParams%Vars%dbdy, i, j)
               end if

               ! Check if point lies in a flux source and if it is, add it
               if (RunParams%nSources >= 1) then ! flux sources
                  do kk = 1, RunParams%nSources
                     srcX = RunParams%FluxSources(kk)%x
                     srcR = RunParams%FluxSources(kk)%Radius

                     x = tileContainer(ttk)%x(i)

                     if (RunParams%isOneD) then
                        R2 = (x - srcX)*(x - srcX)
                     else
                        srcY = RunParams%FluxSources(kk)%y
                        y = tileContainer(ttk)%y(j)
                        R2 = (x - srcX)*(x - srcX) + (y - srcY)*(y - srcY)
                     end if
                     if (R2 <= srcR*srcR) then
                        RunParams%FluxSources(kk)%NumCellsInSrc = &
                           RunParams%FluxSources(kk)%NumCellsInSrc + 1
                     end if
                  end do
               end if
            end do ! i = 1, nXpertile
            if (.not. RunParams%isOneD .and. j /= RunParams%nYpertile) read (52, *)
         end do
         read (52, *)
      end do ! j = 1, nYpertile
      close (52)

      ! Load vertex-centred topographic data
      nYvertices = 1
      if (.not. RunParams%isOneD) nYvertices = nYvertices + RunParams%nYpertile
      read (153, *) ! header line
      do n = 1, nTiles
         do j = 1, nYvertices
            do i = 1, RunParams%nXpertile + 1
               read (153, fmt="(i8, a2, ES18.10E3, a2, ES18.10E3)") ttk, com, &
                  b0, com, &
                  tileContainer(ttk)%bt(i, j)
            end do ! i = 1, nXpertile
            if (.not. RunParams%isOneD .and. j /= RunParams%nYpertile + 1) read (153, *)
         end do
         read (153, *)
      end do ! j = 1, nYpertile
      close (153)

      do n = 1, nTiles
         ttk = grid%ActiveTiles%List(n)
         call ActivateTile(RunParams, grid, ttk, restart=.true.)
      end do

      ! Read in aggregate data. (Only need to do this if restarting.)
      if (RunParams%Restart) then
         MaximumsFilename = RunParams%out_path + RunParams%MaximumsFilename + '.txt'
         
         inquire (file=MaximumsFilename%s, exist=FileExists)
         if (.not. FileExists) then
            call FatalErrorMessage('Unable to locate '//MaximumsFilename%s//'.')
         else
            open (53, file=MaximumsFilename%s, status='OLD', action='read')
         end if

         ! Load max heights
         read (53, *) ! Skip header line
         do n = 1, nTiles
            do j = 1, RunParams%nYpertile
               do i = 1, RunParams%nXpertile
                  read (53, fmt="(i8, a2, 14(ES18.10E3, a2), ES18.10E3)") ttk, com, &
                     x, com, y, com, lat, com, lon, com, &
                     tileContainer(ttk)%Hnmax(i, j, 1), com, &
                     tileContainer(ttk)%Hnmax(i, j, 2), com, &
                     tileContainer(ttk)%umax(i, j, 1), com, &
                     tileContainer(ttk)%umax(i, j, 2), com, &
                     tileContainer(ttk)%emax(i, j, 1), com, &
                     tileContainer(ttk)%emax(i, j, 2), com, &
                     tileContainer(ttk)%dmax(i, j, 1), com, &
                     tileContainer(ttk)%dmax(i, j, 2), com, &
                     tileContainer(ttk)%psimax(i, j, 1), com, &
                     tileContainer(ttk)%psimax(i, j, 2), com, &
                     tileContainer(ttk)%tfirst(i, j)
               end do
               read (53, *)
            end do
            read (53, *)
         end do
         close (53)
      end if ! End 'if RunParams%Restart' block

   end subroutine LoadInitialCondition_txt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if HAVE_NETCDF4
   ! Load initial condition from .nc (NetCDF) data in the files resultfile
   ! (containing all the solution data) and maxfile (containing flow maximums,
   ! which are relevant if restarting).
   subroutine LoadInitialCondition_nc(RunParams, grid, resultfile, maxfile)
      implicit none

      type(RunSet), intent(inout) :: RunParams
      character(len=*), intent(in) :: resultfile
      character(len=*), intent(in) :: maxfile
      type(GridType), target, intent(inout) :: grid

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: activeTiles

      integer :: n
      integer :: ttk
      integer :: ncid, maxncid, nc_status

      integer :: nXpertile, nYpertile, nTiles
      integer :: nXtiles, nYtiles
      integer :: nX
      integer :: nY

      integer :: tile_i, tile_j

      integer :: x_start, y_start
      integer :: x_vertex_start, y_vertex_start
      integer, dimension(:) :: xy_start(2)
      integer, dimension(:) :: nXYpertile(2)
      integer, dimension(:) :: nXY_vertex_pertile(2)

      integer :: tile_left, tile_bottom

      real(kind=wp) :: central_easting, central_northing

      integer, dimension(:), allocatable :: tiles

      real(kind=wp), dimension(:), allocatable :: x, y

      character(len=:), allocatable :: ICversion

      tileContainer => grid%tileContainer
      activeTiles => grid%activeTiles

      ncid = nc_open_read(path=resultfile)

      ! Only get max data if restarting
      maxncid = 0
      if (RunParams%Restart) maxncid = nc_open_read(path=maxfile)

      tiles = nc_get_tiles(ncid)
      nTiles = size(tiles)

      ! Get tile sizes
      call get_nc_att(ncid, 'nXpertile', nXpertile)
      call get_nc_att(ncid, 'nYpertile', nYpertile)
      call get_nc_att(ncid, 'nXtiles', nXtiles)
      call get_nc_att(ncid, 'nYtiles', nYtiles)
      call get_nc_att(ncid, 'version', ICversion)

      if (ICversion /= RunParams%version) call WarningMessage("Results file " // resultfile // " made with different Kestrel version to current installation")

      if (RunParams%Georeference) then
        call get_nc_att(ncid, 'central_easting', central_easting)
        call get_nc_att(ncid, 'central_northing', central_northing)
      else
         central_easting = 0.0_wp
         central_northing = 0.0_wp
      end if

      nX = nXtiles*nXpertile
      nY = nYtiles*nYpertile

      nXYpertile = [nXpertile, nYpertile]
      nXY_vertex_pertile = [nXpertile + 1, nYpertile + 1]

      tile_left = maxval(tiles)
      tile_bottom = tile_left

      allocate (x(nXpertile), y(nYpertile))

      ! Get leftmost and bottommost tile coordinates for determining extent of
      ! simulation data.
      do n = 1, nTiles
         ttk = tiles(n)
         call GridCoords(ttk, grid, tile_i, tile_j)

         if (tile_i < tile_left) tile_left = tile_i
         if (tile_j < tile_bottom) tile_bottom = tile_j
      end do

      ! Read all the solution data in from the NetCDF, tile per tile.
      do n = 1, nTiles
         ttk = tiles(n)
         call AddTile(grid,ttk,RunParams)

         call GridCoords(ttk, grid, tile_i, tile_j)

         x_start = (tile_i - tile_left)*nXpertile + 1
         y_start = (tile_j - tile_bottom)*nYpertile + 1

         xy_start = [x_start, y_start]

         x_vertex_start = (tile_i - tile_left)*(nXpertile + 1) + 1
         y_vertex_start = (tile_j - tile_bottom)*(nYpertile + 1) + 1

         call get_nc_var(ncid, 'x', start=x_start, count=nXpertile, vals=x(:))
         call get_nc_var(ncid, 'y', start=y_start, count=nYpertile, vals=y(:))

         tileContainer(ttk)%x = x - central_easting
         tileContainer(ttk)%y = y - central_northing

         call get_nc_var(ncid, 'w', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%w, :, :))

         call get_nc_var(ncid, 'x_flux', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%rhoHnu, :, :))
         call get_nc_var(ncid, 'y_flux', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%rhoHnv, :, :))

         call get_nc_var(ncid, 'x_velocity', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%u, :, :))
         call get_nc_var(ncid, 'y_velocity', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%v, :, :))
         call get_nc_var(ncid, 'flow_depth', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%Hn, :, :))

         call get_nc_var(ncid, 'Hnpsi', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%Hnpsi, :, :))
         call get_nc_var(ncid, 'solids_fraction', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%psi, :, :))
         call get_nc_var(ncid, 'density', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%rho, :, :))

         call get_nc_var(ncid, 'base_elevation', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%b0, :, :))
         call get_nc_var(ncid, 'elevation_change', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%bt, :, :))
         call get_nc_var(ncid, 'x_slope', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%dbdx, :, :))
         call get_nc_var(ncid, 'y_slope', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%u(RunParams%Vars%dbdy, :, :))

         call get_nc_var(ncid, 'B0_vertex', start=[x_vertex_start, y_vertex_start], count=nXY_vertex_pertile, vals=tileContainer(ttk)%B0(:, :))
         call get_nc_var(ncid, 'Bt_vertex', start=[x_vertex_start, y_vertex_start], count=nXY_vertex_pertile, vals=tileContainer(ttk)%bt(:, :))

         ! Only get max data if restarting
         if (RunParams%Restart) then
            call get_nc_var(maxncid, 'max_depth', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%Hnmax(:, :, 1))
            call get_nc_var(maxncid, 't_max_depth', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%Hnmax(:, :, 2))
            call get_nc_var(maxncid, 'max_speed', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%umax(:, :, 1))
            call get_nc_var(maxncid, 't_max_speed', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%umax(:, :, 2))
            call get_nc_var(maxncid, 'max_erosion', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%emax(:, :, 1))
            call get_nc_var(maxncid, 't_max_erosion', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%emax(:, :, 2))
            call get_nc_var(maxncid, 'max_deposit', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%dmax(:, :, 1))
            call get_nc_var(maxncid, 't_max_deposit', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%dmax(:, :, 2))
            call get_nc_var(maxncid, 'max_solids_frac', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%psimax(:, :, 1))
            call get_nc_var(maxncid, 't_max_solids_frac', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%psimax(:, :, 2))
            call get_nc_var(maxncid, 'inundation_time', start=xy_start, count=nXYpertile, vals=tileContainer(ttk)%tfirst(:, :))
         end if

      end do

      call get_nc_att(ncid, 'DeltaT', RunParams%DeltaT)

      nc_status = nf90_close(ncid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'close nc')

      if (RunParams%Restart) then
         nc_status = nf90_close(maxncid)
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'close maxnc')
      end if

   end subroutine LoadInitialCondition_nc
#endif

   ! When restarting, this routine looks for the last file that was output by
   ! the simulation - i.e. the one to restart from. It does this by consulting
   ! RunParams%InfoFilename.
   subroutine GetLastResultFile(RunParams, LastFileName, ext)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(varString), intent(out) :: LastFileName
      character(len=*), optional, intent(in) :: ext

      type(varString) :: ext_in

      type(varString) :: InfoFile
      logical :: FileExists
      integer :: iostatus ! Status of input/output

      type(varString) :: line
      type(varString) :: label, val
      integer :: LastFile

      ! If there is a user-specified initial condition, then we should be
      ! running from that.
      if (RunParams%InitialCondition%s /= " ") then
         LastFileName = RunParams%InitialCondition
         return
      end if

      if (present(ext)) then
         ext_in%s = ext
      else
         ext_in%s = '.txt'
      end if

      InfoFile = RunParams%out_path + RunParams%InfoFilename

      inquire (file = InfoFile%s, exist = FileExists)

      if (.not. FileExists) then
         call FatalErrorMessage('The requested ' // RunParams%InfoFilename%s // &
                      ' file does not exist in ' // RunParams%basePath%s // '.')
      end if
      open (51, file = InfoFile%s) ! Open input file

      do ! Run through input file
        call ReadFileLine(51, line, iostat=iostatus) ! Read line of input file
         if (iostatus /= 0) exit ! If iostat returns read error, exit loop
         line = line%adjustl() ! Ignore leading blanks
         if (line%len() > 0) then ! if the line isn't blank, process it
            if (line%first_char() /= '%') then ! ignore lines starting with %
               if (line%contains("Last output file = ")) then
                  call line%split('=', label, remain=val)
                  LastFile = val%to_int()
               end if
            end if
         end if
      end do
      close (51)

      LastFileName = GenerateOutputFilename(RunParams%OutDir, LastFile, 6, ext=ext_in%s)

   end subroutine GetLastResultFile

   ! When restarting, this routine is called to read settings from
   ! RunParams%InfoFilename into the RunParams structure, to setup all the
   ! relevant settings from the simulation we wish to restart from.
   subroutine ReadRunInfoToRunParams(RunParams)
      implicit none

      type(RunSet), intent(inout) :: RunParams

      type(varString) :: InfoFile
      logical :: FileExists
      integer :: iostatus ! Status of input/output

      type(varString) :: line, label, val

      integer :: LastFile
      real(kind=wp) :: FileTimeStep, tstart

      type(varString) :: InfoVersion ! Get version string reported in InfoFile

      logical :: found_filetimestep
      logical :: found_lastfile
      logical :: found_tstart

      InfoFile = RunParams%basePath + RunParams%OutDir + RunParams%InfoFilename
      inquire (File=InfoFile%s, Exist=FileExists)

      if (.not. FileExists) then
         call FatalErrorMessage('To restart this simulation, a ' // &
            RunParams%InfoFilename%s // ' file should be present in ' // &
            RunParams%OutDir%s // '.')
      end if
      open (51, file=InfoFile%s) ! Open input file

      found_filetimestep = .false.
      found_lastfile = .false.
      found_tstart = .false.

      FileTimeStep = 0.0_wp
      LastFile = 0
      tstart = 0.0_wp

      do ! Run through info file
         call ReadFileLine(51, line, iostat=iostatus) ! Read line of input file
         if (iostatus /= 0) exit ! If iostat returns read error, exit loop
         line = line%adjustl() ! Ignore leading blanks
         if (line%len() > 0) then ! if the line isn't blank, process it
            if (line%contains("=")) then ! ignore that are not in keyword = value format
               call line%split("=", label, remain=val)
               label = label%to_lower()
               if (label%contains("Kestrel version")) then
                  InfoVersion = val
               else
                  InfoVersion = varString("Unknown")
               end if
               if (label%contains("time step between outputs")) then
                  FileTimeStep = val%to_real()
                  found_filetimestep = .true.
               end if
               if (label%contains("last output file")) then
                  LastFile = val%to_int()
                  found_lastfile = .true.
               end if
               if (label%contains("t start")) then
                  tstart = val%to_real()
                  found_tstart = .true.
               end if
            end if
         end if
      end do
      close (51)

      if (.not. found_filetimestep) call FatalErrorMessage("Could not find 'Time step between outputs' in RunInfo.txt file")
      if (.not. found_lastfile) call FatalErrorMessage("Could not find 'Last output file' in RunInfo.txt file")
      if (.not. found_tstart) call FatalErrorMessage("Could not find 't start' in RunInfo.txt file")

      if (InfoVersion /= RunParams%version) call WarningMessage("Results associated with RunInfo.txt were made using a different Kestrel version to your current installation")

      RunParams%FirstOut = LastFile
      RunParams%DeltaT = FileTimeStep
      RunParams%tstart = tstart

   end subroutine ReadRunInfoToRunParams

end module restart_module
