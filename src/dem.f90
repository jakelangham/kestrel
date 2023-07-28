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


! This module contains routines for including topographic data
! in the model.  This can be in the form of formulae for the 
! topographic elevation, or through external raster data that
! is read into the model.  Using external data requires the GDAL
! library that is accessed through RasterData.cpp
module dem_module

   use, intrinsic :: iso_c_binding
   use set_precision_module, only: wp, pi
   use runsettings_module, only: RunSet, FortranRasterData
   use grid_module, only: GridType
   use utilities_module, only: pair
   use messages_module, only: InfoMessage, WarningMessage, FatalErrorMessage
   use varstring_module, only: varString
   use interpolate_2d_module, only: Bicubic
   use geotiffread_module, only: RasterInfo, GetRasterSection, BuildDEMVRT_srtm, BuildDEMVRT_raster
   use topog_funcs_module, only: TopogFunc

   implicit none

   private
   public :: LoadDEM
   public :: GetHeights

contains

   ! If the user has specified a topography that requires an external datafile,
   ! such as a DEM, we try to load it in preparation for running a simulation.
   ! Note that we need to have read the input file (to populate RunParams) and
   ! initialised the grid before calling this routine.
   subroutine LoadDEM(RunParams, grid)
      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(GridType), intent(inout) :: grid

      select case (RunParams%Topog%s)
         case ('srtm', 'SRTM', 'raster', 'Raster')
            call InfoMessage("Preparing topographic data")
            call PrepareDEM(RunParams, grid%deltaX, grid%deltaY)
      end select
 
   end subroutine LoadDEM

   subroutine PrepareDEM(RunParams, deltaX, deltaY, interp_flag_in)

      implicit none

      type(RunSet), intent(inout) :: RunParams
      real(kind=wp), intent(in) :: deltaX
      real(kind=wp), intent(in) :: deltaY
      integer(kind=c_int), intent(in), optional :: interp_flag_in

      integer(kind=c_int) :: interp_flag = int(2, kind=c_int)

      type(pair) :: res

      real(kind=c_double) :: minE, maxE, minN, maxN

      character(len=1,kind=c_char) :: cpath(len_trim(RunParams%OutDir%s)+1)
      character(len=1,kind=c_char) :: csrtmpath(len_trim(RunParams%SRTMPath%s)+1)

      logical(kind=C_BOOL) :: raster

      type(varString) :: geotif
      character(len=1,kind=c_char) :: cgeotif(len_trim(trim(RunParams%demPath%s))+1+len_trim(RunParams%RasterFile%s) + 1)
      logical FileExists

      if (present(interp_flag_in)) interp_flag = interp_flag_in

      cpath = RunParams%OutDir%to_cstring()
      csrtmpath = RunParams%SRTMPath%to_cstring()

      call DomainExtent(RunParams, deltaX, deltaY, minE, maxE, minN, maxN, res, pad=.True.)

      select case (RunParams%Topog%s)
       case ('srtm')
         raster = .FALSE.

         call BuildDEMVRT_srtm(cpath, &
                          csrtmpath, &
                          RunParams%utmEPSG, &
                          minE, maxE, minN, maxN, &
                          real(res%first, kind=c_double), real(res%second, kind=c_double))

       case ('raster','dem')
         geotif = RunParams%demPath + RunParams%RasterFile
         cgeotif = geotif%to_cstring()

         inquire(file=geotif%s, exist=FileExists)
         if (.not.FileExists) then
            call FatalErrorMessage("geotiff file "//trim(geotif%s)//" does not exist")
         end if
         raster = .TRUE.

         call BuildDEMVRT_raster(cpath, &
                          csrtmpath, &
                          cgeotif, &
                          RunParams%utmEPSG, &
                          RunParams%EmbedRaster, &
                          minE, maxE, minN, maxN, &
                          real(res%first, kind=c_double), real(res%second, kind=c_double))

      end select
      
   end subroutine PrepareDEM

   ! Get the scalar field from the raster and fill into the bandSection array.
   ! If tile lies outside the raster, set to defaultValue and outsideRaster = true
   subroutine GetRasterData(path, fname, tileNW, tileSE, rdata, &
      defaultValue, &
      bandSection, rdX, rdY, rOX, rOY, &
      outsideRaster)
      implicit none

      type(varString), intent(in) :: path
      type(varString), intent(in) :: fname
      type(pair), intent(in) :: tileNW
      type(pair), intent(in) :: tileSE
      type(FortranRasterData), intent(in) :: rdata
      real(kind=wp), intent(in) :: defaultValue ! fill anything outside raster with this
      real(kind=wp), intent(out), dimension(:,:), allocatable :: bandSection
      real(kind=wp), intent(out) :: rdX, rdY, rOX, rOY
      logical, intent(inout) :: outsideRaster

      real(kind=wp) :: imgx, imgy
      integer(kind=c_int) :: xoff, yoff, xsize, ysize
      type(FortranRasterData) :: rasterSection
      integer :: i, j

      rDX = rdata%pixelWidth
      rDY = rdata%pixelHeight
      rOX = rdata%originX
      rOY = rdata%originY

      imgx = (tileNW%first - rOX) / rdX
      imgy = (tileNW%second - rOY) / rdY

      xoff = floor(imgx, kind=c_int) - 2
      yoff = floor(imgy, kind=c_int) - 2

      xsize = ceiling((tileSE%first - tileNW%first)/rdX, kind=c_int) + 5
      ysize = ceiling((tileSE%second - tileNW%second)/rdY, kind=c_int) + 5

      call GetRasterSection(path, fname, xoff, yoff, xsize, ysize, rasterSection)

      rdX = rasterSection%pixelWidth
      rdY = rasterSection%pixelHeight
      rOX = rasterSection%originX
      rOY = rasterSection%originY

      allocate(bandSection(xsize, ysize))
      do i = 1, int(xsize)
         do j = 1, int(ysize)
            bandSection(i, j) = rasterSection%values((j - 1) *  &
               rasterSection%RasterXSize + i)
            if (bandSection(i, j) == rasterSection%nodata) then
               bandSection(i, j) = defaultValue
               outsideRaster = .true.
            end if
         end do
      end do
   end subroutine GetRasterData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine TileHeightData(RunParams, deltaX, deltaY, xtile, ytile, b0)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: deltaX
      real(kind=wp), intent(in) :: deltaY
      real(kind=wp), dimension(:), intent(in) :: xtile
      real(kind=wp), dimension(:), intent(in) :: ytile
      real(kind=wp), dimension(:,:), intent(out) :: b0(RunParams%nXpertile+1,RunParams%nYpertile+1)

      type(varString) :: path
      type(varString) :: fname
      logical :: FileExists

      type(FortranRasterData) :: rasterFull

      integer :: nXpertile, nYpertile

      type(pair) :: tileNW
      type(pair) :: tileSE
      type(pair) :: distxy

      integer :: i, j, k
      integer :: ii, jj

      real(kind=wp) :: imgx, imgy
      real(kind=wp), dimension(:) :: hcorner(4)

      real(kind=wp) :: rasterDeltaX, rasterDeltaY
      real(kind=wp) :: rasterOX, rasterOY
      real(kind=wp), dimension(:,:), allocatable :: Elev

      logical :: outsideRaster

      outsideRaster = .false.

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      path = RunParams%OutDir
      fname = varString("DEM.vrt")

      inquire(file=path%s // fname%s, exist=FileExists)
      if (.not.FileExists) then
         call FatalErrorMessage("dem file " // path%s // fname%s // " not found")
         error stop (1)
      end if

      call RasterInfo(path=path, filename=fname, rasterF=rasterFull)
      rasterDeltaX = rasterFull%pixelWidth
      rasterDeltaY = rasterFull%pixelHeight
      rasterOX = rasterFull%originX
      rasterOY = rasterFull%originY

      tileNW%first  = RunParams%centerUTM%first + (xtile(1) - 0.5 * deltaX)
      tileNW%second = RunParams%centerUTM%second + (ytile(nYpertile+1) + 0.5 * deltaY)

      tileSE%first  = RunParams%centerUTM%first + (xtile(nXpertile+1) + 0.5 * deltaX)
      tileSE%second = RunParams%centerUTM%second + (ytile(1) - 0.5 * deltaY)

      call GetRasterData(path=path, fname=fname, tileNW=tileNW, tileSE=tileSE, rdata=rasterFull, &
         defaultValue=0.0_wp, &
         bandSection=Elev, &
         rdX=rasterDeltaX, rdY=rasterDeltaY, &
         rOX=rasterOX, rOY=rasterOY, outsideRaster=outsideRaster)

      if (size(Elev)>1) then
         do i=1,nXpertile+1
            do j=1,nYpertile+1
               do ii=1,2
                  do jj=1,2
                     k = jj+2*(ii-1)
                     
                     distxy%first  = RunParams%centerUTM%first  + xtile(i)+(real(ii-1,kind=wp)-0.5_wp)*deltaX
                     distxy%second = RunParams%centerUTM%second + ytile(j)+(real(jj-1,kind=wp)-0.5_wp)*deltaY
                     
                     imgx = (distxy%first  - rasterOX)/rasterDeltaX+1
                     imgy = (distxy%second - rasterOY)/rasterDeltaY+1

                     hcorner(k) = Bicubic(imgx,imgy,Elev)

                  end do
               end do
               b0(i,j) = 0.25_wp*sum(hcorner)
            end do
         end do
      else

         do i=1,nXpertile+1
            do j=1,nYpertile+1
               b0(i,j) = -999_wp
            end do
         end do

      end if

      deallocate(Elev)

   end subroutine TileHeightData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine GetHeights(RunParams, grid, tID)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), intent(inout) :: grid
      integer, intent(in) :: tID

      real(kind=wp), dimension(:,:) :: b0(RunParams%nXpertile+1,RunParams%nYpertile+1)

      integer :: i, ib0

      ib0 = RunParams%Vars%b0

      call GetTileHeights(RunParams, grid%deltaX, grid%deltaY, &
         grid%tileContainer(tID)%x_vertex, grid%tileContainer(tID)%y_vertex, b0)
      
      grid%tileContainer(tID)%b0(:,:) = b0(:,:)

      ! These are cell interface values, stored for computing numerical fluxes,
      ! in particular the hydrostatic pressure which depends on w - b
      if (.not. RunParams%isOneD) then
         do i = 1, RunParams%nYpertile
            grid%tileContainer(tID)%uPlusX(ib0,:,i) = 0.5_wp * (b0(:,i) + b0(:,i+1))
            grid%tileContainer(tID)%uMinusX(ib0,:,i) = 0.5_wp * (b0(:,i) + b0(:,i+1))
         end do
         do i = 1, RunParams%nXpertile
            grid%tileContainer(tID)%uPlusY(ib0,i,:) = 0.5_wp * (b0(i,:) + b0(i+1,:))
            grid%tileContainer(tID)%uMinusY(ib0,i,:) = 0.5_wp * (b0(i,:) + b0(i+1,:))
         end do
      else
         grid%tileContainer(tID)%uPlusX(ib0,:,1) = b0(:,1)
         grid%tileContainer(tID)%uMinusX(ib0,:,1) = b0(:,1)
      end if
   end subroutine GetHeights

   subroutine GetTileHeights(RunParams, deltaX, deltaY, x, y, b0)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: deltaX
      real(kind=wp), intent(in) :: deltaY
      real(kind=wp), dimension(:), intent(in) :: x(RunParams%nXpertile+1)
      real(kind=wp), dimension(:), intent(in) :: y(RunParams%nYpertile+1)
      real(kind=wp), dimension(:,:), intent(out) :: b0(RunParams%nXpertile+1,RunParams%nYpertile+1)

      select case (RunParams%Topog%s)

       case ('dem', 'raster', 'srtm')
         call TileHeightData(RunParams, deltaX, deltaY, x, y, b0)
       
       case default
         call TopogFunc(RunParams, x, y, b0)

      end select

   end subroutine GetTileHeights

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine DomainExtent(RunParams, deltaX, deltaY, minE, maxE, minN, maxN, res, pad)

      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: deltaX
      real(kind=wp), intent(in) :: deltaY
      real(kind=c_double), intent(out) :: minE, maxE, minN, maxN
      type(pair), intent(out) :: res
      logical, intent(in), optional :: pad

      type(pair) :: distxy

      logical :: pad_domain

      if (present(pad)) then
         pad_domain = pad
      else
         pad_domain = .FALSE.
      end if

      ! Get domain box NW point -- if pad==True add an additional tile buffer on each side
      if (pad_domain) then
         distxy%first  = (-0.5_wp*(RunParams%nXtiles+2)*RunParams%Xtilesize - 0.5_wp*deltaX)
         distxy%second = (0.5_wp*(RunParams%nYtiles+2)*RunParams%Ytilesize + 0.5_wp*deltaY)
      else
         distxy%first =  RunParams%centerUTM%first  + (-0.5_wp*(RunParams%nXtiles)*RunParams%Xtilesize - 0.5_wp*deltaX)
         distxy%second = RunParams%centerUTM%second + (0.5_wp*(RunParams%nYtiles)*RunParams%Ytilesize + 0.5_wp*deltaY)
      end if

      minE  = RunParams%centerUTM%first + distxy%first
      maxN = RunParams%centerUTM%second + distxy%second

      ! Get domain box SE point -- if pad==True add an additional tile buffer
      if (pad_domain) then
         distxy%first = (0.5_wp*(RunParams%nXtiles+2)*RunParams%Xtilesize + 0.5_wp*deltaX)
         distxy%second = (-0.5_wp*(RunParams%nYtiles+2)*RunParams%Ytilesize - 0.5_wp*deltaY)
      else
         distxy%first = (0.5_wp*(RunParams%nXtiles)*RunParams%Xtilesize + 0.5_wp*deltaX)
         distxy%second = (-0.5_wp*(RunParams%nYtiles)*RunParams%Ytilesize - 0.5_wp*deltaY)
      end if

      maxE  = RunParams%centerUTM%first + distxy%first
      minN = RunParams%centerUTM%second + distxy%second

      res%first = deltaX
      res%second = deltaY

   end subroutine DomainExtent

end module dem_module
