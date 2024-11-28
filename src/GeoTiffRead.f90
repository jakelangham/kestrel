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


! This module defines bindings to C/C++ functions to read
! geotiff files using the GDAL library.
! See also RasterData.cpp and RasterData.h
! See also FortranRasterData in RunSettings.f90
module geotiffread_module

   use, intrinsic :: iso_C_binding
   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage
   use runsettings_module, only: RunSet, FortranRasterData, RasterData
   use varstring_module, only: varString

   implicit none

   private
   public :: RasterInfo
   public :: GetRasterSection
   public :: BuildDEMVRT_raster
   public :: BuildDEMVRT_srtm

! Interfaces to access C memory allocation from Fortran
   interface
      subroutine MallocDouble(ptr, n) bind(C,name="MallocDouble")
         use iso_C_binding
         implicit none
         type(c_ptr), intent(out) :: ptr
         integer(c_long), intent(in) :: n
      end subroutine MallocDouble

      subroutine FreeDouble(ptr) bind(C,name="FreeDouble")
         use iso_C_binding
         implicit none
         type(c_ptr), intent(out) :: ptr
      end subroutine FreeDouble

      subroutine GeoTiffInfo(cgeotifname, raster) bind(C,name="GeoTiffInfo")
         use iso_C_binding
         use runsettings_module, only : RasterData
         implicit none
         type(RasterData) :: raster
         character(len=1,kind=c_char), intent(in) :: cgeotifname(*)
      end subroutine GeoTiffInfo

      subroutine GeoTiffArraySectionRead(cpath, cgeotifname, raster, xoff, yoff, xsize, ysize) bind(C,name="GeoTiffArraySectionRead")
         use iso_C_binding
         use runsettings_module, only : RasterData
         implicit none
         character(len=1,kind=c_char), intent(in) :: cpath(*)
         character(len=1,kind=c_char), intent(in) :: cgeotifname(*)
         type(RasterData) :: raster
         integer(c_int) :: xoff
         integer(c_int) :: yoff
         integer(c_int) :: xsize
         integer(c_int) :: ysize
      end subroutine GeoTiffArraySectionRead

      subroutine BuildDEMVRT_raster(cpath, csrtmpath, cgeotifname, utmEPSG, embed, minE, maxE, minN, maxN, xres, yres) bind(C,name="BuildDEMVRT_raster")
         use iso_C_binding
         implicit none
         character(len=1,kind=c_char), intent(in) :: cpath(*)
         character(len=1,kind=c_char), intent(in) :: csrtmpath(*)
         character(len=1,kind=c_char), intent(in) :: cgeotifname(*)
         integer(kind=c_int), intent(in) :: utmEPSG
         logical(kind=c_bool), intent(in) :: embed
         real(kind=c_double), intent(in) :: minE
         real(kind=c_double), intent(in) :: maxE
         real(kind=c_double), intent(in) :: minN
         real(kind=c_double), intent(in) :: maxN
         real(kind=c_double), intent(in) :: xres
         real(kind=c_double), intent(in) :: yres
      end subroutine BuildDEMVRT_raster

      subroutine BuildDEMVRT_srtm(cpath, csrtmpath, utmEPSG, minE, maxE, minN, maxN, xres, yres) bind(C,name="BuildDEMVRT_srtm")
        use iso_C_binding
        implicit none
        character(len=1,kind=c_char), intent(in) :: cpath(*)
        character(len=1,kind=c_char), intent(in) :: csrtmpath(*)
        integer(kind=c_int), intent(in) :: utmEPSG
        real(kind=c_double), intent(in) :: minE
        real(kind=c_double), intent(in) :: maxE
        real(kind=c_double), intent(in) :: minN
        real(kind=c_double), intent(in) :: maxN
        real(kind=c_double), intent(in) :: xres
        real(kind=c_double), intent(in) :: yres
     end subroutine BuildDEMVRT_srtm

   end interface

contains

   ! Get information about a raster file
   ! Inputs: path - the path of the raster file
   !         filename - the name of the geotif raster file
   ! Returns: rasterF - a FortranRasterData object containing data about the raster file
   subroutine RasterInfo(path, filename, rasterF)

      implicit none
      type(varString), intent(in) :: path
      type(varString), intent(in) :: filename
      type(FortranRasterData), intent(out) :: rasterF

      character(len=1,kind=c_char) :: cpath(path%len()+1)
      character(len=1,kind=c_char) :: cfilename(path%len()+filename%len()+1)

      type(varString) :: filename_full

      type(RasterData) :: raster

      real(kind=wp),pointer :: a(:)
      
      cpath = path%to_cstring()
      filename_full = path + filename
      cfilename = filename_full%to_cstring()

      call GeoTiffInfo(cfilename, raster)

      rasterF%RasterName = filename

      rasterF%ReadSuccess = raster%read_success
      rasterF%RasterXSize = raster%x_size
      rasterF%RasterYSize = raster%y_size
      rasterF%RasterSize = raster%size
      rasterF%pixelWidth = raster%pixel_width
      rasterF%pixelHeight = raster%pixel_height
      rasterF%originX = raster%origin_x
      rasterF%originY = raster%origin_y
      rasterF%rowRotation = raster%row_rotation
      rasterF%columnRotation = raster%column_rotation
      rasterF%EPSG_code = raster%EPSG_code
      rasterF%nodata = raster%no_data
      rasterF%latlon = raster%latlon

      call c_f_pointer(raster%NW,a,[2])
      rasterF%NW%first = a(1)
      rasterF%NW%second = a(2)
      call c_f_pointer(raster%NE,a,[2])
      rasterF%NE%first = a(1)
      rasterF%NE%second = a(2)
      call c_f_pointer(raster%SE,a,[2])
      rasterF%SE%first = a(1)
      rasterF%SE%second = a(2)
      call c_f_pointer(raster%SW,a,[2])
      rasterF%SW%first = a(1)
      rasterF%SW%second = a(2)

      return
   end subroutine RasterInfo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Get data from a section of a raster file
   ! Inputs: path - the path of the raster file
   !         filename - the name of the geotif raster file
   !         Xoff - number of columns to offset from the origin (upper left corner)
   !         Yoff - number of rows to offset from the origin (upper left corner)
   !         Xsize - number of columns to read
   !         Ysize - number of rows to read
   ! Returns: rasterF - a FortranRasterData object containing data about the raster file
   subroutine GetRasterSection(path, filename, Xoff, Yoff, Xsize, Ysize, rasterF)

      implicit none
      type(varString), intent(in) :: path
      type(varString), intent(in) :: filename
      integer(kind=c_int), intent(in) :: Xoff
      integer(kind=c_int), intent(in) :: Yoff
      integer(kind=c_int), intent(in) :: Xsize
      integer(kind=c_int), intent(in) :: Ysize
      type(FortranRasterData), intent(out) :: rasterF

      character(len=1,kind=c_char) :: cpath(path%len()+1)
      character(len=1,kind=c_char) :: cfilename(filename%len()+1)

      type(RasterData) :: raster

      real(kind=wp),pointer :: a(:)

      logical FileExists

      cpath = path%to_cstring()
      cfilename = filename%to_cstring()

      inquire(file=path%s//"/"//filename%s, exist=FileExists)
      if (.not.FileExists) then
         call FatalErrorMessage("DEM file " // path%s//"/"//filename%s // " does not exist")
      end if

      call GeoTiffArraySectionRead(cpath, cfilename, raster, Xoff, Yoff, Xsize, Ysize)

      rasterF%RasterName = filename
      rasterF%ReadSuccess = raster%read_success

      if (rasterF%ReadSuccess == -1) call FatalErrorMessage("Could not read raster file " // rasterF%RasterName%s)
      if (rasterF%ReadSuccess == -2) call FatalErrorMessage("No-data value found in raster file " // rasterF%RasterName%s // ". Check the DEM, fill no-data, and restart simulation.")

      rasterF%RasterXSize = raster%x_size
      rasterF%RasterYSize = raster%y_size
      rasterF%RasterSize = raster%size
      rasterF%pixelWidth = raster%pixel_width
      rasterF%pixelHeight = raster%pixel_height
      rasterF%originX = raster%origin_x
      rasterF%originY = raster%origin_y
      rasterF%EPSG_code = raster%EPSG_code
      rasterF%nodata = raster%no_data

      call c_f_pointer(raster%NW,a,[2])
      rasterF%NW%first = a(1)
      rasterF%NW%second = a(2)
      call c_f_pointer(raster%NE,a,[2])
      rasterF%NE%first = a(1)
      rasterF%NE%second = a(2)
      call c_f_pointer(raster%SE,a,[2])
      rasterF%SE%first = a(1)
      rasterF%SE%second = a(2)
      call c_f_pointer(raster%SW,a,[2])
      rasterF%SW%first = a(1)
      rasterF%SW%second = a(2)

      call c_f_pointer(raster%values,rasterF%values,[XSize*YSize])

      return
   end subroutine GetRasterSection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module geotiffread_module
