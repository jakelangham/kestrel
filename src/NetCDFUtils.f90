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


! This module contains routines that are used to output to NetCDF.
! NetCDF is an optional feature, and requires netCDF libraries.
! If configuration does not specify HAVE_NETCDF then this module is
! not compiled.
!
! The module routines define overloaded routines to define and add variables
! and attributes to a NetCDF file.  There is also a routine (set_nc_crs) that
! writes a coordinate reference system to the NetCDF for georeferencing.
!
! This module is used in Output.f90 if HAVE_NETCDF is set.
#if HAVE_NETCDF4
module netcdf_utils_module

   use, intrinsic :: iso_c_binding
   use netcdf
   use set_precision_module, only: wp
   use messages_module, only: FatalErrorMessage
   use utilities_module, only: Int2String
   use varstring_module, only: varString
   use utm_module, only: ZoneNumberToCentralLongitude

   implicit none

   private
   public :: put_nc_var
   public :: put_nc_att
   public :: define_nc_var_real
   public :: define_nc_var_int
   public :: set_nc_crs
   public :: get_nc_varid
   public :: get_nc_var
   public :: get_nc_var_scalar_int
   public :: get_nc_var_scalar_real
   public :: get_nc_att
   public :: nc_open_read
   public :: nc_get_tiles
   public :: handle_err

   interface put_nc_var
      module procedure :: put_nc_var_scalar_real, put_nc_var_scalar_int, put_nc_var_scalar_str, &
        put_nc_var_vec_real, put_nc_var_vec_int, put_nc_var_vec_str, &
        put_nc_var_arr_real, put_nc_var_arr_int
   end interface

   interface put_nc_att
      module procedure :: put_nc_att_global_int, put_nc_att_global_int_vec, &
         put_nc_att_global_real, put_nc_att_global_real_vec, &
         put_nc_att_global_str, put_nc_att_global_varString, &
         put_nc_att_var_int, put_nc_att_var_int_vec, &
         put_nc_att_var_real, put_nc_att_var_real_vec, &
         put_nc_att_var_str, put_nc_att_var_varString
   end interface put_nc_att

   interface define_nc_var_real
      module procedure :: define_nc_var_scalar_real, define_nc_var_nonscalar_real
   end interface
   interface define_nc_var_int
      module procedure :: define_nc_var_scalar_int, define_nc_var_nonscalar_int
   end interface

   interface get_nc_var
      module procedure :: get_nc_var_scalar_real, get_nc_var_scalar_int, &
         get_nc_var_vec_real, get_nc_var_vec_int, &
         get_nc_var_arr_real, get_nc_var_arr_int
   end interface

   interface get_nc_att
      module procedure :: get_nc_att_int, get_nc_att_real, &
         get_nc_att_int_vec, get_nc_att_real_vec, &
         get_nc_att_str
   end interface

contains

   subroutine put_nc_att_global_int(ncid, name, value)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, NF90_GLOBAL, name, int(value, kind=4))

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att '//name)

      return
   end subroutine put_nc_att_global_int
  
   subroutine put_nc_att_global_int_vec(ncid, name, value)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, dimension(:), intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, NF90_GLOBAL, name, int(value, kind=4))

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att '//name)

      return
   end subroutine put_nc_att_global_int_vec
  
   subroutine put_nc_att_global_real(ncid, name, value)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      real(kind=wp), intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, NF90_GLOBAL, name, real(value, kind=8))

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att '//name)

      return
   end subroutine put_nc_att_global_real
  
   subroutine put_nc_att_global_real_vec(ncid, name, value)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      real(kind=wp), dimension(:), intent(in) :: value
      integer :: nc_status

      nc_status = nf90_put_att(ncid, NF90_GLOBAL, name, real(value, kind=8))

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att '//name)

      return
   end subroutine put_nc_att_global_real_vec
  
   subroutine put_nc_att_global_str(ncid, name, value)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, NF90_GLOBAL, name, value)

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att '//name)

      return
   end subroutine put_nc_att_global_str

   subroutine put_nc_att_global_varString(ncid, name, value)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    type(varString), intent(in) :: value

    integer :: nc_status

    nc_status = nf90_put_att(ncid, NF90_GLOBAL, name, value%s)

    if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att '//name)

    return
 end subroutine put_nc_att_global_varString
  
   subroutine put_nc_att_var_int(ncid, varid, name, value)
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      integer, intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, varid, name, value)

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att_var '//name)

      return
   end subroutine put_nc_att_var_int
  
   subroutine put_nc_att_var_int_vec(ncid, varid, name, value)
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      integer, dimension(:), intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, varid, name, value)

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att_var '//name)

      return
   end subroutine put_nc_att_var_int_vec
  
   subroutine put_nc_att_var_real(ncid, varid, name, value)
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      real(kind=wp), intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, varid, name, real(value, kind=c_double))

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att_var '//name)

      return
   end subroutine put_nc_att_var_real
  
   subroutine put_nc_att_var_real_vec(ncid, varid, name, value)
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      real(kind=wp), dimension(:), intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, varid, name, real(value, kind=c_double))

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att_var '//name)

      return
   end subroutine put_nc_att_var_real_vec
  
   subroutine put_nc_att_var_str(ncid, varid, name, value)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: varid
      character(len=*), intent(in) :: value

      integer :: nc_status

      nc_status = nf90_put_att(ncid, varid, name, value)

      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att_var '//name)

      return
   end subroutine put_nc_att_var_str

   subroutine put_nc_att_var_varString(ncid, varid, name, value)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: varid
    type(varString), intent(in) :: value

    integer :: nc_status

    nc_status = nf90_put_att(ncid, varid, name, value%s)

    if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_put_att_var '//name)

    return
 end subroutine put_nc_att_var_varString
  
   subroutine define_nc_var_scalar_real(ncid, name, var_id, fill_value, deflate, &
                                          units, long_name, standard_name, axis, coordinates, grid_mapping)
  
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: var_id
      real(kind=wp), intent(in), optional :: fill_value
      integer, intent(in), optional :: deflate
      character(len=*), intent(in), optional :: units
      character(len=*), intent(in), optional :: long_name
      character(len=*), intent(in), optional :: standard_name
      character(len=*), intent(in), optional :: axis
      character(len=*), intent(in), optional :: coordinates
      character(len=*), intent(in), optional :: grid_mapping

      integer :: nc_status

      nc_status = nf90_def_var(ncid, name, NF90_DOUBLE, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var '//name)
      if (present(fill_value)) then
         nc_status = nf90_def_var_fill(ncid, var_id, 0, real(fill_value, kind=c_double)) ! Define fill value to be NaN, but fill values will not be written for the variable
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var_fill '//name)
         call put_nc_att(ncid, var_id, '_FillValue', fill_value)
      end if

      if (present(deflate)) then
         nc_status = nf90_def_var_deflate(ncid, var_id, 0, 1, deflate)
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var_deflate '//name)
      end if

      if (present(long_name)) call put_nc_att(ncid, var_id, 'long_name', long_name)
      if (present(standard_name)) then
         call put_nc_att(ncid, var_id, 'standard_name', standard_name)
      else
         call put_nc_att(ncid, var_id, 'standard_name', name)
      end if
      if (present(units)) call put_nc_att(ncid, var_id, 'units', units)
      if (present(axis)) call put_nc_att(ncid, var_id, 'axis', axis)
      if (present(coordinates)) call put_nc_att(ncid, var_id, 'coordinates', coordinates)
      if (present(grid_mapping)) call put_nc_att(ncid, var_id, 'grid_mapping', grid_mapping)

      return
   end subroutine define_nc_var_scalar_real
  
   subroutine define_nc_var_nonscalar_real(ncid, name, size_ids, var_id, fill_value, deflate, &
                                             units, long_name, standard_name, axis, coordinates, grid_mapping)
  
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, dimension(:), intent(in) :: size_ids
      integer, intent(out) :: var_id
      real(kind=wp), intent(in), optional :: fill_value
      integer, intent(in), optional :: deflate
      character(len=*), intent(in), optional :: units
      character(len=*), intent(in), optional :: long_name
      character(len=*), intent(in), optional :: standard_name
      character(len=*), intent(in), optional :: axis
      character(len=*), intent(in), optional :: coordinates
      character(len=*), intent(in), optional :: grid_mapping

      integer :: nc_status

      nc_status = nf90_def_var(ncid, name, NF90_DOUBLE, size_ids, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var '//name)
      if (present(fill_value)) then
         nc_status = nf90_def_var_fill(ncid, var_id, 0, real(fill_value, kind=c_double)) ! Define fill value to be NaN, but fill values will not be written for the variable
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var_fill '//name)
      end if

      if (present(deflate)) then
         nc_status = nf90_def_var_deflate(ncid, var_id, 0, 1, deflate)
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var_deflate '//name)
      end if

      if (present(long_name)) call put_nc_att(ncid, var_id, 'long_name', long_name)
      if (present(standard_name)) then
         call put_nc_att(ncid, var_id, 'standard_name', standard_name)
      else
         call put_nc_att(ncid, var_id, 'standard_name', name)
      end if
      if (present(units)) call put_nc_att(ncid, var_id, 'units', units)
      if (present(axis)) call put_nc_att(ncid, var_id, 'axis', axis)
      if (present(coordinates)) call put_nc_att(ncid, var_id, 'coordinates', coordinates)
      if (present(grid_mapping)) call put_nc_att(ncid, var_id, 'grid_mapping', grid_mapping)

      return
   end subroutine define_nc_var_nonscalar_real
  
   subroutine define_nc_var_scalar_int(ncid, name, var_id, fill_value, deflate, &
                                         units, long_name, standard_name, axis, coordinates, grid_mapping)
  
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: var_id
      integer, intent(in), optional :: fill_value
      integer, intent(in), optional :: deflate
      character(len=*), intent(in), optional :: units
      character(len=*), intent(in), optional :: long_name
      character(len=*), intent(in), optional :: standard_name
      character(len=*), intent(in), optional :: axis
      character(len=*), intent(in), optional :: coordinates
      character(len=*), intent(in), optional :: grid_mapping

      integer :: nc_status

      nc_status = nf90_def_var(ncid, name, NF90_INT, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var '//name)
      if (present(fill_value)) then
         nc_status = nf90_def_var_fill(ncid, var_id, 0, fill_value) ! Define fill value to be NaN, but fill values will not be written for the variable
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var_fill '//name)
         call put_nc_att(ncid, var_id, '_FillValue', fill_value)
      end if

      if (present(deflate)) then
         nc_status = nf90_def_var_deflate(ncid, var_id, 0, 1, deflate)
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var_deflate '//name)
      end if

      if (present(long_name)) call put_nc_att(ncid, var_id, 'long_name', long_name)
      if (present(standard_name)) then
         call put_nc_att(ncid, var_id, 'standard_name', standard_name)
      else
         call put_nc_att(ncid, var_id, 'standard_name', name)
      end if
      if (present(units)) call put_nc_att(ncid, var_id, 'units', units)
      if (present(axis)) call put_nc_att(ncid, var_id, 'axis', axis)
      if (present(coordinates)) call put_nc_att(ncid, var_id, 'coordinates', coordinates)
      if (present(grid_mapping)) call put_nc_att(ncid, var_id, 'grid_mapping', grid_mapping)

      return
   end subroutine define_nc_var_scalar_int
  
   subroutine define_nc_var_nonscalar_int(ncid, name, size_ids, var_id, fill_value, deflate, &
                                            units, long_name, standard_name, axis, coordinates, grid_mapping)
  
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, dimension(:), intent(in) :: size_ids
      integer, intent(out) :: var_id
      integer, intent(in), optional :: fill_value
      integer, intent(in), optional :: deflate
      character(len=*), intent(in), optional :: units
      character(len=*), intent(in), optional :: long_name
      character(len=*), intent(in), optional :: standard_name
      character(len=*), intent(in), optional :: axis
      character(len=*), intent(in), optional :: coordinates
      character(len=*), intent(in), optional :: grid_mapping

      integer :: nc_status

      nc_status = nf90_def_var(ncid, name, NF90_INT, size_ids, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var '//name)
      if (present(fill_value)) then
         nc_status = nf90_def_var_fill(ncid, var_id, 0, fill_value) ! Define fill value to be NaN, but fill values will not be written for the variable
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var_fill '//name)
      end if

      if (present(deflate)) then
         nc_status = nf90_def_var_deflate(ncid, var_id, 0, 1, deflate)
         if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_def_var_deflate '//name)
      end if

      if (present(long_name)) call put_nc_att(ncid, var_id, 'long_name', long_name)
      if (present(standard_name)) then
         call put_nc_att(ncid, var_id, 'standard_name', standard_name)
      else
         call put_nc_att(ncid, var_id, 'standard_name', name)
      end if
      if (present(units)) call put_nc_att(ncid, var_id, 'units', units)
      if (present(axis)) call put_nc_att(ncid, var_id, 'axis', axis)
      if (present(coordinates)) call put_nc_att(ncid, var_id, 'coordinates', coordinates)
      if (present(grid_mapping)) call put_nc_att(ncid, var_id, 'grid_mapping', grid_mapping)

      return
   end subroutine define_nc_var_nonscalar_int
  
   subroutine get_nc_varid(ncid, name, varid)
  
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: varid

      integer :: nc_status

      nc_status = nf90_inq_varid(ncid, name, varid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, 'nf90_inq_varid '//name)

      return
   end subroutine get_nc_varid
  
   subroutine set_nc_crs(ncid, UTM_zone_number, hemisphere, EPSG, crs_id, include_wkt)
      integer, intent(in) :: ncid
      integer, intent(in) :: UTM_zone_number
      character(len=1), intent(in) :: hemisphere
      integer, intent(in) :: EPSG
      integer, intent(out) :: crs_id
      logical, intent(in), optional :: include_wkt

      type(varString) :: crs_wkt

      logical :: wkt

      wkt = .FALSE.
      if (present(include_wkt)) wkt = include_wkt

      call define_nc_var_int(ncid, "crs", crs_id)
      call put_nc_att(ncid, crs_id, "grid_mapping_name", 'transverse_mercator')
      call put_nc_att(ncid, crs_id, "false_easting", 500000)
      if (hemisphere=='N') then
         call put_nc_att(ncid, crs_id, "false_northing", 0)
      else
         call put_nc_att(ncid, crs_id, "false_northing", 10000000)
      end if
      call put_nc_att(ncid, crs_id, "latitude_of_projection_origin", 0)
      call put_nc_att(ncid, crs_id, "longitude_of_central_meridian", ZoneNumberToCentralLongitude(UTM_zone_number))
      call put_nc_att(ncid, crs_id, "scale_factor_at_central_meridian", 0.9996_wp)

      if (wkt) then
         crs_wkt = varString('')
         crs_wkt = crs_wkt + 'PROJCRS["WGS 84 / UTM zone '+Int2String(UTM_zone_number) + hemisphere + '",\n'
         crs_wkt = crs_wkt + '  BASEGEOGCRS["WGS 84",\n'
         crs_wkt = crs_wkt + '    DATUM["World Geodetic System 1984",\n'
         crs_wkt = crs_wkt + '      ELLIPSOID["WGS 84",6378137,298.257223563,\n'
         crs_wkt = crs_wkt + '        LENGTHUNIT["metre",1]]],\n'
         crs_wkt = crs_wkt + '    PRIMEM["Greenwich",0,\n'
         crs_wkt = crs_wkt + '      ANGLEUNIT["degree",0.0174532925199433]],\n'
         crs_wkt = crs_wkt + '    ID["EPSG",4326]],\n'
         crs_wkt = crs_wkt + '  CONVERSION["UTM zone "'+Int2String(UTM_zone_number) + hemisphere + '",\n'
         crs_wkt = crs_wkt + '    METHOD["Transverse Mercator",\n'
         crs_wkt = crs_wkt + '      ID["EPSG",9807]],\n'
         crs_wkt = crs_wkt + '    PARAMETER["Latitude of natural origin",0,\n'
         crs_wkt = crs_wkt + '      ANGLEUNIT["degree",0.0174532925199433],\n'
         crs_wkt = crs_wkt + '      ID["EPSG",8801]],\n'
         crs_wkt = crs_wkt + '    PARAMETER["Longitude of natural origin",' + Int2String(ZoneNumberToCentralLongitude(UTM_zone_number)) +',\n'
         crs_wkt = crs_wkt + '      ANGLEUNIT["degree",0.0174532925199433],\n'
         crs_wkt = crs_wkt + '      ID["EPSG",8802]],\n'
         crs_wkt = crs_wkt + '    PARAMETER["Scale factor at natural origin",0.9996,\n'
         crs_wkt = crs_wkt + '      SCALEUNIT["unity",1],\n'
         crs_wkt = crs_wkt + '      ID["EPSG",8805]],\n'
         crs_wkt = crs_wkt + '    PARAMETER["False easting",500000,\n'
         crs_wkt = crs_wkt + '      LENGTHUNIT["metre",1],\n'
         crs_wkt = crs_wkt + '      ID["EPSG",8806]],\n'
         if (hemisphere=='N') then
            crs_wkt = crs_wkt + '    PARAMETER["False northing",0,\n'
         else
            crs_wkt = crs_wkt + '    PARAMETER["False northing",10000000,\n'
         end if
         crs_wkt = crs_wkt + '      LENGTHUNIT["metre",1],\n'
         crs_wkt = crs_wkt + '      ID["EPSG",8807]]],\n'
         crs_wkt = crs_wkt + '  CS[Cartesian,2],\n'
         crs_wkt = crs_wkt + '    AXIS["(E)",east,\n'
         crs_wkt = crs_wkt + '      ORDER[1],\n'
         crs_wkt = crs_wkt + '      LENGTHUNIT["metre",1]],\n'
         crs_wkt = crs_wkt + '    AXIS["(N)",north,\n'
         crs_wkt = crs_wkt + '      ORDER[2],\n'
         crs_wkt = crs_wkt + '      LENGTHUNIT["metre",1]],\n'
         crs_wkt = crs_wkt + '  ID["EPSG",'+Int2String(EPSG) + ']]'
         call put_nc_att(ncid, crs_id, "crs_wkt", crs_wkt%s)
      end if

      return
  
   end subroutine set_nc_crs
  
   subroutine put_nc_var_scalar_real(ncid, var_id, val)
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      real(kind=wp), intent(in) :: val
      
      integer :: nc_status
      integer :: nc_status_name
      character(len=NF90_MAX_NAME) :: var_name

      nc_status = nf90_put_var(ncid, var_id, [real(val, kind=c_double)], start=[1], count=[1])
      if (nc_status /= NF90_NOERR) then
         nc_status_name = nf90_inquire_variable(ncid, var_id, name=var_name)
         call handle_err(nc_status, 'nf90_put_var '//trim(var_name))
      end if
      return
   end subroutine put_nc_var_scalar_real
  
   subroutine put_nc_var_scalar_int(ncid, var_id, val)
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      integer, intent(in) :: val
      
      integer :: nc_status
      integer :: nc_status_name
      character(len=NF90_MAX_NAME) :: var_name
   
      nc_status = nf90_put_var(ncid, var_id, [val], start=[1], count=[1])
      if (nc_status /= NF90_NOERR) then
         nc_status_name = nf90_inquire_variable(ncid, var_id, name=var_name)
         call handle_err(nc_status, 'nf90_put_var '//trim(var_name))
      end if
      return
   end subroutine put_nc_var_scalar_int
  
   subroutine put_nc_var_scalar_str(ncid, var_id, val)
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      character(len=*), intent(in) :: val
      
      integer :: nc_status
      integer :: nc_status_name
      character(len=NF90_MAX_NAME) :: var_name
   
      nc_status = nf90_put_var(ncid, var_id, [val])
      if (nc_status /= NF90_NOERR) then
         nc_status_name = nf90_inquire_variable(ncid, var_id, name=var_name)
         call handle_err(nc_status, 'nf90_put_var '//trim(var_name))
      end if
      return
   end subroutine put_nc_var_scalar_str
  
   subroutine put_nc_var_vec_real(ncid, var_id, vals, start, count)
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      real(kind=wp), dimension(:), intent(in) :: vals
      integer, dimension(:), intent(in) :: start(1)
      integer, dimension(:), intent(in) :: count(1)

      integer :: nc_status
      integer :: nc_status_name
      character(len=NF90_MAX_NAME) :: var_name

      nc_status = nf90_put_var(ncid, var_id, real(vals(:), kind=c_double), start=start, count=count)
      if (nc_status /= NF90_NOERR) then
         nc_status_name = nf90_inquire_variable(ncid, var_id, name=var_name)
         call handle_err(nc_status, 'nf90_put_var '//trim(var_name))
      end if
      return
   end subroutine put_nc_var_vec_real
  
   subroutine put_nc_var_vec_int(ncid, var_id, vals, start, count)
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      integer, dimension(:), intent(in) :: vals
      integer, dimension(:), intent(in) :: start(1)
      integer, dimension(:), intent(in) :: count(1)

      integer :: nc_status
      integer :: nc_status_name
      character(len=NF90_MAX_NAME) :: var_name

      nc_status = nf90_put_var(ncid, var_id, vals(:), start=start, count=count)
      if (nc_status /= NF90_NOERR) then
         nc_status_name = nf90_inquire_variable(ncid, var_id, name=var_name)
         call handle_err(nc_status, 'nf90_put_var '//trim(var_name))
      end if
      return
   end subroutine put_nc_var_vec_int
  
   subroutine put_nc_var_vec_str(ncid, var_id, vals, start, count)
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      character(len=*), dimension(:), intent(in) :: vals
      integer, dimension(:), intent(in) :: start(1)
      integer, dimension(:), intent(in) :: count(1)

      integer :: nc_status
      integer :: nc_status_name
      character(len=NF90_MAX_NAME) :: var_name

      nc_status = nf90_put_var(ncid, var_id, vals(:), start=start, count=count)
      if (nc_status /= NF90_NOERR) then
         nc_status_name = nf90_inquire_variable(ncid, var_id, name=var_name)
         call handle_err(nc_status, 'nf90_put_var '//trim(var_name))
      end if
      return
   end subroutine put_nc_var_vec_str
  
   subroutine put_nc_var_arr_real(ncid, var_id, vals, start, count)
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      real(kind=wp), dimension(:, :), intent(in) :: vals
      integer, dimension(:), intent(in) :: start(2)
      integer, dimension(:), intent(in) :: count(2)

      integer :: nc_status
      integer :: nc_status_name
      character(len=NF90_MAX_NAME) :: var_name

      nc_status = nf90_put_var(ncid, var_id, real(vals(:, :), kind=c_double), start=start, count=count)
      if (nc_status /= NF90_NOERR) then
         nc_status_name = nf90_inquire_variable(ncid, var_id, name=var_name)
         call handle_err(nc_status, 'nf90_put_var '//trim(var_name))
      end if
      return
   end subroutine put_nc_var_arr_real
  
   subroutine put_nc_var_arr_int(ncid, var_id, vals, start, count)
      integer, intent(in) :: ncid
      integer, intent(in) :: var_id
      integer, dimension(:, :), intent(in) :: vals
      integer, dimension(:), intent(in) :: start(2)
      integer, dimension(:), intent(in) :: count(2)

      integer :: nc_status
      integer :: nc_status_name
      character(len=NF90_MAX_NAME) :: var_name

      nc_status = nf90_put_var(ncid, var_id, vals(:, :), start=start, count=count)
      if (nc_status /= NF90_NOERR) then
         nc_status_name = nf90_inquire_variable(ncid, var_id, name=var_name)
         call handle_err(nc_status, 'nf90_put_var '//trim(var_name))
      end if
      return
   end subroutine put_nc_var_arr_int

   function nc_open_read(path) result(ncid)
      character(len=*), intent(in) :: path
      integer :: ncid

      integer :: nc_status

      nc_status = nf90_open(path=path, mode=nf90_nowrite, ncid=ncid)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' open '//path)

   end function nc_open_read

   function nc_get_tiles(ncid) result(tiles)
      integer, intent(in) :: ncid
      integer, dimension(:), allocatable :: tiles

      integer :: nTiles
      integer :: nc_status

      nc_status = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'tiles', len=nTiles)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inquire_attribute tiles')

      allocate (tiles(nTiles))
      nc_status = nf90_get_att(ncid, NF90_GLOBAL, 'tiles', tiles)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_att tiles')

      return
   end function nc_get_tiles

   subroutine get_nc_var_scalar_real(ncid, var_name, vals)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      real(kind=wp), intent(out) :: vals

      integer :: nc_status
      integer :: var_id
      real(kind=c_double) :: nc_vals

      nc_status = nf90_inq_varid(ncid, var_name, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inq_varid '//var_name)

      nc_status = nf90_get_var(ncid, var_id, nc_vals)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_var '//var_name)
      vals = real(nc_vals, kind=wp)
      return

   end subroutine get_nc_var_scalar_real

   subroutine get_nc_var_scalar_int(ncid, var_name, vals)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      integer, intent(out) :: vals

      integer :: nc_status
      integer :: var_id

      nc_status = nf90_inq_varid(ncid, var_name, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inq_varid '//var_name)

      nc_status = nf90_get_var(ncid, var_id, vals)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_var '//var_name)
      return

   end subroutine get_nc_var_scalar_int


   subroutine get_nc_var_vec_real(ncid, var_name, start, count, vals)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      integer, intent(in) :: start
      integer, intent(in) :: count
      real(kind=wp), dimension(:), intent(out) :: vals(count)

      integer :: nc_status
      integer :: var_id
      real(kind=c_double), dimension(:) :: nc_vals(count)

      nc_status = nf90_inq_varid(ncid, var_name, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inq_varid '//var_name)

      nc_status = nf90_get_var(ncid, var_id, nc_vals(:), start=[start], count=[count])
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_var '//var_name)

      vals = real(nc_vals, kind=wp)
      return

   end subroutine get_nc_var_vec_real

   subroutine get_nc_var_vec_int(ncid, var_name, start, count, vals)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      integer, intent(in) :: start
      integer, intent(in) :: count
      integer, dimension(:), intent(out) :: vals(count)

      integer :: nc_status
      integer :: var_id

      nc_status = nf90_inq_varid(ncid, var_name, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inq_varid '//var_name)

      nc_status = nf90_get_var(ncid, var_id, vals(:), start=[start], count=[count])
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_var '//var_name)
      return

   end subroutine get_nc_var_vec_int

   subroutine get_nc_var_arr_real(ncid, var_name, start, count, vals)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      integer, dimension(:), intent(in) :: start(2)
      integer, dimension(:), intent(in) :: count(2)
      real(kind=wp), dimension(:, :), intent(out) :: vals(count(1), count(2))

      integer :: nc_status
      integer :: var_id
      real(kind=c_double), dimension(:, :) :: nc_vals(count(1), count(2))

      nc_status = nf90_inq_varid(ncid, var_name, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inq_varid '//var_name)

      nc_status = nf90_get_var(ncid, var_id, nc_vals(:, :), start=start, count=count)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_var '//var_name)

      vals = real(nc_vals, kind=wp)
      return

   end subroutine get_nc_var_arr_real

   subroutine get_nc_var_arr_int(ncid, var_name, start, count, vals)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      integer, dimension(:), intent(in) :: start(2)
      integer, dimension(:), intent(in) :: count(2)
      integer, dimension(:, :), intent(out) :: vals(count(1), count(2))

      integer :: nc_status
      integer :: var_id

      nc_status = nf90_inq_varid(ncid, var_name, var_id)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inq_varid '//var_name)

      nc_status = nf90_get_var(ncid, var_id, vals(:, :), start=start, count=count)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_var '//var_name)
      return

   end subroutine get_nc_var_arr_int


   subroutine get_nc_att_int(ncid, att_name, val)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: att_name
      integer, intent(out) :: val

      integer :: nc_status

      nc_status = nf90_get_att(ncid, NF90_GLOBAL, att_name, val)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_att '//att_name)
      return
   end subroutine get_nc_att_int

   subroutine get_nc_att_real(ncid, att_name, val)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: att_name
      real(kind=wp), intent(out) :: val

      integer :: nc_status
      real(kind=c_double) :: nc_val

      nc_status = nf90_get_att(ncid, NF90_GLOBAL, att_name, nc_val)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_att '//att_name)
      val = real(nc_val, kind=wp)
      return
   end subroutine get_nc_att_real

   subroutine get_nc_att_int_vec(ncid, att_name, val)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: att_name
      integer, dimension(:), allocatable, intent(out) :: val
      
      integer :: n
      integer :: nc_status

      nc_status = nf90_inquire_attribute(ncid, NF90_GLOBAL, att_name, len=n)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inquire_att '//att_name)

      allocate(val(n))

      nc_status = nf90_get_att(ncid, NF90_GLOBAL, att_name, val)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_att '//att_name)
      return
   end subroutine get_nc_att_int_vec

   subroutine get_nc_att_real_vec(ncid, att_name, val)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: att_name
      real(kind=wp), dimension(:), allocatable, intent(out) :: val
      
      integer :: n
      integer :: nc_status
      real(kind=c_double), dimension(:), allocatable :: nc_val

      nc_status = nf90_inquire_attribute(ncid, NF90_GLOBAL, att_name, len=n)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_inquire_att '//att_name)

      allocate(nc_val(n), val(n))

      nc_status = nf90_get_att(ncid, NF90_GLOBAL, att_name, nc_val)
      if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_att '//att_name)

      val = real(nc_val, kind=wp)
      return
   end subroutine get_nc_att_real_vec

   subroutine get_nc_att_str(ncid, att_name, val)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: att_name
    character(len=:), allocatable, intent(out) :: val

    integer :: nc_status
    integer :: length

    nc_status = nf90_inquire_attribute(ncid, NF90_GLOBAL, att_name, len=length)
    if (nc_status /= NF90_NOERR) then
        val = "Unknown"
        return
    end if

    allocate(character(len=length) :: val)

    nc_status = nf90_get_att(ncid, NF90_GLOBAL, att_name, val)
    if (nc_status /= NF90_NOERR) call handle_err(nc_status, ' nf90_get_att '//att_name)
    return
 end subroutine get_nc_att_str

   subroutine handle_err(status, extra)
      integer, intent(in) :: status
      character(len=*), intent(in), optional :: extra

      type(varString) :: msg

      if (status .ne. nf90_noerr) then
         msg = varString(trim(nf90_strerror(status)))
         if (present(extra)) then
            msg = msg + ' ' + extra
         end if
         call FatalErrorMessage(msg%s)
      end if
   end subroutine handle_err


end module netcdf_utils_module
#endif
