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


! This module defines an interface to cUTM.cpp that contains a class proj_transformer
! that uses the proj library to convert between WGS84 and WGS84 / UTM coordinates
!
! See cUTM.h and cUTM.cpp
module utm_module

    use, intrinsic :: iso_c_binding
    use set_precision_module, only: wp
    use utilities_module, only : pair
    
    implicit none

    private
    public :: LatLonToZoneNumber, ZoneNumberToCentralLongitude, LatLonToUtmEpsg
    public :: new, delete, proj_transformer
   
    ! Fortran interface to proj_transformer
    ! Use a Fortran type to represent a C++ class in an opaque manner
    type proj_transformer
        private
        type(c_ptr) :: object = c_null_ptr ! Pointer to C++ proj_transformer object
    contains
        final :: proj_transformer__delete ! Destructor in Fortran

        procedure :: proj_transformer_wgs84_to_utm_f ! Transform from wgs84 to UTM with coordinates passed as reals
        procedure :: proj_transformer_wgs84_to_utm_pair ! Transform from wgs84 to UTM with coordinates passed as type(pair)
        generic, public :: wgs84_to_utm => proj_transformer_wgs84_to_utm_f, proj_transformer_wgs84_to_utm_pair ! Overloaded interface to wgs84_to_utm
        procedure :: proj_transformer_utm_to_wgs84_f ! Transform from UTM to WGS84 with coordinates passed as reals
        procedure :: proj_transformer_utm_to_wgs84_pair ! Transform from UTM to WGS84 with coordinates passed as type(pair)
        generic, public :: utm_to_wgs84 => proj_transformer_utm_to_wgs84_f, proj_transformer_utm_to_wgs84_pair ! Overloaded interface to utm_to_wgs84
    end type proj_transformer
    interface
        ! Bindings to C functions
        function proj_transformer__new_c(utm_code) result(this) bind(C, name="proj_transformer__new")
            import
            implicit none
            type(c_ptr) :: this
            integer(kind=c_int), value :: utm_code
        end function proj_transformer__new_c

        subroutine proj_transformer__delete_c(self) bind(C,name="proj_transformer__delete")
            import
            implicit none
            type(c_ptr), value :: self
        end subroutine proj_transformer__delete_c

        function proj_transformer__wgs84_to_utm_c(this, latitude, longitude) result(en) bind(C,name="proj_transformer__wgs84_to_utm")
            import
            implicit none
            type(c_ptr), intent(in), value :: this
            real(kind=c_double), intent(in), value :: latitude
            real(kind=c_double), intent(in), value :: longitude
            type(c_ptr) :: en
        end function proj_transformer__wgs84_to_utm_c

        function proj_transformer__utm_to_wgs84_c(this, easting, northing) result(latlon) bind(C,name="proj_transformer__utm_to_wgs84")
            import
            type(c_ptr), intent(in), value :: this
            real(kind=c_double), intent(in), value :: easting
            real(kind=c_double), intent(in), value :: northing
            type(c_ptr) :: latlon
        end function proj_transformer__utm_to_wgs84_c

        function LatLonToZoneNumber(latitude, longitude) result(zone_number) bind(C,name="latlon_to_zone_number")
            import
            real(kind=c_double), intent(in), value :: latitude
            real(kind=c_double), intent(in), value :: longitude
            integer(kind=c_int) :: zone_number
        end function LatLonToZoneNumber

        function ZoneNumberToCentralLongitude(zone_number) result(lon) bind(C, name="zone_number_to_central_longitude")
            import
            integer(kind=c_int), intent(in), value :: zone_number
            integer(kind=c_int) :: lon
        end function ZoneNumberToCentralLongitude

        function LatLonToUtmEpsg(latitude, longitude) result(utm_code) bind(C, name="latlon_to_utm_epsg")
            import
            real(kind=c_double), intent(in), value :: latitude
            real(kind=c_double), intent(in), value :: longitude
            integer(kind=c_int) :: utm_code
        end function LatLonToUtmEpsg
    end interface

    interface new
        module procedure proj_transformer__new
    end interface new
    interface delete
        module procedure proj_transformer__delete
    end interface delete

contains
    subroutine proj_transformer__new(this, utm_code)
        implicit none
        type(proj_transformer), intent(out) :: this
        integer(kind=c_int), intent(in) :: utm_code
        this%object = proj_transformer__new_c(utm_code)
    end subroutine proj_transformer__new

    subroutine proj_transformer__delete(this)
        implicit none
        type(proj_transformer), intent(inout) :: this
        call proj_transformer__delete_c(this%object)
        this%object = c_null_ptr
    end subroutine proj_transformer__delete

    function proj_transformer_wgs84_to_utm_f(this, latitude, longitude) result(EN)
        implicit none
        class(proj_transformer), intent(in) :: this
        real(kind=c_double), intent(in) :: latitude
        real(kind=c_double), intent(in) :: longitude
        type(pair) :: EN

        type(c_ptr) :: EN_ptr
        real(kind=c_double), pointer :: a(:)

        EN_ptr = proj_transformer__wgs84_to_utm_c(this%object, latitude, longitude)

        call c_f_pointer(EN_ptr,a,[2])

        EN%first = a(1)
        EN%second = a(2)

    end function proj_transformer_wgs84_to_utm_f

    function proj_transformer_wgs84_to_utm_pair(this, latlon) result(EN)
        implicit none
        class(proj_transformer), intent(in) :: this
        type(pair), intent(in) :: latlon
        type(pair) :: EN

        type(c_ptr) :: EN_ptr
        real(kind=c_double), pointer :: a(:)

        EN_ptr = proj_transformer__wgs84_to_utm_c(this%object, latlon%first, latlon%second)

        call c_f_pointer(EN_ptr,a,[2])

        EN%first = a(1)
        EN%second = a(2)

    end function proj_transformer_wgs84_to_utm_pair

    function proj_transformer_utm_to_wgs84_f(this, easting, northing) result(latlon)
        implicit none
        class(proj_transformer), intent(in) :: this
        real(kind=c_double), intent(in) :: easting
        real(kind=c_double), intent(in) :: northing
        type(pair) :: latlon

        type(c_ptr) :: latlon_ptr
        real(kind=c_double), pointer :: a(:)

        latlon_ptr = proj_transformer__utm_to_wgs84_c(this%object, easting, northing)

        call c_f_pointer(latlon_ptr,a,[2])

        latlon%first = a(1)
        latlon%second = a(2)

    end function proj_transformer_utm_to_wgs84_f

    function proj_transformer_utm_to_wgs84_pair(this, EN) result(latlon)
        implicit none
        class(proj_transformer), intent(in) :: this
        type(pair), intent(in) :: EN
        type(pair) :: latlon

        type(c_ptr) :: latlon_ptr
        real(kind=c_double), pointer :: a(:)

        latlon_ptr = proj_transformer__utm_to_wgs84_c(this%object, EN%first, EN%second)

        call c_f_pointer(latlon_ptr,a,[2])

        latlon%first = a(1)
        latlon%second = a(2)

    end function proj_transformer_utm_to_wgs84_pair
    
end module utm_module
