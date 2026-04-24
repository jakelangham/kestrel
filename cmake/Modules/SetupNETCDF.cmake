IF(NETCDF_FOUND)
    return()
ENDIF()

set (NETCDF_F90 "YES")
find_package(NetCDF)
    