
IF(GDAL_FOUND)
    RETURN()
ENDIF()

find_package(GDAL CONFIG)

IF(NOT GDAL_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GDAL IMPORTED_TARGET gdal)
ENDIF()

IF(NOT GDAL_FOUND)
    MESSAGE(FATAL ERROR " Could not locate required GDAL library.")
ELSE()
    MESSAGE(STATUS " Found GDAL library.")
ENDIF()