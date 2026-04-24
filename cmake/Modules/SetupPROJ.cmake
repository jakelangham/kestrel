
if(PROJ_FOUND)
    return()
endif()

IF(NOT TIFF_FOUND)

    find_package(TIFF)

    # If not found, try pkg-config
    IF(NOT TIFF_FOUND)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(TIFF libtiff-4)
    ELSE()
        MESSAGE(STATUS " FOUND TIFF USING FIND_PACKAGE")
    ENDIF()
ENDIF()

IF(NOT TIFF_FOUND)
    MESSAGE(FATAL ERROR " Could not locate TIFF library, required dependency of proj.")
ELSE()
    MESSAGE(STATUS " Found TIFF library.")
ENDIF(NOT TIFF_FOUND)

# Try find_package for cmake build
find_package(PROJ CONFIG)

# If not found, try pkg-config
IF(NOT PROJ_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(PROJ proj)
ENDIF()
IF(NOT PROJ_FOUND)
    MESSAGE(FATAL ERROR " Could not locate required proj library.")
ELSE()
    MESSAGE(STATUS " Found proj library.")
ENDIF()
