######################################################
# Determine and set the Fortran compiler flags we want 
######################################################

####################################################################
# Make sure that the default build type is RELEASE if not specified.
####################################################################
# Set Fortran compiler flags explicitly by toolchain instead of relying on the
# legacy custom flag-check helper, which fails on modern CMake/GCC setups.

# Make sure the build type is uppercase
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, PROFILE, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are DEBUG, PROFILE, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "PROFILE")
    SET (CMAKE_BUILD_TYPE PROFILE CACHE STRING
      "Choose the type of build, options are DEBUG, PROFILE, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "TESTING")
    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
      "Choose the type of build, options are DEBUG, PROFILE, RELEASE, or TESTING."
      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, PROFILE, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, PROFILE, RELEASE, or TESTING")
ENDIF()

IF(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    SET(CMAKE_Fortran_FLAGS "-ffree-line-length-0 -fno-range-check -cpp -Wall")
    SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fcheck=bounds -Wall")
    SET(CMAKE_Fortran_FLAGS_PROFILE "-O2 -g -pg -Wall -fcheck=all")
    SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -g -fcheck=all -Wall")
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -g -fcheck=all -Wall -Wextra")
ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(CMAKE_Fortran_FLAGS "-heap-arrays -cpp")
    SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -check bounds -warn all")
    SET(CMAKE_Fortran_FLAGS_PROFILE "-O2 -pg")
    SET(CMAKE_Fortran_FLAGS_TESTING "-O2")
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -g -traceback")
ELSEIF(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    SET(CMAKE_Fortran_FLAGS "-Mfreeform -Mbackslash")
    SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -Mbounds")
    SET(CMAKE_Fortran_FLAGS_PROFILE "-O2 -pg")
    SET(CMAKE_Fortran_FLAGS_TESTING "-O2")
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -g")
ELSE()
    SET(CMAKE_Fortran_FLAGS "")
    SET(CMAKE_Fortran_FLAGS_DEBUG "-g")
    SET(CMAKE_Fortran_FLAGS_PROFILE "-O2")
    SET(CMAKE_Fortran_FLAGS_TESTING "-O2")
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
ENDIF()

# Preserve the CMake cache values for the chosen build type.
SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" CACHE STRING "Fortran release flags" FORCE)
SET(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}" CACHE STRING "Fortran testing flags" FORCE)
SET(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE}" CACHE STRING "Fortran profile flags" FORCE)
SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" CACHE STRING "Fortran debug flags" FORCE)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" CACHE STRING "Fortran general flags" FORCE)

