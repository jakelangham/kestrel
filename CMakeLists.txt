cmake_minimum_required(VERSION 3.20)

project(kestrel LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# There is an error in CMAKE with this flag for pgf90. Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")



####
# Define the actual files and folders that make up the build
###

# Define the executable name
SET(KESTRELEXE kestrel)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(HEADER_FILES 
    ${SRC}/cUTM.h
    ${SRC}/RasterData.h
)

set(KESTREL_src
    ${SRC}/cUTM.cpp
    ${SRC}/RasterData.cpp
    ${SRC}/SetPrecision.f90
    ${SRC}/varStringClass.f90
	${SRC}/Messages.f90
    ${SRC}/utilities.f90
    ${SRC}/Interp2d.f90
	${SRC}/utm.f90
    ${SRC}/RunSettings.f90
    ${SRC}/GeoTiffRead.f90
    ${SRC}/DomainSettings.f90
	${SRC}/InitConds.f90
    ${SRC}/Closures.f90
    ${SRC}/TopogFuncs.f90
    ${SRC}/Parameters.f90
    ${SRC}/Limiters.f90
	${SRC}/SolverSettings.f90
    ${SRC}/OutputSettings.f90
    ${SRC}/TopogSettings.f90
    ${SRC}/Input.f90
    ${SRC}/Grid.f90
    ${SRC}/Equations.f90
    ${SRC}/MorphodynamicRHS.f90
    ${SRC}/HydraulicRHS.f90
    ${SRC}/Redistribute.f90
    ${SRC}/dem.f90
    ${SRC}/UpdateTiles.f90
    ${SRC}/SetSources.f90
    ${SRC}/NetCDFUtils.f90
    ${SRC}/Output.f90
    ${SRC}/Restart.f90
    ${SRC}/TimeStepper.f90
    ${SRC}/cversion.cpp
    ${SRC}/version.f90
    ${SRC}/main.f90
)

# Define the executable in terms of the source files

add_executable(kestrel ${KESTREL_src})

target_include_directories(kestrel
  PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SRC}>
)

# (optional but helpful) put all generated .mod/.smod in one place:
set_target_properties(kestrel PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/all
)

option(WITH_NETCDF "Include NetCDF" ON)

option(WITH_OPENMP "Include OpenMP" ON)

option(DEBUG "Build with -g and the compiler's run-time back-trace flags" OFF)

INCLUDE(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("filesystem" FILESYSTEM_FOUND)
if (NOT FILESYSTEM_FOUND)
    MESSAGE(STATUS "C++ does not support filesystem -- looking for boost")
    INCLUDE(${CMAKE_MODULE_PATH}/SetupBOOST.cmake)
    if(Boost_FOUND)
        add_compile_definitions("HAVE_BOOST_FILESYSTEM")
        target_link_libraries(${KESTRELEXE} PUBLIC Boost::filesystem)
    else()
        MESSAGE(ERROR FATAL "Boost not found")
    endif()
else()
    MESSAGE(STATUS "C++ supports filesystem")
    add_compile_definitions("HAVE_FILESYSTEM")
    target_link_libraries(${KESTRELEXE} PUBLIC stdc++fs)
endif()


INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# Setup libraries
if(WITH_NETCDF)
    MESSAGE(STATUS ${NETCDF_INCLUDES})
    IF(NOT DEFINED NETCDF_INCLUDES)
        INCLUDE(${CMAKE_MODULE_PATH}/SetupNETCDF.cmake)
        IF (NOT NETCDF_FOUND)
            MESSAGE(ERROR FATAL " Could not find NetCDF.  Try setting -DNETCDF_INCLUDES and -DNETCDF_LIBRARIES")
        ENDIF()
    ELSE()
        set(NETCDF_FOUND TRUE)
        # target_link_libraries (kestrel PUBLIC "${NETCDF_LIBRARIES}")
    ENDIF()
    # include_directories(${NETCDF_INCLUDES})
    target_include_directories(kestrel PRIVATE "${NETCDF_INCLUDES}")
    target_link_libraries (kestrel PRIVATE "${NETCDF_LIBRARIES}")
    add_definitions("-DHAVE_NETCDF4")
endif()

IF(NOT DEFINED PROJ_LIBRARIES OR NOT DEFINED PROJ_INCLUDE_DIRS)
    INCLUDE(${CMAKE_MODULE_PATH}/SetupPROJ.cmake)
    target_include_directories(kestrel PRIVATE "${PROJ_INCLUDE_DIRS}")
    target_link_libraries(kestrel PRIVATE "${PROJ_LIBRARIES}")
ELSE()
    MESSAGE(STATUS "USING PROJ FROM ${PROJ_LIBRARIES}")
    target_include_directories(kestrel PRIVATE "${PROJ_INCLUDE_DIRS}")
    target_link_libraries(kestrel PRIVATE "${PROJ_LIBRARIES}")
ENDIF()

if (NOT DEFINED GDAL_INCLUDE_DIRS)
    INCLUDE(${CMAKE_MODULE_PATH}/SetupGDAL.cmake)
    if (TARGET GDAL::GDAL)
        target_link_libraries(kestrel PRIVATE GDAL::GDAL)
    elseif(TARGET PkgConfig::GDAL)
        target_link_libraries(kestrel PRIVATE PkgConfig::GDAL)
    endif()
else()
    target_include_directories(kestrel PRIVATE "${GDAL_INCLUDE_DIRS}")
    target_link_libraries(kestrel PRIVATE "${GDAL_LIBRARIES}")
endif()

if(WITH_OPENMP)
    # Prefer CMake's OpenMP imported targets (modern, portable)
    find_package(OpenMP REQUIRED COMPONENTS Fortran CXX)

    if(OpenMP_Fortran_FOUND)
        message(STATUS "OpenMP Fortran flags: ${OpenMP_Fortran_FLAGS}")
        set(OPENMP_FOUND TRUE)

        # Link the Fortran OpenMP target to the executable
        target_link_libraries(${KESTRELEXE} PRIVATE OpenMP::OpenMP_Fortran)

        # Also link the CXX OpenMP target if available
        if(OpenMP_CXX_FOUND)
            target_link_libraries(${KESTRELEXE} PRIVATE OpenMP::OpenMP_CXX)
        endif()
    else()
        message(FATAL_ERROR "OpenMP not found. Try setting -DOpenMP_Fortran_FLAGS manually or checking your compiler support.")
    endif()
endif()

if(DEBUG)
  # force Debug if the user hasn’t chosen another type
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
  endif()

  # -------- GNU gfortran --------------------------------------
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-g -O0 -fbacktrace -fcheck=all -Wall)
    add_link_options    (-g)

  # -------- Intel oneAPI ifx / ifort --------------------------
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options(-g -O0 -traceback -check all -warn all)
    add_link_options    (-traceback)

  # -------- NVIDIA HPC (nvfortran) ----------------------------
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    add_compile_options(-g -O0 -traceback -Mbounds -Mchkptr)
    add_link_options    (-traceback)
  endif()
endif()

# Set up testing
set(KESTREL_EXECUTABLE "${CMAKE_BINARY_DIR}/bin/kestrel")
set(KESTREL_SERIAL_EXECUTABLE "${CMAKE_BINARY_DIR}/bin/kestrel_serial")
set(KESTREL_PARAL_EXECUTABLE "${CMAKE_BINARY_DIR}/bin/kestrel_paral")

enable_testing()
add_subdirectory(tests)

add_custom_target(purge_mods
  COMMAND ${CMAKE_COMMAND} -E rm -f
          ${CMAKE_SOURCE_DIR}/src/*.mod
          ${CMAKE_SOURCE_DIR}/src/*.smod
          ${CMAKE_SOURCE_DIR}/src/*/*.mod
          ${CMAKE_SOURCE_DIR}/src/*/*.smod
  COMMENT "Removing stray Fortran .mod/.smod files from the source tree"
)

add_custom_target(distclean
  COMMAND ${CMAKE_COMMAND} -E echo "Running distclean..."
  # purge stray Fortran modules (same commands as purge_mods)
  COMMAND ${CMAKE_COMMAND} -E rm -f
          ${CMAKE_SOURCE_DIR}/src/*.mod
          ${CMAKE_SOURCE_DIR}/src/*.smod
          ${CMAKE_SOURCE_DIR}/src/*/*.mod
          ${CMAKE_SOURCE_DIR}/src/*/*.smod
  # remove build artifacts
  COMMAND ${CMAKE_COMMAND} -E rm -f
          ${CMAKE_BINARY_DIR}/CMakeCache.txt
  COMMAND ${CMAKE_COMMAND} -E rm -rf
          ${CMAKE_BINARY_DIR}/CMakeFiles
          ${CMAKE_BINARY_DIR}/cmake_install.cmake
          ${CMAKE_BINARY_DIR}/CMakeScripts
          ${CMAKE_BINARY_DIR}/Testing
          ${CMAKE_BINARY_DIR}/modules
          ${CMAKE_BINARY_DIR}/_deps
  COMMENT "Full distclean: removing build artifacts and stale Fortran modules"
)

