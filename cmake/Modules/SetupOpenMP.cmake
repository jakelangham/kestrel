IF(OPENMP_FOUND)
    return()
ENDIF()

set (OPENMP_F90 "YES")
find_package(OpenMP_Fortran)