IF(NOT Boost_FOUND)
    find_package(Boost 1.46.9 REQUIRED COMPONENTS "filesystem")
ENDIF()