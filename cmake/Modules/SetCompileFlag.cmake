#############################################################################
# Given a list of flags, this function will try each, one at a time,
# and choose the first flag that works.  If no flags work, then nothing
# will be set, unless the REQUIRED key is given, in which case an error
# will be given.
# 
# Call is:
# SET_COMPILE_FLAG(FLAGVAR FLAGVAL (Fortran|C|CXX) <REQUIRED> flag1 flag2...)
# 
# For example, if you have the flag CMAKE_C_FLAGS and you want to add
# warnings and want to fail if this is not possible, you might call this
# function in this manner:
# SET_COMPILE_FLAGS(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" C REQUIRED
#                   "-Wall"     # GNU
#                   "-warn all" # Intel
#                  )
# The optin "-Wall" will be checked first, and if it works, will be
# appended to the CMAKE_C_FLAGS variable.  If it doesn't work, then
# "-warn all" will be tried.  If this doesn't work then checking will
# terminate because REQUIRED was given.  
#
# The reasong that the variable must be given twice (first as the name then
# as the value in quotes) is because of the way CMAKE handles the passing
# of variables in functions; it is difficult to extract a variable's
# contents and assign new values to it from within a function.
#############################################################################

INCLUDE(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXCompilerFlag.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckFortranCompilerFlag.cmake)

FUNCTION(SET_COMPILE_FLAG FLAGVAR FLAGVAL LANG)

    # Make a variable holding the flags. Filter out REQUIRED if it is there.
    SET(FLAG_REQUIRED FALSE)
    SET(FLAG_FOUND FALSE)
    UNSET(FLAGLIST)
    FOREACH (var ${ARGN})
        STRING(TOUPPER "${var}" UP)
        IF(UP STREQUAL "REQUIRED")
            SET(FLAG_REQUIRED TRUE)
        ELSE()
            SET(FLAGLIST ${FLAGLIST} "${var}")
        ENDIF(UP STREQUAL "REQUIRED")
    ENDFOREACH (var ${ARGN})

    # Now, loop over each flag. Some compilers accept compound options such as
    # "-O2 -g -Wall" as a single command-line item, but CMake's flag checks need
    # to evaluate each option individually before accepting the whole group.
    FOREACH(flag ${FLAGLIST})

        SET(FLAG_WORKS FALSE)
        IF(LANG STREQUAL "Fortran")
            SET(FLAG_PARTS "${flag}")
            STRING(REPLACE " " ";" FLAG_PARTS "${flag}")
            SET(FLAG_PARTS_OK TRUE)
            FOREACH(part IN LISTS FLAG_PARTS)
                IF("${part}" STREQUAL "")
                    CONTINUE()
                ENDIF()
                SET(PART_WORKS FALSE)
                CHECK_Fortran_COMPILER_FLAG("${part}" PART_WORKS)
                IF(NOT PART_WORKS)
                    SET(FLAG_PARTS_OK FALSE)
                    BREAK()
                ENDIF()
            ENDFOREACH()
            IF(FLAG_PARTS_OK)
                SET(FLAG_WORKS TRUE)
            ENDIF()
        ELSEIF(LANG STREQUAL "C")
            CHECK_C_COMPILER_FLAG("${flag}" FLAG_WORKS)
        ELSEIF(LANG STREQUAL "CXX")
            CHECK_CXX_COMPILER_FLAG("${flag}" FLAG_WORKS)
        ELSE()
            MESSAGE(FATAL_ERROR "Unknown language in SET_COMPILE_FLAGS: ${LANG}")
        ENDIF()

        IF(FLAG_WORKS)
            SET(${FLAGVAR} "${FLAGVAL} ${flag}" CACHE STRING
                 "Set the ${FLAGVAR} flags" FORCE)
            SET(FLAG_FOUND TRUE)
            BREAK()
        ENDIF()

    ENDFOREACH(flag ${FLAGLIST})

    IF(FLAG_REQUIRED AND NOT FLAG_FOUND)
        MESSAGE(FATAL_ERROR "No compile flags were found")
    ENDIF()

ENDFUNCTION()