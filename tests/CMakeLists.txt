find_program(JULIA_EXECUTABLE julia REQUIRED)

# path to your julia test script (adjust as needed)
set(JULIA_TEST_SCRIPT ${CMAKE_SOURCE_DIR}/tests/runall.jl)

# add_test(NAME tests
#          # Force color output if Julia supports it; fall back to default behavior
#          COMMAND ${JULIA_EXECUTABLE} 
#             --color=yes 
#             --project=${CMAKE_SOURCE_DIR}/tests/julia 
#             ${JULIA_TEST_SCRIPT}
# )

# add_test(NAME "tests:1d"
#          # Force color output if Julia supports it; fall back to default behavior
#          COMMAND ${JULIA_EXECUTABLE} 
#             --color=yes 
#             --project=${CMAKE_SOURCE_DIR}/tests/julia 
#             ${JULIA_TEST_SCRIPT} 1d
# )

macro(add_test_suite testname test_suite)
    add_test(NAME ${testname} COMMAND ${JULIA_EXECUTABLE} 
        --color=yes 
        --project=${CMAKE_SOURCE_DIR}/tests ${JULIA_TEST_SCRIPT} ${test_suite}
    )
    set_tests_properties(${testname} PROPERTIES LABELS "${test_suite}")
    list(APPEND _all_julia_tests ${testname})
endmacro()

# add_test_suite(tests all)
add_test_suite("tests:1d" 1d)
add_test_suite("tests:2d" 2d)
add_test_suite("tests:noflow" noflow)
add_test_suite("tests:identical" identical)
add_test_suite("tests:netcdf" netcdf)

set_tests_properties(${_all_julia_tests} PROPERTIES
  PASS_REGULAR_EXPRESSION "Results: .*passed"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests
  ENVIRONMENT "TERM=xterm-256color"
)

## Simple `tests` target: build the `kestrel` executable, then run ctest 
# Invoke with: `cmake --build build --target tests-1d` etc

add_custom_target(tests-1d
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -V -R "tests:1d"
    DEPENDS kestrel
)

add_custom_target(tests-2d
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -V -R "tests:2d"
    DEPENDS kestrel
)

add_custom_target(tests-noflow
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -V -R "tests:noflow"
    DEPENDS kestrel
)

add_custom_target(tests-identical
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -V -R "tests:identical"
    DEPENDS kestrel
)

add_custom_target(tests-netcdf
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -V -R "tests:netcdf"
    DEPENDS kestrel
)

# If there is a `tests/julia` subdirectory with extra CMake helpers, include it.
if(EXISTS ${CMAKE_SOURCE_DIR}/tests/julia/CMakeLists.txt)
    add_subdirectory(${CMAKE_SOURCE_DIR}/tests/julia)
endif()