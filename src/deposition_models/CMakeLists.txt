# Collect the deposition sources and export them to the parent
# Do not create a lib here
set(DEPOSITION_SOURCES
  ${CMAKE_SOURCE_DIR}/src/DepositionInterface.f90
  ${CMAKE_SOURCE_DIR}/src/DepositionTypes.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Deposition_NoDeposition.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Deposition_SimpleHinderedSettling.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Deposition_SpearmanManningHinderedSettling.f90
  PARENT_SCOPE
)