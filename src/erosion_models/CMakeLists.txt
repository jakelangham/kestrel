# Collect the erosion sources and export them to the parent
# Do not create a lib here
set(EROSION_SOURCES
  ${CMAKE_SOURCE_DIR}/src/ErosionInterface.f90
  ${CMAKE_SOURCE_DIR}/src/ErosionTypes.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Erosion_Fluid.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Erosion_Mixed.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Erosion_Simple.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Erosion_Granular.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Erosion_NoErosion.f90
  PARENT_SCOPE
)