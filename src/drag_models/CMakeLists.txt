# Collect the drag sources and export them to the parent
# Do not create a lib here
set(DRAG_SOURCES
  ${CMAKE_SOURCE_DIR}/src/DragInterface.f90
  ${CMAKE_SOURCE_DIR}/src/DragTypes.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Drag_Chezy.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Drag_Manning.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Drag_Variable.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Drag_Coulomb.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Drag_Pouliquen.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Drag_Voellmy.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Drag_Edwards2019.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Drag_Powerlaw.f90
  PARENT_SCOPE
)