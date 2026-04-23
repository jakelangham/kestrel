# Collect the drag sources and export them to the parent
# Do not create a lib here
set(EROSION_TRANSITION_SOURCES
  ${CMAKE_SOURCE_DIR}/src/ErosionTransitionInterface.f90
  ${CMAKE_SOURCE_DIR}/src/ErosionTransitionTypes.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/ErosionTransition_NoTransition.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/ErosionTransition_SmoothTransition.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/ErosionTransition_StepTransition.f90
  PARENT_SCOPE
)