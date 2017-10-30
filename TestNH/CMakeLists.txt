# CMake script for TestNH package
# Authors:
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 07/04/2010

# Executables of TestNH
add_executable (testnh TestNH.cpp)
#add_executable (stepnh StepNH.cpp)
add_executable (mapnh MultinomialClustering.cpp MapNH.cpp)
add_executable (partnh PartNH.cpp)
add_executable (randnh RandNH.cpp)

set (testnh-targets
  testnh
  # stepnh
  mapnh
  partnh
  randnh
  )

foreach (target ${testnh-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC})
    set_target_properties (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED})
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${testnh-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
