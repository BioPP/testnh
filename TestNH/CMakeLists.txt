# CMake script for TestNH package
# Authors:
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 07/04/2010

# Executables of TestNH
add_library(bppTools bppTools.cpp)

#add_executable (testnh TestNH.cpp)
#add_executable (stepnh StepNH.cpp)
#add_executable (mapnh MultinomialClustering.cpp MapNH.cpp bppTools.cpp)
add_executable (mapnh MapNH.cpp)
# add_executable (partnh PartNH.cpp)
# add_executable (randnh RandNH.cpp)

set (testnh-targets
  # testnh
  # stepnh
  mapnh
  # partnh
  # randnh
  )

if (BUILD_STATIC)
  target_link_libraries (bppTools ${BPP_LIBS_STATIC} Eigen3::Eigen)
  set_target_properties (bppTools PROPERTIES LINK_SEARCH_END_STATIC TRUE)
else (BUILD_STATIC)
  target_link_libraries (bppTools ${BPP_LIBS_SHARED} Eigen3::Eigen)
  set_target_properties (bppTools PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif (BUILD_STATIC)

foreach (target ${testnh-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC} bppTools)
    set_target_properties (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED} bppTools)
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${testnh-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
