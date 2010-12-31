# CMake script for TestNH package
# Author: Julien Dutheil
# Created: 07/04/2010

ADD_EXECUTABLE(testnh TestNH.cpp)
TARGET_LINK_LIBRARIES(testnh ${LIBS})

ADD_EXECUTABLE(stepnh StepNH.cpp)
TARGET_LINK_LIBRARIES(stepnh ${LIBS})

ADD_EXECUTABLE(mapnh ChiClustering.cpp MultinomialClustering.cpp MapNH.cpp)
TARGET_LINK_LIBRARIES(mapnh ${LIBS})

ADD_EXECUTABLE(partnh PartNH.cpp)
TARGET_LINK_LIBRARIES(partnh ${LIBS})

# Install progs
INSTALL(TARGETS
  testnh
  stepnh
  mapnh
  partnh
  DESTINATION bin)
