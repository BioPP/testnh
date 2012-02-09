# CMake script for TestNH package
# Author: Julien Dutheil
# Created: 07/04/2010

ADD_EXECUTABLE(testnh TestNH.cpp)
TARGET_LINK_LIBRARIES(testnh ${LIBS})

#ADD_EXECUTABLE(stepnh StepNH.cpp)
#TARGET_LINK_LIBRARIES(stepnh ${LIBS})

ADD_EXECUTABLE(mapnh MultinomialClustering.cpp MapNH.cpp)
TARGET_LINK_LIBRARIES(mapnh ${LIBS})

ADD_EXECUTABLE(partnh PartNH.cpp)
TARGET_LINK_LIBRARIES(partnh ${LIBS})

ADD_EXECUTABLE(randnh RandNH.cpp)
TARGET_LINK_LIBRARIES(randnh ${LIBS})

# Install progs
INSTALL(TARGETS
  testnh
  #  stepnh
  mapnh
  partnh
  randnh
  DESTINATION bin)
