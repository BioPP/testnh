# CMake script for TestNH package
# Author: Julien Dutheil
# Created: 07/04/2010

ADD_EXECUTABLE(testnh TestNH.cpp)
TARGET_LINK_LIBRARIES(testnh ${LIBS})

ADD_EXECUTABLE(stepnh StepNH.cpp)
TARGET_LINK_LIBRARIES(stepnh ${LIBS})

# Install progs
INSTALL(TARGETS
  testnh
  stepnh
  DESTINATION bin)
