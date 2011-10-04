# CMake script for TestNH
# Author: Julien Dutheil
# Created: 04/10/2011

IF(MAN)
  INSTALL(FILES testnh.1.gz DESTINATION share/man/man1)
  INSTALL(FILES mapnh.1.gz DESTINATION share/man/man1)
  INSTALL(FILES partnh.1.gz DESTINATION share/man/man1)
  INSTALL(FILES randnh.1.gz DESTINATION share/man/man1)
ENDIF(MAN)
