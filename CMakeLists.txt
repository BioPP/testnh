# CMake script for TestNH package
# Authors:
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 07/04/2010

cmake_minimum_required (VERSION 2.8.12)
project (testnh CXX)

# Compile options
set (private-compile-options -std=c++11 -Wall -Weffc++ -Wshadow -Wconversion)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NOT DOC_COMPRESS)
  SET(DOC_COMPRESS gzip CACHE STRING
      "Set program for compressing documentation."
      FORCE)
ENDIF(NOT DOC_COMPRESS)

IF(NOT DOC_COMPRESS_EXT)
  SET(DOC_COMPRESS_EXT gz CACHE STRING
      "Set extension of compressed documentation."
      FORCE)
ENDIF(NOT DOC_COMPRESS_EXT)


IF(NO_DEP_CHECK)
  MESSAGE(STATUS "Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)

#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE(STATUS "Static linkage requested.")
  list (APPEND private-compile-options -static -static-libgcc)
ENDIF()

#build man pages?
IF(NOT DEFINED MAN)
  SET(MAN UNIX)
ENDIF(NOT DEFINED MAN)

#find executables for documentation
FIND_PROGRAM(MAKEINFO_EXE NAMES makeinfo)
FIND_PROGRAM(PDFTEX_EXE NAMES pdftex)
FIND_PROGRAM(NROFF_EXE NAMES nroff)

# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

find_package (bpp-phyl 9.1.3 REQUIRED)

IF(MAKEINFO_EXE)
  MESSAGE("-- Found makeinfo here: ${MAKEINFO_EXE}")
  MESSAGE("   Adding targets: info, html")

  IF(INFO)
    SET(ADD_INFO_TO "ALL")
    MESSAGE("   Adding target info to target all")
  ENDIF(INFO)

  ADD_CUSTOM_TARGET(info
    ${ADD_INFO_TO}
    COMMAND ${MAKEINFO_EXE} testnh.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )

  ADD_CUSTOM_TARGET(html
    COMMAND ${MAKEINFO_EXE} --html --css-ref=http://www.w3.org/StyleSheets/Core/Oldstyle testnh.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
ENDIF(MAKEINFO_EXE)

IF(PDFTEX_EXE)
  MESSAGE("-- Found pdftex here: ${PDFTEX_EXE}")
  MESSAGE("   Adding target: pdf")

  ADD_CUSTOM_TARGET(pdf
    COMMAND ${PDFTEX_EXE} testnh.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "doc/testnh.pdf;doc/testnh.aux;doc/testnh.cp;doc/testnh.fn;doc/testnh.info;doc/testnh.ky;doc/testnh.log;doc/testnh.pg;doc/testnh.toc;doc/testnh.tp;doc/testnh.vr")
ENDIF(PDFTEX_EXE)

IF(NROFF_EXE)
  MESSAGE("-- Found nroff here: ${NROFF_EXE}")
  MESSAGE("   Adding targets: man")

  ADD_CUSTOM_TARGET(man
    ALL
    COMMAND cp testnh.1.txt testnh.1
    COMMAND gzip -f testnh.1
    COMMAND cp mapnh.1.txt mapnh.1
    COMMAND gzip -f mapnh.1
    COMMAND cp partnh.1.txt partnh.1
    COMMAND gzip -f partnh.1
    COMMAND cp randnh.1.txt randnh.1
    COMMAND gzip -f randnh.1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
    )
  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "man/testnh.1.gz;man/mapnh.1.gz;man/partnh.1.gz;man/randnh.1.gz")
ENDIF(NROFF_EXE)

# Subdirectories
add_subdirectory (TestNH)
add_subdirectory (doc)
add_subdirectory (man)

ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "testnh")
SET(CPACK_PACKAGE_VENDOR "Bastien Boussau and Julien Dutheil")
SET(CPACK_PACKAGE_VERSION "1.1.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Non-homogeneous model testing")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 ".*\\\\.git"
 ".*\\\\.gz"
 ".*\\\\.zip"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 ".*\\\\.dmg"
 ".*\\\\.sh"
 ".*\\\\..*\\\\.swp"
 ".*stamp"
 "\\\\.sh"
 "examples/.*/.*\\\\BIC.*"
 "examples/.*/.*\\\\AIC.*"
 "examples/.*/.*\\\\.out"
 "examples/.*/.*\\\\.profile.*"
 "examples/.*/.*\\\\.message.*"
 "examples/.*/.*\\\\.ml.*"
 "examples/.*/.*\\\\.counts.*"
 "examples/.*/.*\\\\.cluster.*"
 "examples/.*/.*\\\\.params.*"
 "examples/.*/.*\\\\.infos"
 "examples/.*/.*\\\\.csv"
 "examples/.*/.*\\\\.R"
 "examples/.*/.*\\\\.txt"
 "examples/.*/treefile.*"
 "examples/.*/Figures/.*"
 "examples/.*/rst"
 "examples/.*/rst1"
 "examples/.*/rub"
 "examples/.*/lnf"
 "examples/.*/mlc"
 "examples/.*/4fold.nuc"
 "examples/.*/2NG.*"
 "examples/.*/PAML.*/rst"
 "examples/.*/PAML.*/rst1"
 "examples/.*/PAML.*/rub"
 "doc/testnh/"
 "doc/testnh\\\\.info"
 "doc/testnh\\\\.toc"
 "doc/testnh\\\\.vr"
 "doc/testnh\\\\.tp"
 "doc/testnh\\\\.log"
 "doc/testnh\\\\.fn"
 "doc/testnh\\\\.ky"
 "doc/testnh\\\\.cp"
 "doc/testnh\\\\.pg"
 "doc/testnh\\\\.aux"
 "doc/testnh\\\\.pdf"
 "man/.*\\\\.1.gz"
 "examples/.*/PAML.*/lnf"
 "examples/.*/PAML.*/mlc"
 "examples/.*/PAML.*/4fold.nuc"
 "examples/.*/PAML.*/2NG.*"
 "TestNH/\\\\..*"
 "TestNH/testnh"
 "TestNH/mapnh"
 "TestNH/partnh"
 "TestNH/randnh"
 "debian/tmp"
 "debian/testnh"
 "debian/testnh\\\\.substvars"
 "debian/testnh\\\\.debhelper"
 "debian/debhelper\\\\.log"
 "install_manifest.txt"
 "DartConfiguration.tcl"
 ${CPACK_SOURCE_IGNORE_FILES}
)
IF (MACOS)
  SET(CPACK_GENERATOR "Bundle")
ENDIF()

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
# 'clean' is not (yet) a first class target. However, we need to clean the directories before building the sources:
IF("${CMAKE_GENERATOR}" MATCHES "Make")
  ADD_CUSTOM_TARGET(make_clean
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
  )
  ADD_DEPENDENCIES(dist make_clean)

  ADD_CUSTOM_TARGET(make_clean_man
    COMMAND rm -f *.${DOC_COMPRESS_EXT}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
  )
  ADD_DEPENDENCIES(dist make_clean_man)
ENDIF()

IF(NOT NO_DEP_CHECK)
IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(origdist COMMAND cp ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz ../${CPACK_DEBSOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(origdist dist)
ADD_CUSTOM_TARGET(deb dpkg-buildpackage -uc -us -i${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(deb origdist)
ADD_DEPENDENCIES(deb info man)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info man)

ENDIF(UNIX)

ENDIF(NOT NO_DEP_CHECK)
