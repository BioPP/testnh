# SPDX-FileCopyrightText: The Bio++ Development Group
#
# SPDX-License-Identifier: CECILL-2.1

# Executables of TestNH
add_executable(testnh TestNH.cpp)
add_executable(mapnh MapNH.cpp)
add_executable(clustnh MultinomialClustering.cpp ClustNH.cpp)
#add_executable (stepnh StepNH.cpp)
add_executable(partnh PartNH.cpp)
# add_executable (randnh RandNH.cpp)

set(testnh-targets
    testnh
    clustnh
    mapnh
    # stepnh
    partnh
    # randnh
)

foreach(target ${testnh-targets})
    # Link (static or shared)
    if(BUILD_STATIC)
        target_link_libraries(${target} ${BPP_LIBS_STATIC})
        set_target_properties(${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
        target_link_options(${target} PRIVATE -static-libgcc -static-libstdc++)
    else(BUILD_STATIC)
        target_link_libraries(${target} ${BPP_LIBS_SHARED})
        set_target_properties(
            ${target}
            PROPERTIES
                INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
        )
    endif(BUILD_STATIC)
endforeach(target)

install(TARGETS ${testnh-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
