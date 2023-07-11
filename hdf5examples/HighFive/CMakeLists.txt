project                     (sbn-diyhighfive)
cmake_minimum_required      (VERSION 2.8)

option                      (diy_thread     "Enable diy threading"                 OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -ggdb")# -fPermissive")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 11)

# OSX flags
if                          (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions           (-DMAC_OSX)
  set                       (CMAKE_MACOSX_RPATH on)

    # --- following RPATH settings are for Sierra w/ Clang, hopefully they don't hurt other versions
    # ref: https://cmake.org/Wiki/CMake_RPATH_handling
    # use, i.e. don't skip, the full RPATH for the build tree
    set                     (CMAKE_SKIP_BUILD_RPATH            false)
    # when building, don't use the install RPATH already (but later on when installing)
    set                     (CMAKE_BUILD_WITH_INSTALL_RPATH    false)
    # set RPATH to install path
    set                     (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set                     (CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
    # the RPATH to be used when installing, but only if it's not a system directory
    list                    (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
                                                               "${CMAKE_INSTALL_PREFIX}/lib"
                                                               isSystemDir)
    if                      ("${isSystemDir}" STREQUAL         "-1")
      set                   (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    endif                   ()
endif                       (${CMAKE_SYSTEM_NAME} MATCHES      "Darwin")

 # MPI
find_package                (MPI REQUIRED)
    set                     (libraries ${libraries} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})

#find_package                (HDF5)
set                     (libraries ${libraries} ${HDF5_C_LIBRARIES})


add_definitions         (-DDIY_NO_THREADS)
# DIY
find_path                   (DIY_INCLUDE_DIRS diy/types.h)

set (MY_PUBLIC_HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GenericBlock.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/opts.h
    )

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${DIY_INCLUDE_DIRS}
                            ${HDF5_INCLUDE_DIRS}
                            ${HIGHFIVE_DIR}/include 
                            SYSTEM ${MPI_INCLUDE_PATH}
                            )

# Figure out correct library suffix
if                          (APPLE)
  set                       (CMAKE_SHARED_LIBRARY_SUFFIX ".dylib")
  message                   ("Shared library suffix set to dylib")
else                        ()
  set                       (CMAKE_SHARED_LIBRARY_SUFFIX ".so")
  message                   ("Shared library suffix set to so")
endif                       ()


# Test if user has supplied -DHIGHFIVE_DIR                        
if (NOT HIGHFIVE_DIR)
    message                 (FATAL_ERROR " Need HighFive in order to build program --- check -DHIGHFIVE_DIR")
endif                       ()





find_file(HAVE_H5_H
    NAMES H5File.hpp
    HINTS "${HIGHFIVE_DIR}/include/highfive" 
    )
if(NOT HAVE_H5_H)
    message(FATAL_ERROR "Specified HighFive dir not found --- check -DHIGHFIVE_DIR")
endif()

# If all found, add build target
add_subdirectory        (src)

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES
    "Docker;benchmark;test;.git;~;build;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
