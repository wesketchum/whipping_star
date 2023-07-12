project                     (sbn-diyhighfive)
cmake_minimum_required      (VERSION 3.9)

option                      (diy_thread     "Enable diy threading"                 OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-march=native -Wall -O2 -DNDEBUG")# -fPermissive")
#set (CMAKE_CXX_FLAGS_DEBUG " -fno-omit-frame-pointer -fsanitize=address ${CMAKE_CXX_FLAGS_DEBUG}")
#set (CMAKE_LINKER_FLAGS_DEBUG "-fno-omit-frame-pointer -fsanitize=address ${CMAKE_LINKER_FLAGS_DEBUG})")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

if( supported )
    message(STATUS "IPO / LTO enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()


## AddressSanitize
set(CMAKE_C_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C++ compiler during AddressSanitizer builds."
    FORCE)

 #LeakSanitizer
set(CMAKE_C_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during LeakSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C++ compiler during LeakSanitizer builds."
    FORCE)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)



set(CMAKE_CXX_STANDARD 17)

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
    set                     (libraries ${libraries} ${MPI_C_LIBRARIES})# ${MPI_CXX_LIBRARIES})

#find_package                (HDF5)
set                     (libraries ${libraries} ${HDF5_C_LIBRARIES})
#set                     (libraries ${libraries} ${HDF5_C_LIBRARIES} ${NOMAD_INCLUDE_DIRS}/../lib/libnomad.so)

find_package(GSL REQUIRED)
set                     (libraries ${libraries} ${GSL_LIBRARIES})


add_definitions         (-DDIY_NO_THREADS)
# DIY
find_path                   (DIY_INCLUDE_DIRS diy/types.h)
#find_path                   (NOMAD_INCLUDE_DIRS nomad.hpp)
#find_path                   (SGTE_INCLUDE_DIRS sgtelib.hpp)

set (MY_PUBLIC_HEADERS 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GenericBlock.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/opts.h
    )

# HighFive
if                          (HighFive_INCLUDE_DIR)
    include_directories     (SYSTEM ${HighFive_INCLUDE_DIR})
    message                 ("HighFive_INCLUDE_DIR = " ${HighFive_INCLUDE_DIR})
endif                   ()

# MFA
if                          (MFA_INCLUDE_DIR)
    include_directories     (SYSTEM ${MFA_INCLUDE_DIR})
    message                 ("MFA_INCLUDE_DIR = " ${MFA_INCLUDE_DIR})
    add_definitions         (-DMFA_NO_WEIGHTS)
    add_definitions         (-DMFA_SERIAL)
    message                 ("Weights are disabled for MFA")
endif  ()

# LBFGSB
if                          (LBFGSB_INCLUDE_DIR)
    include_directories     (SYSTEM ${LBFGSB_INCLUDE_DIR})
    message                 ("LBFGSB_INCLUDE_DIR = " ${LBFGSB_INCLUDE_DIR})
    set                     (libraries ${libraries} ${LBFGSB_LIBRARIES})
endif  ()

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${DIY_INCLUDE_DIRS}
                            ${GSL_INCLUDE_DIRS}
                            ${HDF5_INCLUDE_DIRS}
			    ${EIGEN3_INCLUDE_DIRS}
			    ${LBFGS_INCLUDE_DIRS}
                            SYSTEM 
			    ${MPI_INCLUDE_PATH}
                            )
                            #${SGTE_INCLUDE_DIRS}
                            #${NOMAD_INCLUDE_DIRS}
                        

set (ROOT_LIBS "-I${ROOTSYS}/include -L${ROOTSYS}/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_LIBS}")
message("CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# Figure out correct library suffix
if                          (APPLE)
  set                       (CMAKE_SHARED_LIBRARY_SUFFIX ".dylib")
  message                   ("Shared library suffix set to dylib")
else                        ()
  set                       (CMAKE_SHARED_LIBRARY_SUFFIX ".so")
  message                   ("Shared library suffix set to so")
endif                       ()

find_file(HAVE_H5_H
    NAMES H5File.hpp
    HINTS "${HighFive_INCLUDE_DIR}/highfive" 
    )
if(NOT HAVE_H5_H)
    message(FATAL_ERROR "Specified HighFive dir not found --- check -DHighFive_INCLUDE_DIR")
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
