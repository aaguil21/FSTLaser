# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(FstTracking)
set(CMAKE_INSTALL_PREFIX $ENV{HOME}/WorkSpace/STAR/ForwardSiliconTracker/FstTracking)

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# include IstTracking
# include_directories("${PROJECT_SOURCE_DIR}/src/IstTracking")
# add_subdirectory(src/IstTracking)

# include FstUtil
include_directories("${PROJECT_SOURCE_DIR}/src/FstUtil")
add_subdirectory(src/FstUtil)

# include FstNoiseStudy
include_directories("${PROJECT_SOURCE_DIR}/src/FstNoiseStudy")
add_subdirectory(src/FstNoiseStudy)

# include FstClusterMaker
include_directories("${PROJECT_SOURCE_DIR}/src/FstClusterMaker")
add_subdirectory(src/FstClusterMaker)

# include FstQAStudy
include_directories("${PROJECT_SOURCE_DIR}/src/FstQAStudy")
add_subdirectory(src/FstQAStudy)

# include FstTracking
include_directories("${PROJECT_SOURCE_DIR}/src/FstTracking")
add_subdirectory(src/FstTracking)

# add the executable
# add_executable(FstTracking ${PROJECT_SOURCE_DIR}/src/IstTracking/IstTracking.cxx)
# target_link_libraries(FstTracking)

# add the install targets
# install (TARGETS FstTracking DESTINATION bin)


