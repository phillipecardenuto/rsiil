cmake_minimum_required(VERSION 2.6)
project(Vole)

set(Vole_CMFD ON)
set(Vole_CMFD_Ground_Truth ON)
set(Vole_Shell ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


include(VoleHelperMacros)
include(VoleModuleMacros)

if(WIN32)
	add_definitions(-DNOMINMAX) # prevents min- and max-macro brainf*ck in VS 2003
endif(WIN32)

# Global options
option(VOLE_CMAKE_DEBUG_OUTPUT "Show additonal cmake debug output." OFF)
mark_as_advanced(VOLE_CMAKE_DEBUG_OUTPUT)
option(VOLE_CHECK_DEPENDENCIES "Do additional dependency check with nice error messages." ON)
mark_as_advanced(VOLE_CHECK_DEPENDENCIES)
option(VOLE_SINGLE_TARGETS "Build single targets for each command." OFF)
set(VOLE_MODULE_DIRECTORIES "" CACHE STRING "Semicolon-separated list of module directories for vole, e.g. the forensics or the color code.")
# maybe there is a better place for these two options:
option(Boost_NO_SYSTEM_PATHS "Suppress searching in system paths or other locations outside of BOOST_ROOT." OFF)
set(BOOST_ROOT "" CACHE STRING "The preferred installation prefix for searching for Boost. Set this if the module has problems finding the proper Boost installation.")

# Global variables
set(VOLE_PACKAGE_LIST "" CACHE INTERNAL "Package list" FORCE)
set(VOLE_EXTERNAL_SUBDIRECTORIES "" CACHE INTERNAL "External subdirectories" FORCE)
set(VOLE_MODULE_LIST "" CACHE INTERNAL "Module list" FORCE)

include(VoleFindPackages)
# include(VoleFindExternalModules)

add_subdirectory(common core/common)

include_directories(common shell cmfd ground_truth_db)
add_subdirectory(cmfd cmfd)
add_subdirectory(ground_truth_db ground_truth_db)
add_subdirectory(shell core/shell)

include(VoleCheckDependencies)

INCLUDE("CMakeModules/TargetDoc.cmake" OPTIONAL)

