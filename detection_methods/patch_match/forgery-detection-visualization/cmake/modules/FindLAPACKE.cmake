# - Find the LAPACKE library
#
# Usage:
#   FIND_PACKAGE(LAPACKE [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   LAPACK_FOUND               ... true if LAPACKE is found on the system
#   LAPACK_LIBRARIES           ... full path to LAPACKE library
#   LAPACK_INCLUDES            ... LAPACKE include directory
#

IF(NOT LAPACKE_ROOT AND ENV{LAPACKEDIR})
  SET(LAPACKE_ROOT $ENV{LAPACKEDIR})
ENDIF()

# Check if we can use PkgConfig
FIND_PACKAGE(PkgConfig)

#Determine from PKG
IF(PKG_CONFIG_FOUND AND NOT LAPACKE_ROOT)
  PKG_CHECK_MODULES( PKG_LAPACKE QUIET "lapacke")
ENDIF()

IF(LAPACKE_ROOT)
    #find libs
    FIND_LIBRARY(
        LAPACKE_LIB
        NAMES "lapacke" "LAPACKE" "liblapacke"
        PATHS ${LAPACKE_ROOT}
        PATH_SUFFIXES "lib" "lib64"
        DOC "LAPACKE Library"
        NO_DEFAULT_PATH
        )
    FIND_LIBRARY(
        LAPACK_LIB
        NAMES "lapack" "LAPACK" "liblapack"
        PATHS ${LAPACKE_ROOT}
        PATH_SUFFIXES "lib" "lib64"
        DOC "LAPACK Library"
        NO_DEFAULT_PATH
        )
    FIND_PATH(
        LAPACKE_INCLUDES
        NAMES "lapacke.h"
        PATHS ${LAPACKE_ROOT}
        PATH_SUFFIXES "include"
        DOC "LAPACKE Include Directory"
        NO_DEFAULT_PATH
        )

ELSE()
    FIND_LIBRARY(
        LAPACKE_LIB
        NAMES "lapacke" "liblapacke"
        PATHS
        ${PKG_LAPACKE_LIBRARY_DIRS}
        ${LIB_INSTALL_DIR}
        /usr/lib64
        /usr/lib
        /usr/local/lib64
        /usr/local/lib
        /sw/lib
        /opt/local/lib
	/home/pariasm/local/lib
        DOC "LAPACKE Library"
        )
    FIND_LIBRARY(
       LAPACK_LIB
        NAMES "lapack" "liblapack"
        PATHS
        ${PKG_LAPACKE_LIBRARY_DIRS}
        ${LIB_INSTALL_DIR}
        /usr/lib64
        /usr/lib
        /usr/local/lib64
        /usr/local/lib
        /sw/lib
        /opt/local/lib
	/home/pariasm/local/lib
        DOC "LAPACK Library"
        )
    FIND_PATH(
        LAPACKE_INCLUDES
        NAMES "lapacke.h"
        PATHS
        ${PKG_LAPACKE_INCLUDE_DIRS}
        ${INCLUDE_INSTALL_DIR}
        /usr/include
        /usr/local/include
        /sw/include
        /opt/local/include
	/home/pariasm/local/include
        DOC "LAPACKE Include Directory"
        )
ENDIF(LAPACKE_ROOT)

SET(LAPACK_LIBRARIES ${LAPACKE_LIB} ${LAPACK_LIB})
SET(LAPACK_INCLUDE_DIR ${LAPACKE_INCLUDES})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LAPACK DEFAULT_MSG
  LAPACK_INCLUDE_DIR LAPACK_LIBRARIES)

MARK_AS_ADVANCED(LAPACK_INCLUDES LAPACK_LIBRARIES)

if(LAPACK_FOUND)
	message(STATUS "LAPACKE C API found.")
else()
	if(LAPACKE_FIND_REQUIRED)
		message(FATAL_ERROR "Required LAPACKE API not found. Please specify library location." )
	else()
		message(STATUS "LAPACKE API not found. Please specify library location." )
	endif()
endif()
