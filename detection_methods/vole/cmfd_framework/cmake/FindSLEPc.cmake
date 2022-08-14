# -*- mode: cmake -*-
#
#  This file is part of the Feel library
#
#  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
#       Date: 2010-01-22
#
#  Copyright (C) 2010 Universit� Joseph Fourier
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 3.0 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
find_path (SLEPC_DIR include/slepc.h
  HINTS ENV SLEPC_DIR
  PATHS
  "/net/cv/lib64/slepcdir"
  "/usr/lib/slepcdir/3.0.0" # Debian
  $ENV{HOME}/slepc
  DOC "SLEPc Directory"
)

SET(SLEPC_INCLUDE_DIR "${SLEPC_DIR}/include/")
#INCLUDE(CheckIncludeFile)
#CHECK_INCLUDE_FILE(${SLEPC_INCLUDE_DIR}/slepc.h HAVE_SLEPC_H)
FIND_LIBRARY(SLEPC_LIB_SLEPC slepc PATHS "${SLEPC_DIR}/lib")
SET(SLEPC_LIBRARIES ${SLEPC_LIB_SLEPC} CACHE STRING "SLEPc libraries" FORCE)
#if (HAVE_SLEPC_H AND SLEPC_DIR AND SLEPC_LIBRARIES )
if (SLEPC_DIR AND SLEPC_LIBRARIES)
  set(HAVE_SLEPC 1)
  set(SLEPC_FOUND ON)
endif()
set(SLEPC_INCLUDES ${SLEPC_INCLUDE_DIR} CACHE STRING "SLEPc include path" FORCE)
MARK_AS_ADVANCED(SLEPC_DIR SLEPC_LIB_SLEPC SLEPC_INCLUDES SLEPC_LIBRARIES)
