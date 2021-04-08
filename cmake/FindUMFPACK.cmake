#
# Copyright (C) 2012 - 2014 by Matthias Maier
#
# This file is part of Gascoigne 3D
#
# Gascoigne 3D is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#

#
# Try to find the UMFPACK library
#
# This module exports
#
#   UMFPACK_LIBRARIES
#   UMFPACK_INCLUDE_DIRS
#   UMFPACK_LINKER_FLAGS
#

INCLUDE(FindPackageHandleStandardArgs)
INCLUDE(CheckCXXSourceCompiles)

FIND_PACKAGE(LAPACK)

#
# Two macros to make life easier:
#
MACRO(FIND_UMFPACK_PATH _comp _file)
  STRING(TOLOWER ${_comp} _comp_lowercase)
  STRING(TOUPPER ${_comp} _comp_uppercase)
  FIND_PATH(${_comp}_INCLUDE_DIR ${_file}
    HINTS
      ${${_comp_uppercase}_DIR}
      ${SUITESPARSE_DIR}/${_comp}
      ${UMFPACK_DIR}/../${_comp}
      ${UMFPACK_DIR}/${_comp}
      ${UMFPACK_DIR}
    PATH_SUFFIXES
      ${_comp_lowercase} include/${_comp_lowercase} include Include ${_comp}/Include suitesparse
    )
ENDMACRO()

MACRO(FIND_UMFPACK_LIBRARY _comp _name)
  STRING(TOUPPER ${_comp} _comp_uppercase)
  FIND_LIBRARY(${_comp}_LIBRARY
    NAMES ${_name}
    HINTS
      ${${_comp_uppercase}_DIR}
      ${SUITESPARSE_DIR}/${_comp}
      ${UMFPACK_DIR}/../${_comp}
      ${UMFPACK_DIR}/${_comp}
      ${UMFPACK_DIR}
    PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib Lib ${_comp}/Lib
    )
  IF(NOT "${ARGN}" STREQUAL "REQUIRED")
    IF(${_comp}_LIBRARY MATCHES "-NOTFOUND")
      SET(${_comp}_LIBRARY "")
      UNSET(${_comp}_LIBRARY CACHE)
    ENDIF()
  ENDIF()
ENDMACRO()

#
# Search for include directories:
#
FIND_UMFPACK_PATH(UMFPACK umfpack.h)
FIND_UMFPACK_PATH(AMD amd.h)

#
# Well, recent versions of UMFPACK include SuiteSparse_config.h, if so,
# ensure that we'll find these headers as well.
#
IF(NOT UMFPACK_INCLUDE_DIR MATCHES "-NOTFOUND")
  FILE(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_SUITESPARSE_STRING
    REGEX "#include \"SuiteSparse_config.h\"")
  IF(NOT "${UMFPACK_SUITESPARSE_STRING}" STREQUAL "")
    FIND_UMFPACK_PATH(SuiteSparse_config SuiteSparse_config.h)
    LIST(APPEND _required SuiteSparse_config_INCLUDE_DIR)
  ENDIF()
ENDIF()

#
# Link against everything we can find to avoid underlinkage:
#
FIND_UMFPACK_LIBRARY(UMFPACK umfpack REQUIRED)
FIND_UMFPACK_LIBRARY(AMD amd REQUIRED)
FIND_UMFPACK_LIBRARY(CHOLMOD cholmod)
FIND_UMFPACK_LIBRARY(COLAMD colamd)
FIND_UMFPACK_LIBRARY(CCOLAMD ccolamd)
FIND_UMFPACK_LIBRARY(CAMD camd)
FIND_UMFPACK_LIBRARY(SuiteSparse_config suitesparseconfig)

#
# Test whether libsuitesparseconfig.xxx can be used for shared library
# linkage. If not, exclude it from the command line.
#
LIST(APPEND CMAKE_REQUIRED_LIBRARIES
  "-shared"
  ${SuiteSparse_config_LIBRARY}
  )
CHECK_CXX_SOURCE_COMPILES("extern int SuiteSparse_version (int[3]);
  void foo(int bar[3]) { SuiteSparse_version(bar);}"
  LAPACK_SUITESPARSECONFIG_WITH_PIC
  )
SET(CMAKE_REQUIRED_LIBRARIES)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(UMFPACK DEFAULT_MSG
  UMFPACK_LIBRARY
  AMD_LIBRARY
  UMFPACK_INCLUDE_DIR
  AMD_INCLUDE_DIR
  ${_required}
  LAPACK_FOUND
  )

MARK_AS_ADVANCED(
  AMD_INCLUDE_DIR AMD_LIBRARY atlas_LIBRARY blas_LIBRARY CAMD_LIBRARY
  CHOLMOD_LIBRARY CCOLAMD_LIBRARY COLAMD_LIBRARY
  SuiteSparse_config_INCLUDE_DIR SuiteSparse_config_LIBRARY
  UMFPACK_INCLUDE_DIR UMFPACK_LIBRARY
  )

IF(UMFPACK_FOUND)
  SET(UMFPACK_LIBRARIES
    ${UMFPACK_LIBRARY}
    ${CHOLMOD_LIBRARY}
    ${CCOLAMD_LIBRARY}
    ${COLAMD_LIBRARY}
    ${CAMD_LIBRARY}
    ${AMD_LIBRARY}
    ${LAPACK_LIBRARIES}
    )
  IF(LAPACK_SUITESPARSECONFIG_WITH_PIC OR NOT BUILD_SHARED_LIBS)
    LIST(APPEND UMFPACK_LIBRARIES ${SuiteSparse_config_LIBRARY})
  ENDIF()

  FIND_LIBRARY(rt_LIBRARY NAMES rt)
  MARK_AS_ADVANCED(rt_LIBRARY)
  IF(NOT rt_LIBRARY MATCHES "-NOTFOUND")
    LIST(APPEND UMFPACK_LIBRARIES ${rt_LIBRARY})
  ENDIF()

  SET(UMFPACK_INCLUDE_DIRS
    ${UMFPACK_INCLUDE_DIR}
    ${AMD_INCLUDE_DIR}
    ${SuiteSparse_config_INCLUDE_DIR}
    )

  SET(UMFPACK_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
ENDIF()
