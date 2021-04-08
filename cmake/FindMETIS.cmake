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
# Try to find the (serial) METIS library
#
# This module exports
#
#   METIS_LIBRARIES
#   METIS_INCLUDE_DIRS
#

INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(METIS_INCLUDE_DIR metis.h
  HINTS
    ${METIS_DIR}
  PATH_SUFFIXES
    metis include/metis include
  )

FIND_LIBRARY(METIS_LIBRARY
  NAMES metis
  HINTS
    ${METIS_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
  )

FIND_LIBRARY(PARMETIS_LIBRARY
  NAMES parmetis
  HINTS
    ${METIS_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(METIS DEFAULT_MSG
  METIS_LIBRARY
  METIS_INCLUDE_DIR
  )

MARK_AS_ADVANCED(
  METIS_LIBRARY
  PARMETIS_LIBRARY
  METIS_INCLUDE_DIR
  )

IF(METIS_FOUND)
  #
  # Sanity check: Only include parmetis library if it is in the same
  # directory as the metis library...
  #
  GET_FILENAME_COMPONENT(_path1 "${PARMETIS_LIBRARY}" PATH)
  GET_FILENAME_COMPONENT(_path2 "${ETIS_LIBRARY}" PATH)
  IF( NOT PARMETIS_LIBRARY MATCHES "-NOTFOUND"
      AND "${_path1}" STREQUAL "${_path2}" )
    SET(METIS_LIBRARIES ${PARMETIS_LIBRARY})
  ENDIF()

  LIST(APPEND METIS_LIBRARIES
    ${METIS_LIBRARY}
    ${MPI_C_LIBRARIES} # for good measure
    )
  SET(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
ENDIF()
