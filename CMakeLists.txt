set(CMAKE_CXX_FLAGS_RELEASE "_UNSET" CACHE STRING "")
PROJECT(GASCOIGNE CXX C)

#
# The name of the library:
#
SET(GASCOIGNE_LIBRARY GascoigneStd)




#### CMAKE CONFIG
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
cmake_policy(SET CMP0042 NEW)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)


#### Import all Options with standard values that are also used in dependent projects
INCLUDE(CMakeOptions.txt)

#### Directory Settings

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
#SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../lib"  CACHE PATH "output path for libraries")
SET (GASCOIGNE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#### Options for Gascoigne Core Library
#### not necessarily used in derived projects
OPTION(WITH_UMFPACK      "Build with support for UMFPACK" ON)
OPTION(WITH_UMFPACK_LONG "Long-integer indices in UMFPACK >4GB matrices" OFF)


#
# Include directories:
#
INCLUDE(CMakeLists_Includes.txt)

#
# Setup compiler flags and compile definitions:
#
SET(GASCOIGNE_CXX_FLAGS
  "-std=c++11 -pedantic -Wall -Wno-sign-compare -Wold-style-cast -Wno-unknown-pragmas" CACHE STRING "Configuration flags"
  )


OPTION(WITH_THREADS      "Build with support for openmp" OFF)
IF(WITH_THREADS)

  # require OPENMP
  
  find_package(OpenMP REQUIRED)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  SET(GASCOIGNE_CXX_FLAGS "${GASCOIGNE_CXX_FLAGS} -fopenmp")

endif()


#
# Search for metis:
#
FIND_PACKAGE(METIS REQUIRED)

LIST(APPEND GASCOIGNE_INCLUDE_DIRS ${METIS_INCLUDE_DIRS})
LIST(APPEND GASCOIGNE_LIBRARIES ${METIS_LIBRARIES})

#
# Search for UMFPACK:
#
IF(WITH_UMFPACK)
  FIND_PACKAGE(UMFPACK REQUIRED)

  LIST(APPEND GASCOIGNE_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIRS})
  LIST(APPEND GASCOIGNE_LIBRARIES ${UMFPACK_LIBRARIES})
  SET(GASCOIGNE_LINKER_FLAGS "${GASCOIGNE_LINKER_FLAGS} ${UMFPACK_LINKER_FLAGS}")
  LIST(APPEND GASCOIGNE_DEFINITIONS -D __WITH_UMFPACK__)

  if (WITH_UMFPACK_LONG)
    LIST(APPEND GASCOIGNE_DEFINITIONS -D __WITH_UMFPACK_LONG__)
  ENDIF()
ENDIF()

#
# Search for blas:
#
FIND_PACKAGE(BLAS REQUIRED)

LIST(APPEND GASCOIGNE_LIBRARIES ${BLAS_LIBRARIES})
SET(GASCOIGNE_LINKER_FLAGS "${GASCOIGNE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")

#
# Check whether tr1/unordered_map is available:
#
INCLUDE(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <tr1/unordered_map>
  int main()
  {
  std::tr1::unordered_map<int, int> foo;
  return 0;
  }
  "
  HAVE_TR1_UNORDERED_MAP
  )

IF(HAVE_TR1_UNORDERED_MAP)
  LIST(APPEND GASCOIGNE_DEFINITIONS -D __NEWER_THAN_GCC_4_2__)
ENDIF()

SEPARATE_ARGUMENTS(GASCOIGNE_CXX_FLAGS)
ADD_COMPILE_OPTIONS(${GASCOIGNE_CXX_FLAGS})
ADD_DEFINITIONS(${GASCOIGNE_DEFINITIONS})

SET(GASCOIGNE_LIBS "GascoigneStd${BUILD_ENDING}")

#
# Code Quality
#

OPTION(USE_IWYU "Use include what you use for include optimization" OFF)
IF(USE_IWYU)
    find_program(IWYU include-what-you-use HINTS /usr/bin/)
    if(NOT IWYU)
        message("IWYU not found!")
    endif()

    set_property(TARGET ${GASCOIGNE_LIBRARY} 
        PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU})
endif()

OPTION(USE_CLANG_FORMAT "Formats the codebase" ON)
if(USE_CLANG_FORMAT)
  FILE(GLOB_RECURSE ALL_CXX_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/**/*.cc ${CMAKE_CURRENT_SOURCE_DIR}/**/*.h ${CMAKE_CURRENT_SOURCE_DIR}/**/*.xx ${CMAKE_CURRENT_SOURCE_DIR}/**/*.cu)
  # Adding clang-format target if executable is found
  find_program(CLANG_FORMAT "clang-format")
  if(CLANG_FORMAT)
    add_custom_target(
      clang-format ALL
      COMMAND ${CLANG_FORMAT}
      -i
      -style=Mozilla
      ${ALL_CXX_SOURCE_FILES}
      )
  else()
    message(STATUS "Could not find clang-format")
  endif()
endif()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(doc)
