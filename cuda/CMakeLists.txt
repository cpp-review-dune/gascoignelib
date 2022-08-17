OPTION(USE_CLANG_TIDY "Use clang-tidy warnings" OFF)
if(USE_CLANG_TIDY)
    set( CMAKE_CXX_CLANG_TIDY "clang-tidy;-checks=*,-fuchsia-*")
endif()

SET(CUDA_SDK_ROOT_DIR /usr/local/cuda)
SET(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda)
find_package(CUDA 10.1 REQUIRED)

FILE(GLOB_RECURSE GASCOIGNE_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/../src/**/*.cc)
FILE(GLOB_RECURSE GASCUDA_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/*.cc ${CMAKE_CURRENT_SOURCE_DIR}/*.cu)
LIST(APPEND GASCOIGNE_SOURCE ${GASCUDA_SOURCE})

list(APPEND CUDA_NVCC_FLAGS "-arch=sm_60") #Volta Architecture
list(APPEND CUDA_NVCC_FLAGS "-maxrregcount=64")
#list(APPEND CUDA_NVCC_FLAGS "--ptxas-options=-v")

CUDA_ADD_LIBRARY( ${GASCOIGNE_LIBRARY} ${GASCOIGNE_SOURCE} SHARED )

target_compile_options(${GASCOIGNE_LIBRARY} PUBLIC
	-std=c++11
	-pedantic
	-Wextra
	-Wconversion
	-Wold-style-cast
)

set_target_properties(${GASCOIGNE_LIBRARY} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${GASCOIGNE_LIBRARY} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

IF(USE_CUDA)
    LIST(APPEND GASCOIGNE_DEFINITIONS -DUSE_CUDA)
ENDIF(USE_CUDA)

target_compile_definitions(${GASCOIGNE_LIBRARY} PUBLIC
    ${GASCOIGNE_DEFINITIONS}
)

target_include_directories(${GASCOIGNE_LIBRARY} SYSTEM PUBLIC
	${CUDA_INCLUDE_DIRS}
	${GASCOIGNE_INCLUDE_DIRS}
	${GASCOIGNE_DIR}/cuda
	${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${GASCOIGNE_LIBRARY}
	${GASCOIGNE_LIBRARIES}
	${CUDA_LIBRARIES}
	${CUDA_cusparse_LIBRARY}
	${CUDA_CUBLAS_LIBRARIES}
)

OPTION(USE_IWYU "Use include what you use for include optimization" OFF)
IF(USE_IWYU)
    find_program(IWYU include-what-you-use HINTS /usr/bin/)
    if(NOT IWYU)
        message("IWYU not found!")
    endif()

    set_property(TARGET ${GASCOIGNE_LIBRARY} 
        PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU})
endif()