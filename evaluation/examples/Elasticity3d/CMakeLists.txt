SET(ELASTICITY_SOURCE 
	${CMAKE_CURRENT_SOURCE_DIR}/elasticityproblem.h
	${CMAKE_CURRENT_SOURCE_DIR}/local.h
)

set(ELASTICITY_3D Elasticity3d)

ADD_EXECUTABLE(${ELASTICITY_3D} ${ELASTICITY_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)

target_compile_definitions(${ELASTICITY_3D} PUBLIC
	${GASCOIGNE_DEFINITIONS}
)

target_include_directories(${ELASTICITY_3D} PUBLIC
	${GASCOIGNE_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${ELASTICITY_3D} 
	GascoigneStd
)

set(ALL_CXX_SOURCE_FILES ${ALL_CXX_SOURCE_FILES} ${ELASTICITY_SOURCE} PARENT_SCOPE)

