set(ADAPTIVE_ELASTICITY_2D AdaptiveElasticity2d)

SET(ELASTICITY_SOURCE 
	${CMAKE_CURRENT_SOURCE_DIR}/elasticityproblem.h
	${CMAKE_CURRENT_SOURCE_DIR}/local.h
)

ADD_EXECUTABLE(${ADAPTIVE_ELASTICITY_2D} ${ELASTICITY_SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)

target_compile_definitions(${ADAPTIVE_ELASTICITY_2D} PUBLIC
	${GASCOIGNE_DEFINITIONS}
)

target_include_directories(${ADAPTIVE_ELASTICITY_2D} PUBLIC
	${GASCOIGNE_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${ADAPTIVE_ELASTICITY_2D} 
	GascoigneStd
)

set(ALL_CXX_SOURCE_FILES ${ALL_CXX_SOURCE_FILES} ${ELASTICITY_SOURCE} PARENT_SCOPE)

