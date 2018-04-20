INCLUDE_DIRECTORIES(${GASCOIGNE_INCLUDE_DIRS})

IF(NOT "${EXAMPLE_OUTPUT_PATH}" STREQUAL "")
  SET(EXECUTABLE_OUTPUT_PATH "${EXAMPLE_OUTPUT_PATH}")
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GASCOIGNE_CXX_FLAGS}")
SET(COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}")

SUBDIRS(
   fsi3dbench
	 Vanka
	 Solid3D
	 SolidEuler
	 Vein
   DG
   fast3d_neu
   Example1
   Example2
   Example3
   Example4
   Example5
   Example6
   Example7
   Example8
   Example9
   MeshTest
  Laplace2D 
  Laplace3D 
  NavierStokes2D 
  NavierStokes3D 
  NavierStokesFace
  Heat2D 
  ConvectionDiffusion2D 
#  IntegrationTest
  Functionals
  MeshInterpolator
#   PeriodicStokes2D
)

LINK_LIBRARIES(${GASCOIGNE_LIBS})
