SET(NAVIERSTOKES2D_SRC
main1.cc 
problemdescriptor1.cc 
benchmarkdirichletdata.cc
)

IF(BUILD_EXAMPLES)
  ADD_EXECUTABLE(NavierStokes2D NAVIERSTOKES2D_SRC)
ENDIF(BUILD_EXAMPLES)


