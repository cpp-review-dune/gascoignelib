SET(NAVIERSTOKES2D_SRC
  benchmarkdirichletdata.cc
  main1.cc 
  problemdescriptor1.cc 
)

IF(BUILD_EXECUTABLES)
  ADD_EXECUTABLE(NavierStokes2D NAVIERSTOKES2D_SRC)
ENDIF(BUILD_EXECUTABLES)


