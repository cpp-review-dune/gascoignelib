SET(NAVIERSTOKES3D_SRC
  benchmark3ddirichletdata.cc
  main1.cc
  problemdescriptor1.cc 
)

IF(BUILD_EXECUTABLES)
  ADD_EXECUTABLE(NavierStokes3D ${NAVIERSTOKES3D_SRC})
ENDIF(BUILD_EXECUTABLES)



