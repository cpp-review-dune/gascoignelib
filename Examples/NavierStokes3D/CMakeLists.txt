SET(NAVIERSTOKES3D_SRC
main1.cc
problemdescriptor1.cc 
benchmark3ddirichletdata.cc
)

IF(BUILD_EXAMPLES)
  ADD_EXECUTABLE(NavierStokes3D NAVIERSTOKES3D_SRC)
ENDIF(BUILD_EXAMPLES)



