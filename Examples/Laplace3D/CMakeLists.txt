SET(LAPLACE3D_SRC
  main1.cc 
  problemdescriptor1.cc 
  polynomialexactsolution3d.cc
)

IF(BUILD_EXECUTABLES)
  ADD_EXECUTABLE(Laplace3D ${LAPLACE3D_SRC})
ENDIF(BUILD_EXECUTABLES)

