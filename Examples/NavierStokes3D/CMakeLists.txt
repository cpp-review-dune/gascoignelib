SET(NAVIERSTOKES3D_SRC
  benchmark3ddirichletdata.cc
  main1.cc
  problemdescriptor1.cc 
)

SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/NavierStokes3D)
ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)

ADD_EXECUTABLE(NavierStokes3D ${NAVIERSTOKES3D_SRC})



