SET(LAPLACE3D_SRC
  main1.cc 
  problemdescriptor1.cc 
  polynomialexactsolution3d.cc
)

SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/Laplace3D)
IF(COPY_PARAMFILES)
  ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)
ENDIF(COPY_PARAMFILES)

ADD_EXECUTABLE(Laplace3D ${LAPLACE3D_SRC})

