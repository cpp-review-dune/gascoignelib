SET(HEAT2D_SRC
  localequation.cc 
  localinitialcondition.cc
  main1.cc 
  problemdescriptor1.cc 
)

SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/Heat2D)
ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)

ADD_EXECUTABLE(Heat2D ${HEAT2D_SRC})


