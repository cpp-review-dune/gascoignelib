SET(IntegrationTest_SRC
  main.cc 
)

SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/IntegrationTest)
ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)

ADD_EXECUTABLE(IntegrationTest IntegrationTest_SRC)



