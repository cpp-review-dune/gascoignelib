SET(IntegrationTest_SRC
  main.cc 
)

SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/IntegrationTest)
IF(COPY_PARAMFILES)
  ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)
ENDIF(COPY_PARAMFILES)

ADD_EXECUTABLE("IntegrationTest${BUILD_ENDING}" IntegrationTest_SRC)



