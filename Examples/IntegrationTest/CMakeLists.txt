SET(IntegrationTest_SRC
  main.cc 
)

IF(BUILD_EXECUTABLES)
  ADD_EXECUTABLE(IntegrationTest IntegrationTest_SRC)
ENDIF(BUILD_EXECUTABLES)



