SET(HEAT2D_SRC
  localequation.cc 
  localinitialcondition.cc
  main1.cc 
  problemdescriptor1.cc 
)

IF(BUILD_EXECUTABLES)
  ADD_EXECUTABLE(Heat2D HEAT2D_SRC)
ENDIF(BUILD_EXECUTABLES)


