SET(HEAT2D_SRC
main1.cc 
problemdescriptor1.cc 
localequation.cc 
localinitialcondition.cc
)

IF(BUILD_EXAMPLES)
  ADD_EXECUTABLE(Heat2D HEAT2D_SRC)
ENDIF(BUILD_EXAMPLES)


