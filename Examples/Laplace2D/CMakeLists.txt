SET(LAPLACE2D_1_SRC
  main1.cc 
  problemdescriptor1.cc 
  polynomialexactsolution.cc 
)

SET(LAPLACE2D_2_SRC
  main2.cc 
  neumanndatabyexactsolution.cc
  problemdescriptor2.cc 
  polynomialexactsolution.cc
)

SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/Laplace2D)
ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)

ADD_EXECUTABLE(Laplace2D-1 ${LAPLACE2D_1_SRC})
ADD_EXECUTABLE(Laplace2D-2 ${LAPLACE2D_2_SRC}) 



