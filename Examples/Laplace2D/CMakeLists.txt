SET(LAPLACE2D_SRC_1
main1.cc 
problemdescriptor1.cc 
polynomialexactsolution.cc 
)

SET(LAPLACE2D_SRC_2
main2.cc 
problemdescriptor2.cc 
polynomialexactsolution.cc
neumanndatabyexactsolution.cc
)

IF(BUILD_EXAMPLES)
  ADD_EXECUTABLE(Laplace2D-1 LAPLACE2D_SRC_1)
  ADD_EXECUTABLE(Laplace2D-2 LAPLACE2D_SRC_2) 
ENDIF(BUILD_EXAMPLES)



