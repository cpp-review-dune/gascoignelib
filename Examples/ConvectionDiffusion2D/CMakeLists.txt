SET(CONVECTIONDIFFUSION2D_SRC
main1.cc 
problemdescriptor1.cc
localdirichletdata.cc
localequation.cc 
localloop.cc
localmeshinterpretor.cc
)

IF(BUILD_EXAMPLES)
  ADD_EXECUTABLE(ConvectionDiffusion2D CONVECTIONDIFFUSION2D_SRC)
ENDIF(BUILD_EXAMPLES)

