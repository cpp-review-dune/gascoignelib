SET(CONVECTIONDIFFUSION2D_SRC
  localdirichletdata.cc
  localequation.cc 
  localloop.cc
  main1.cc 
  problemdescriptor1.cc
)

IF(BUILD_EXECUTABLES)
  ADD_EXECUTABLE(ConvectionDiffusion2D ${CONVECTIONDIFFUSION2D_SRC})
ENDIF(BUILD_EXECUTABLES)

