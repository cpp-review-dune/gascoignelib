SET(CONVECTIONDIFFUSION2D_SRC
  localdirichletdata.cc
  localequation.cc 
  localloop.cc
  main1.cc 
  problemdescriptor1.cc
)

SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/ConvectionDiffusion2D)
ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param .)

ADD_EXECUTABLE(ConvectionDiffusion2D ${CONVECTIONDIFFUSION2D_SRC})

