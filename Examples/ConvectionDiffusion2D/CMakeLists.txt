SET(CONVECTIONDIFFUSION2D_SRC
  localdirichletdata.cc
  localequation.cc 
  localloop.cc
  main1.cc 
  problemdescriptor1.cc
)

SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/ConvectionDiffusion2D)
IF(COPY_PARAMFILES)
  ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param .)
ENDIF(COPY_PARAMFILES)

ADD_EXECUTABLE(ConvectionDiffusion2D ${CONVECTIONDIFFUSION2D_SRC})

