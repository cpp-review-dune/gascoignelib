SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/ConvectionDiffusion2D)
ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp ${_SRC}/*.gup ${_SRC}/*.bup .)

ADD_EXECUTABLE(ConvectionDiffusion2D main.cc local.cc localequation.cc)

