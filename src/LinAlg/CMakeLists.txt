INCLUDE_DIRECTORIES(
${STD_SRC}/Common
${STD_SRC}/Mesh
${STD_SRC}/MeshInterpretor
${STD_SRC}/Solver
)

SET(LINALG_SRC 
columnstencil.cc
columndiagstencil.cc
twinstencil.cc
simplematrix.cc
nodesparsestructureadaptor.cc
sparsestructureadaptor.cc
pointilu.cc
pointmatrix.cc
simpleilu.cc
simplesparsestructureadaptor.cc
entrymatrix.cc
sparsestructure.cc
umfilu.cc
)

ADD_LIBRARY(LinAlg ${LINALG_SRC})
