SET(LINALG_SRC 
  columndiagstencil.cc
  columnstencil.cc
  entrymatrix.cc
  nodesparsestructureadaptor.cc
  pointilu.cc
  pointmatrix.cc
  simpleilu.cc
  simplematrix.cc
  simplesparsestructureadaptor.cc
  sparsestructureadaptor.cc
  sparsestructure.cc
  twinstencil.cc
  umfilu.cc
)

ADD_LIBRARY(LinAlg LINALG_SRC)
