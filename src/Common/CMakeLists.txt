SET(COMMON_SRC 
  backup.cc
  cginfo.cc
  compose_name.cc
  compvector.cc
  dataformathandler.cc
  extrapolator.cc
  filescanner.cc
  fixarray.cc
  linescanner.cc
  nlinfo.cc
  nmatrix.cc
  nvector.cc
  set2vec.cc
  stlio.cc
  stopwatch.cc
  stringutil.cc
  timer.cc
  vecalgo.cc
)

ADD_LIBRARY(Common COMMON_SRC)
