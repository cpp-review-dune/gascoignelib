SET(SOLVER_SRC 
  basicloop.cc
  cg.cc
  cuthillmckee.cc
  ilupermutate.cc
  monitor.cc
  multilevelsolverdata.cc
  newtonstdmultilevelsolver.cc
  pi.cc
  solverdata.cc
  starter.cc
  stdiomanager.cc
  stdloop.cc
  stdmultilevelsolver.cc
  stdsolver.cc
  stdtimeloop.cc
  stdtimesolver.cc
  timeinfo.cc
)

ADD_LIBRARY(Solver SOLVER_SRC)
