INCLUDE_DIRECTORIES(
${STD_SRC}/Common
${STD_SRC}/Mesh
${STD_SRC}/MeshInterpretor
${STD_SRC}/Problems
${STD_SRC}/LinAlg
)

SET(SOLVER_SRC 
mginterpolator.cc
stdsolver.cc
stdmultilevelsolver.cc
stdloop.cc
timeinfo.cc
stdtimeloop.cc
cg.cc
cuthillmckee.cc
gmres.cc
ilupermutate.cc
monitor.cc
multilevelsolverdata.cc
newtonstdmultilevelsolver.cc
pi.cc
stdiomanager.cc
solverdata.cc
stdtimesolver.cc
starter.cc
mginterpolatornested.cc
)

ADD_LIBRARY(Solver ${SOLVER_SRC})
