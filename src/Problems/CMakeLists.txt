INCLUDE_DIRECTORIES(
${STD_SRC}/Common
${STD_SRC}/Mesh
${STD_SRC}/MeshInterpretor
)

SET(PROBLEMS_SRC 
boundarymanager.cc
constantboundaryfunctional.cc
constantrighthandside.cc
dirichletdatabycolor.cc
domainmeanfunctional.cc
exactsolution.cc
functionalmanager.cc
functionalmanagerconstruct.cc
laplace.cc
laplace3d.cc
navierstokes3d.cc
navierstokes.cc
navierstokesgls.cc
navierstokesgls3d.cc
pointfunctional.cc
residualfunctional.cc
glsstabilization.cc
stokes3d.cc
stokes.cc
stokesgls.cc
stokesgls3d.cc
usefullfunctionsbd.cc
)

ADD_LIBRARY(Problems ${PROBLEMS_SRC})
