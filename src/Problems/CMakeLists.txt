SET(PROBLEMS_SRC 
  boundarymanager.cc
  constantboundaryfunctional.cc
  constantrighthandside.cc
  dirichletdatabycolor.cc
  domainmeanfunctional.cc
  exactsolution.cc
  functionalmanager.cc
  functionalmanagerconstruct.cc
  glsstabilization.cc
  laplace3d.cc
  laplace.cc
  navierstokes3d.cc
  navierstokes.cc
  navierstokesgls3d.cc
  navierstokesgls.cc
  pointfunctional.cc
  residualfunctional.cc
  stokes3d.cc
  stokes.cc
  stokesgls3d.cc
  stokesgls.cc
  usefullfunctionsbd.cc
)

ADD_LIBRARY(Problems PROBLEMS_SRC)
