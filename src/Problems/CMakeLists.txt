SET(PROBLEMS_SRC 
  laplace3d.cc
  laplace2d.cc
  stokes3d.cc
  stokes2d.cc
  stokesgls3d.cc
  stokesgls2d.cc
  navierstokes3d.cc
  navierstokes2d.cc
  navierstokesgls3d.cc
  navierstokesgls2d.cc

  boundarymanager.cc
  constantboundaryfunctional.cc
  constantrighthandside.cc
  dirichletdatabycolor.cc
  domainmeanfunctional.cc
  exactsolution.cc
  glsstabilization.cc
  pointfunctional.cc
  residualfunctional.cc
  usefullfunctionsbd.cc
)

ADD_LIBRARY(Problems PROBLEMS_SRC)
