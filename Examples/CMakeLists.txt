SUBDIRS(
  Laplace2D 
  Laplace3D 
  NavierStokes2D 
  NavierStokes3D 
  Heat2D 
  ConvectionDiffusion2D 
  Control
  IntegrationTest
  NoName
  IntegrationTest
)

LINK_LIBRARIES(${GASCOIGNE_LIBS})
