SUBDIRS(
  Laplace2D 
  Laplace3D 
  NavierStokes2D 
  NavierStokes3D 
  Heat2D 
  ConvectionDiffusion2D 
  IntegrationTest
  Functionals
  MeshInterpolator
#  LLPS
#  Dung
)

LINK_LIBRARIES(${GASCOIGNE_LIBS})
