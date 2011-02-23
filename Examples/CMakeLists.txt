SUBDIRS(
   Example1
   Example2
   Example3
   Example4
   Example5
   Example6
   Example7
   Example8
   Example9
   MeshTest
  Laplace2D 
  Laplace3D 
  NavierStokes2D 
  NavierStokes3D 
  NavierStokesFace
  Heat2D 
  ConvectionDiffusion2D 
#  IntegrationTest
  Functionals
  MeshInterpolator
   PeriodicStokes2D
)

LINK_LIBRARIES(${GASCOIGNE_LIBS})
