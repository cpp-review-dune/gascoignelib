
INCLUDE_DIRECTORIES(${GASCOIGNE_INCLUDE_DIRS})


SUBDIRS(
  visco
  Laplace2D 
  Laplace3D 
  NavierStokes2D 
  NavierStokes3D 
  NavierStokesFace
  Heat2D 
  ConvectionDiffusion2D 
  Functionals
  MeshInterpolator
)

LINK_LIBRARIES(${GASCOIGNE_LIBS})
