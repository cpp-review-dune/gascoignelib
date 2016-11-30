INCLUDE_DIRECTORIES(${GASCOIGNE_INCLUDE_DIRS})


SUBDIRS(
  Laplace2D 
  NavierStokes2D 
  NavierStokes3D 
  NavierStokesInstat
  seaice
  dg
)

LINK_LIBRARIES(${GASCOIGNE_LIBRARY})
