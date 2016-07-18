INCLUDE_DIRECTORIES(${GASCOIGNE_INCLUDE_DIRS})


SUBDIRS(
  bench
  Laplace2D 
  NavierStokes2D 
  NavierStokes3D 
  NavierStokesInstat
)

LINK_LIBRARIES(${GASCOIGNE_LIBRARY})
