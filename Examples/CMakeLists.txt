SUBDIRS(
  NavierStokes2D  bench Burger Nichtlinear Laplace_t Problems Crank
)

LINK_LIBRARIES(${GASCOIGNE_LIBRARIES} )
