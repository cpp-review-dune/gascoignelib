SUBDIRS(
  NavierStokes2D  bench Transport_split   DG_ice SUPG DG_test

)

LINK_LIBRARIES(${GASCOIGNE_LIBRARIES} )
