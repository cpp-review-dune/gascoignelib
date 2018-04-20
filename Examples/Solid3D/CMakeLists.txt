ADD_EXECUTABLE("Gascoigne_Solid3D${BUILD_ENDING}"    solid.cc    loop.cc
main3d.cc  finehierarchicalmesh3d.cc boundarysolid.cc localmeshagent.cc localhierarchicalmesh3d.cc solvers.cc vanka_matrix_vector.cc)

SET_TARGET_PROPERTIES(
  Gascoigne_Solid3D
  PROPERTIES
  COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
  )

