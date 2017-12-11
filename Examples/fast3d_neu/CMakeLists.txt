ADD_EXECUTABLE("fast2d_neu" 
 alebase.cc  ale.cc  alediscretization.cc  chi.cc  fsi.cc    loop.cc
main.cc  solvers.cc finehierarchicalmesh3d.cc boundaryfsi.cc
)
ADD_EXECUTABLE("fast3d_neu" 
 alebase.cc  ale.cc  alediscretization.cc  chi.cc  fsi.cc    loop.cc
main3d.cc  solvers.cc finehierarchicalmesh3d.cc boundaryfsi.cc
)

SET_TARGET_PROPERTIES(
  fast3d_neu
  PROPERTIES
  COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
  )
SET_TARGET_PROPERTIES(
  fast2d_neu
  PROPERTIES
  COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
  )
