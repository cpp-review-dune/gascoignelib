ADD_EXECUTABLE("Burger_gekoppelt" main.cc loop.cc myequation.cc solvers.cc)



SET_TARGET_PROPERTIES(
  Burger_gekoppelt
  PROPERTIES
  COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
  )
