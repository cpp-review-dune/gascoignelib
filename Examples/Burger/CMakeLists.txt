ADD_EXECUTABLE("Burger" main.cc loop.cc myequation.cc solvers.cc)



SET_TARGET_PROPERTIES(
  Burger
  PROPERTIES
  COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
  )
