ADD_EXECUTABLE("ChickenBurger" main.cc loop.cc myequation.cc solvers.cc)




SET_TARGET_PROPERTIES(
  ChickenBurger
  PROPERTIES
  COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
  )
