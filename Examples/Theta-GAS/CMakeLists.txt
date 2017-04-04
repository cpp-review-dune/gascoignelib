ADD_EXECUTABLE("Theta" 
loop.cc
main.cc
myequation.cc
solvers.cc
)

ADD_EXECUTABLE("ThetaExact" 
loop.cc
main_exact.cc
myequation.cc
solvers.cc
)
