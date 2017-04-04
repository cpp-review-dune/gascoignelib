ADD_EXECUTABLE("ThetaNS"
ns.cc   loop.cc main.cc solvers.cc functionals.cc)

ADD_EXECUTABLE("ThetaNSExact"
ns.cc   loop.cc main_exact.cc solvers.cc functionals.cc)
