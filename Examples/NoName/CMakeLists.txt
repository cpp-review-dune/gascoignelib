#
PROJECT(Local)
#
INCLUDE($ENV{GASCOIGNE}/GascoigneStd/CMakeGlobal.txt)
#
LINK_LIBRARIES(${GASCOIGNE_LIBS})
#
ADD_EXECUTABLE(local 

main.cc 

problemdescriptor.cc 
localdirichletdata.cc
localequation.cc

localmesh2d.cc
localmeshconstructor.cc

localloop.cc
localmultilevelsolver.cc
localsolver.cc
localmeshinterpretor.cc
localintegrator.cc
)

