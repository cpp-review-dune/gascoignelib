
SET(MESH_SRC
main_mesh.cc
)
LINK_LIBRARIES(Solver Problems MeshInterpretor LinAlg Mesh Common ${EXT_LIBS})

ADD_EXECUTABLE(mesh MESH_SRC)

SET(COMMON_SRC
localequation.cc
localmeshagent.cc
localtimeloop.cc
)

SET(INITIAL_SRC
${COMMON_SRC}
localinitialcondition.cc
main_initial.cc
problemdescriptorinitial.cc
)
ADD_EXECUTABLE(initial INITIAL_SRC)

SET(FORWARD_SRC
${COMMON_SRC}
main_forward.cc
problemdescriptorforward.cc
)
ADD_EXECUTABLE(forward FORWARD_SRC)

SET(TERMINAL_SRC
${COMMON_SRC}
localterminalcondition.cc
main_terminal.cc
problemdescriptorterminal.cc
)
ADD_EXECUTABLE(terminal TERMINAL_SRC)

SET(BACKWARD_SRC
${COMMON_SRC}
main_backward.cc
backwardequation.cc
problemdescriptorbackward.cc
)
ADD_EXECUTABLE(backward BACKWARD_SRC)
