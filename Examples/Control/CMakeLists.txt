
SET(MESH_SRC
main_mesh.cc
)

ADD_EXECUTABLE(Control-Mesh MESH_SRC)

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
ADD_EXECUTABLE(Control-Initial INITIAL_SRC)

SET(FORWARD_SRC
${COMMON_SRC}
main_forward.cc
problemdescriptorforward.cc
)
ADD_EXECUTABLE(Control-Forward FORWARD_SRC)

SET(TERMINAL_SRC
${COMMON_SRC}
localterminalcondition.cc
main_terminal.cc
problemdescriptorterminal.cc
)
ADD_EXECUTABLE(Control-Terminal TERMINAL_SRC)

SET(BACKWARD_SRC
${COMMON_SRC}
main_backward.cc
backwardequation.cc
problemdescriptorbackward.cc
)
ADD_EXECUTABLE(Control-Backward BACKWARD_SRC)
