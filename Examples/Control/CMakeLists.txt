SET(CONTROL_COMMON_SRC
  localequation.cc
  localmeshagent.cc
  localtimeloop.cc
)

SET(CONTROL_MESH_SRC
  main_mesh.cc
)

SET(CONTROL_INITIAL_SRC
  ${CONTROL_COMMON_SRC}
  localinitialcondition.cc
  main_initial.cc
  problemdescriptorinitial.cc
)

SET(CONTROL_FORWARD_SRC
  ${CONTROL_COMMON_SRC}
  main_forward.cc
  problemdescriptorforward.cc
)

SET(CONTROL_TERMINAL_SRC
  ${CONTROL_COMMON_SRC}
  localterminalcondition.cc
  main_terminal.cc
  problemdescriptorterminal.cc
)

SET(CONTROL_BACKWARD_SRC
  ${CONTROL_COMMON_SRC}
  backwardequation.cc
  main_backward.cc
  problemdescriptorbackward.cc
)

IF(BUILD_EXECUTABLES)
  ADD_EXECUTABLE(Control-Mesh     ${CONTROL_MESH_SRC})
  ADD_EXECUTABLE(Control-Initial  ${CONTROL_INITIAL_SRC})
  ADD_EXECUTABLE(Control-Forward  ${CONTROL_FORWARD_SRC})
  ADD_EXECUTABLE(Control-Terminal ${CONTROL_TERMINAL_SRC})
  ADD_EXECUTABLE(Control-Backward ${CONTROL_BACKWARD_SRC})
ENDIF(BUILD_EXECUTABLES)

