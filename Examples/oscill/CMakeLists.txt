ADD_EXECUTABLE("oscill"
umfsplit.cc alebase.cc  ale.cc  alediscretization.cc  chi.cc  fsi.cc div.cc loop.cc
main.cc  solvers.cc
)
SET_TARGET_PROPERTIES(
    oscill
    PROPERTIES
    COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
    )

ADD_EXECUTABLE("box3d"
umfsplit.cc alebase.cc  ale.cc  alediscretization.cc  chi.cc  fsi.cc div.cc loop.cc
main3d.cc  solvers.cc
)
SET_TARGET_PROPERTIES(
    box3d
    PROPERTIES
    COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
    )

ADD_EXECUTABLE("para_box3d"
umfsplit.cc alebase.cc  ale.cc  alediscretization.cc  chi.cc  fsi.cc div.cc loop.cc
mainBox3d.cc  solvers.cc
)
SET_TARGET_PROPERTIES(
    para_box3d
    PROPERTIES
    COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
    )

ADD_EXECUTABLE("paraoscill"
umfsplit.cc alebase.cc  ale.cc  alediscretization.cc  chi.cc  fsi.cc  div.cc  loop.cc
mainPara.cc  solvers.cc
)
SET_TARGET_PROPERTIES(
    paraoscill
    PROPERTIES
    COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
    )

ADD_EXECUTABLE("paraoscill_bench"
umfsplit.cc alebase.cc  ale.cc  alediscretization.cc  chi.cc  fsi.cc  div.cc  loop.cc
benchPara.cc  solvers.cc
)
SET_TARGET_PROPERTIES(
    paraoscill_bench
    PROPERTIES
    COMPILE_DEFINITIONS "${GASCOIGNE_DEFINITIONS}"
    )
