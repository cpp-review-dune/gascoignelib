#-----------------------------#
PROJECT(LOCAL)
#-----------------------------#

#
# .gascoignerc contains the user depending settings
#

INCLUDE($ENV{HOME}/.gascoignerc)

#
# CMakeGlobal.txt is intentended for use in derived projects
#

INCLUDE(${GASCOIGNE_STD}/CMakeGlobal.txt)

LINK_LIBRARIES(${GASCOIGNE_LIBS})

IF(COPY_PARAMFILES)
  ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)
ENDIF(COPY_PARAMFILES)

ADD_EXECUTABLE(NoName main.cc)
