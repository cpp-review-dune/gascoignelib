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

SET(LOCAL_SRC
  main.cc 
)

ADD_EXECUTABLE(Local ${LOCAL_SRC})
