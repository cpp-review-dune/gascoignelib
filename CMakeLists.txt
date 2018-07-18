PROJECT(GASCOIGNESTD CXX C)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#
# .gascoignerc_newcmake contains the user dependent settings
#

INCLUDE($ENV{HOME}/.gascoignerc_newcmake)
#
# CMakeGlobal.txt is intended for use in derived projects
#

INCLUDE(CMakeGlobal.txt)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(Examples)

