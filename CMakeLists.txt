#-----------------------------#
PROJECT(GASCOIGNESTD)
#-----------------------------#

#
# .gascoignerc contains the user dependent settings
#

INCLUDE($ENV{HOME}/.gascoignerc)

#
# CMakeGlobal.txt is intended for use in derived projects
#

INCLUDE(CMakeGlobal.txt)

#
# the following is local
#

SUBDIRS(src test Examples)
