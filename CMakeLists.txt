#-----------------------------#
PROJECT(GASCOIGNESTD)
#-----------------------------#

#
# .gascoignerc contains the user depending settings
#

INCLUDE($ENV{HOME}/.gascoignerc)

#
# CMakeGlobal.txt is intentended for use in derived projects
#

INCLUDE(CMakeGlobal.txt)

#
# the following is local
#

SUBDIRS(src Examples)
