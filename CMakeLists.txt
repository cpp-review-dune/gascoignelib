#-----------------------------#
PROJECT(GASCOIGNESTD)
#-----------------------------#

#
# CMakeGlobal.txt is intentended for use in derived projects
#

INCLUDE(CMakeGlobal.txt)

#
# the following is local
#

OPTION(BUILD_EXAMPLES "Build the executable examples" ON)

SUBDIRS(src Examples)
