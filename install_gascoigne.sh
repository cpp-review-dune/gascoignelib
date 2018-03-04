#! /bin/sh

BUILDDIR=build

## check that we are in the Gascoigne Directory
if [ ! -f LICENSE.TXT ]; then
    echo "Call installation script only in main Gascoigne directory"
    exit;
fi

## check that is the right file
firstline=`head -1 LICENSE.TXT`
if [ ! "${firstline}" == '###GASCOIGNE3D###' ]; then
    echo "Call installation script only in main Gascoigne directory"
    echo "File LICENSE.TXT is not the Gascoigne 3d library license file"
    exit;
fi



## create directory to build gascoigne
if [ -d $BUILDDIR ]; then
    echo "error: the build directiory '${BUILDDIR}' already exists. Remove it!"
    exit;
fi

echo "gascoigne: make build directory: ${BUILDDIR}"
mkdir $BUILDDIR

echo "gascoigne: configure ..."
cd $BUILDDIR
cmake ../ 1>/dev/null

echo "gascoigne: make ..."
make -j4
