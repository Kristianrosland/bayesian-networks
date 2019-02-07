#! /bin/bash

function usage() {
   echo "USAGE: ./configure.sh <SCIPDIR>"
   echo "   where <SCIPDIR> is the path to the SCIP directory"
   echo "   For example"
   echo "      ./configure.sh /opt/scipoptsuite-3.0.0/scip-3.0.0"
   echo "   or"
   echo "      ./configure.sh ../../ziboptsuite-2.1.1/scip-2.1.1"
   exit -1
}

function errmsg() {
   echo $1
   exit -1
}

### Check there is the correct number of input arguments
if [ $# -ne 1 ]; then
   usage
fi

### Check the input argument is a valid directory
SCIPDIR=$1
if [ ! -d $SCIPDIR ]; then
   errmsg "ERROR: Couldn't find \"$SCIPDIR\""
fi

### Create the link to the SCIP directory
echo "Creating link to $SCIPDIR"
if [ -e scip ]; then
   errmsg "ERROR: Directory \"scip\" already exists"
fi
ln -s "$SCIPDIR" scip
if [ $? -ne 0 ]; then
   errmsg "ERROR: Couldn't create link to \"$SCIPDIR\""
fi

### Copy the linear ordering files across
echo "Copying linear ordering files"
if [ ! -e scip/examples/LOP/src/cons_linearordering.c ]; then
   errmsg "ERROR: Couldn't find file \"$SCIPDIR/examples/LOP/src/cons_linearordering.c\""
fi
cp scip/examples/LOP/src/cons_linearordering.c src/
if [ $? -ne 0 ]; then
   errmsg "ERROR: Couldn't copy \"$SCIPDIR/examples/LOP/src/cons_linearordering.c\""
fi
if [ ! -e scip/examples/LOP/src/cons_linearordering.h ]; then
   errmsg "ERROR: Couldn't find file \"$SCIPDIR/examples/LOP/src/cons_linearordering.h\""
fi
cp scip/examples/LOP/src/cons_linearordering.h src/
if [ $? -ne 0 ]; then
   errmsg "ERROR: couldn't copy \"$SCIPDIR/examples/LOP/src/cons_linearordering.h\""
fi

### No problems occurred
echo "SUCCEEDED"
