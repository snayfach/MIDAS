#!/bin/bash

REALME=`realpath $0`
MIDAS=`dirname ${REALME}`
EBSROOT=`dirname ${MIDAS}`
export PYTHONPATH=$PYTHONPATH:${MIDAS}:${MIDAS}/smelter
export PATH=$PATH:${MIDAS}/scripts
export MIDAS_DB=${EBSROOT}/IGGdb/v1.0.0

if [ ! -d ${MIDAS_DB} ]; then
   echo "Directory ${MIDAS_DB} not found.  Please provide it, or adjust ${REALME} for your install."
   exit 1
fi

python3 ${MIDAS}/scripts/run_midas.py $*
