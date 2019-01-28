#!/bin/bash

set -e

if [ $# -eq 0 ]; then
    ARGS=*.py
else
    ARGS=$*
fi

exec pylint3 --init-hook='import sys; sys.path.append("smelter")' ${ARGS}



