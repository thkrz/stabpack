#!/bin/sh -e
export PYTHONPATH="/home/thk32is/work/src/ssat:$PYTHONPATH"
cd examples && ../bin/ssat ssat.cfg "$1"
