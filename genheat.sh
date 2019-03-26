#!/bin/sh
if [ -z "$1" ]
then
    echo "Usage: genheat <infile>"
else
    mpirun -np 8 ./cv "$1"
    python3 cv.py "$1.out" 1024
fi
