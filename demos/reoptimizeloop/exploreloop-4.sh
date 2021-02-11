#!/bin/bash

cd $(dirname $0)

export SCUBA_DATAPATH=../../data
export OMP_NUM_THREADS=1
exe=../../bin/exploreloops

$exe ./exploreloop-4.in >exploreloop-4.log
