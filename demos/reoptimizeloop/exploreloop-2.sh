#!/bin/bash

cd $(dirname $0)

export SCUBA_DATAPATH=../../data
export OMP_NUM_THREADS=1
exe=../../bin/exploreloops

$exe ./exploreloop-2.in >exploreloop-2.log
