#!/bin/bash
cd $(dirname $0)
#export OMP_NUM_THREADS=1
#select correct path fr=or shared data and executable
export SCUBA_DATAPATH=../../data
exec=../../bin/SCUBA-RunSD
sed -e '/^#/d' sdwithrmsdrestraint.in.uncommented >sdwithrmsdrestraint.in
${exec} sdwithrmsdrestraint.in  >sdrmsdrestrained.log
