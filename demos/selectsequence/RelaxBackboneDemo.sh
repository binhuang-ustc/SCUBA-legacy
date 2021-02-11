#!/bin/sh
cd $(dirname $0)
export SCUBA_DATAPATH=../../data
export OMP_NUM_THREADS=1
#First, 2000 steps of SD simulation with fircition coeffcient =5 and Tr=0.2
# the resulting structure in in relax1.pdb
../../bin/SCUBA-RunSD sdrelax1.in >sdrelax1.log
#then, 5000 steps of SD simulation with fircition coeffcient =1 and Tr=0.5
# the resulting structure in in relax2.pdb
../../bin/SCUBA-RunSD sdrelax2.in >sdrelax2.log
