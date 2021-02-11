#!/bin/bash
cd $(dirname $0)
# uncomment the next line to use only one thread. This is needed if you want to reproduce ref_optimized_with_seed_xxx.pdb   
export OMP_NUM_THREADS=1
#select the correct path for shared data and executable
export SCUBA_DATAPATH=../../data
exec=../../bin/SCUBA-RunSD
# change the random number seed as you wish to obtain different optimized backbones
seed=23333
for task in 1-1 1-2 2-1 2-2 2-3
do
sed -e '/^#/d; s/seedvalue/'${seed}'/' sd${task}.in.uncommented >sd${task}.in
${exec} sd${task}.in  >${task}.log
done
mv out2-3.pdb optimized_with_seed_${seed}.pdb
