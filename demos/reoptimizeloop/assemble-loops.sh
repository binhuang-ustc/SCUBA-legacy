#!/bin/bash

cd $(dirname $0)

export SCUBA_DATAPATH=../../data

outputdir="./resampled-confs"
rm -rf $outputdir
mkdir $outputdir

numconfs=10
pdbprefix="conf-"
#in the command line parameters, 10 means 10 lowest energy nonredundant conformations
# for each loop will be considered. 
#The 2 means loop configuations with length 2 residues shorter thant the optimum
# loop length will be considered if their energies are lower than the 10th lowest
#energy of the optimum loop length
#The 0.1 means the criterion for nonredundancy is RMSD larger #than 0.1*looplength
../../bin/processloopexplorationresult assemble 10 2 0.1 $numconfs $pdbprefix $outputdir sampledloops-1.dat sampledloops-2.dat sampledloops-3.dat sampledloops-4.dat sampledloops-5.dat
