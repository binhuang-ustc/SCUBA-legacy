#!/bin/sh
for task in 1-1 1-2 2-1 2-2 2-3
do
rm -f sd${task}.in out${task}.pdb ${task}.log 
done
rm -f optimized_with_seed_*.pdb
