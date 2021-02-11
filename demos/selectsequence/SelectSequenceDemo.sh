#!/bin/sh
cd $(dirname $0)
#you need to first download the data files needed to run ABACUS2 sequence deign from https://github.com/binhuang-ustc/ABACUS2-data.git
#
export SCUBA_DATAPATH=../../data

# select amino acid sequence with default ABACUS2 parameters
# On command line, -in specify the input backbone structure, -out names the output, 
# -n specify the number of sequences to design
#-log specify to store  the output sequences and  ABACUS2 energies
#../../bin/ABACUS-DesignSeq -in target.pdb -out seqdesign.pdb -log design.seq -n 1 >designseq.log


#The following  commands elext sequence with other parameters
#reference energy of alanine at non-buried helix position increased by 0.5
../../bin/ABACUS-DesignSeq -in target.pdb -out seqdesign.pdb -log design.seq -n 1 -para ./abacus.param_ref1 >designseq.log


#reference energy of alanine at non-buried helix position increased by 0.5, of Ile and Val at non-buried strand position increased by 0.9, and Thr at strand positions increased by 1.2
#../../bin/ABACUS-DesignSeq -in target.pdb -out seqdesign.pdb -log design.seq -n 1 -para ./abacus.param_ref2 >designseq.log


#Reference energies of polar residues (except Thr) at exposed positions on secondary structure elements decreased by 0.8. Ile i(Val) at non-buried strand positions increased by 1.0(1.5). Thr at non-buried strand positions increasde by 1.2. 
#../../bin/ABACUS-DesignSeq -in target.pdb -out seqdesign.pdb -log design.seq -n 1 -para ./abacus.param_ref3 >designseq.log

# shrink the Van der Waals radii of side chain atoms when calculate clash interaction
#../../bin/ABACUS-DesignSeq -in target.pdb -out seqdesign.pdb -log design.seq -n 1 -para ./abacus.param_shrinkSC >designseq.log


