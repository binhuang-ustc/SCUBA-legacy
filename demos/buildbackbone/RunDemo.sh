#!/bin/bash

cd $(dirname $0)

export SCUBA_DATAPATH=../../data

# build initial backbone according to sketch definition file
../../bin/buildbackbone sketch_definition.txt initial.pdb 11329 >/dev/null
