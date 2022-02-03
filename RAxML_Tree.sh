#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=4:mem=12gb
#PBS -N NF_Tree_Coordinator
#PBS -j oe

## A bit redundant, but just to ensure the files are in the correct location. 
Project_Dir="/path/to/project/directory"
cd $Project_Dir

module load nextflow/20.10.0

echo "Starting: `date`"
nextflow run RAxML_Tree.nf -c RAxML_Tree.config --profile imperial
echo "Finished: `date`"
