#!/bin/bash

#PBS -N original_test
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=48gb
#PBS -l walltime=1:00:00
#PBS -q debug

WD=/N/u/rtraborn/Carbonate/scratch/DaphniaVariantCall

cd $WD

time ./original_pipeline_novoalign.sh

exit
