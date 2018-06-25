#!/bin/bash

#PBS -N hisat2_test
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=48gb
#PBS -l walltime=1:00:00
#PBS -q debug

WD=/N/u/rtraborn/Carbonate/scratch/DaphniaVariantCall/hisat2_align_test

cd $WD

time ./hisat2_pipeline_test.sh

exit
