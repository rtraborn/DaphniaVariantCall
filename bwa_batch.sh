#!/bin/bash

#PBS -N bwa_test
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=48gb
#PBS -l walltime=1:00:00
#PBS -q debug

WD=/N/u/rtraborn/Carbonate/scratch/development/DaphniaVariantCall

cd $WD

time ./bwa_pipeline.sh

exit
