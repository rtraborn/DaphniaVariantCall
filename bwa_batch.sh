#!/bin/bash

#PBS -N bwa_test
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=48gb
#PBS -l walltime=1:00:00
#PBS -q debug

WD=/N/dc2/projects/daph_gene/DaphniaVariantCall

cd $WD

time ./bwa_pipeline.sh

exit
