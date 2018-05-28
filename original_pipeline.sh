#!/usr/bin/bash

######### Local paths
WD=/N/u/rtraborn/Carbonate/scratch/DaphniaVariantCall

######### Loading carbonate-specific resources (change if not using this cluster) #######
module load java
#module load trimmomatic #this loads version 0.36
module load fastqc

######### Paths to binaries #########
novoindex=/N/soft/rhel6/novoalign/novocraft/novoindex
novoalign=/N/soft/rhel6/novoalign/novocraft/novoalign
samtools=/N/soft/rhel6/samtools/1.3.1/bin/samtools
picard='java -jar /N/soft/rhel6/picard/2.8.1/picard.jar'
Trimmomatic='java -jar /N/dc2/projects/daphpops/Software/Trimmomatic-0.36/trimmomatic-0.36.jar'
GATK='java -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar'
bamUtil=/N/soft/rhel6/bamUtil/1.0.13/bam

# please install ngsutils (see INSTALL.txt for instructions)
# then, provide the path to the fastqutils binary 
fastqutils=/N/u/rtraborn/Carbonate/scratch/DaphniaVariantCall/software/ngsutils/bin/fastqutils

####### Path to adapters #######
adapterTrim=../adapters/Bioo_Adapters.fa

###### Misc names ######
assemblyName=PA42_with_mt.fasta
novoIndexName=PA42_with_mt.ndx
sequenceDict=PA42_with_mt.dict

##### Sequence names #####
clone1_R1=KAP-00074_CGCTGATC_L008_R2_001.fastq
clone1_R2=KAP-00074_CGCTGATC_L008_R2_001.fastq

cd $WD

echo "Copying sample files (from clone KAP-00074) to a new directory called fastq/."
mkdir fastq
cp /N/dc2/projects/daphpops/Population_samples/KAP2013/140501/Sample_KAP-00074/*.fastq fastq

echo "Creating a symbolic link to the PA42 assembly"
mkdir assembly
cd assembly
ln -s /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.fasta $assemblyName

echo "Making a Novoalign index file."
$novoindex $novoIndexName $assemblyName

echo "Making a samtools index file"
$samtools faidx $assemblyName

echo "Making a dictionary file of a reference."
# The Java version should not matter as long as it works. (ed: ?)
$picard CreateSequenceDictionary R=$assemblyName O=$sequenceDict

# 0. Performing fastqc on sample

echo "Performing fastqc on read pairs for this sample."
cd ../fastq
fastqc $clone_R1 $clone_R2 

# 1. After preparing the FASTA file of adapter sequences, trim adapter sequences from sequence reads.

echo "Trimming adapter sequences from sequence reads."
$Trimmomatic PE $clone1_R1 $clone1_R2 KAP-00074_15lanes_R1-paired.fastq KAP-00074_15lanes_R1-unpaired.fastq KAP-00074_15lanes_R2-paired.fastq KAP-00074_15lanes_R2-unpaired.fastq HEADCROP:3 ILLUMINACLIP:$adapterTrim:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:30

# The additional information on Trimmomatic arguments was omitted because it is available in the Trimmomatic documentation.

echo "Splitting the fastq files into 8 parts. We are doing this because we are using the free version of Novoalign, which does not permit multi-threading."
$fastqutils split KAP-00074_15lanes_R1-paired.fastq KAP-00074_15lanes_R1-paired 8 &
$fastqutils split KAP-00074_15lanes_R2-paired.fastq KAP-00074_15lanes_R2-paired 8 &
$fastqutils split KAP-00074_15lanes_R1-unpaired.fastq KAP-00074_15lanes_R1-unpaired 8 &
$fastqutils split KAP-00074_15lanes_R2-unpaired.fastq KAP-00074_15lanes_R2-unpaired 8 &

wait

# 3. Map reads to the reference sequence.

$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f KAP-00074_15lanes_R1-paired.1.fastq KAP-00074_15lanes_R2-paired.1.fastq > KAP-00074_PA42_with_mt-paired.1.sam &
$novoalign -d ../assembly/$novoIndexName  -r None -o Sam -f KAP-00074_15lanes_R1-paired.2.fastq KAP-00074_15lanes_R2-paired.2.fastq > KAP-00074_PA42_with_mt-paired.2.sam &
$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f KAP-00074_15lanes_R1-paired.3.fastq KAP-00074_15lanes_R2-paired.3.fastq > KAP-00074_PA42_with_mt-paired.3.sam &
$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f KAP-00074_15lanes_R1-paired.4.fastq KAP-00074_15lanes_R2-paired.4.fastq > KAP-00074_PA42_with_mt-paired.4.sam &
$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f KAP-00074_15lanes_R1-paired.5.fastq KAP-00074_15lanes_R2-paired.5.fastq > KAP-00074_PA42_with_mt-paired.5.sam &
$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f KAP-00074_15lanes_R1-paired.6.fastq KAP-00074_15lanes_R2-paired.6.fastq > KAP-00074_PA42_with_mt-paired.6.sam &
$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f KAP-00074_15lanes_R1-paired.7.fastq KAP-00074_15lanes_R2-paired.7.fastq > KAP-00074_PA42_with_mt-paired.7.sam &
$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f KAP-00074_15lanes_R1-paired.8.fastq KAP-00074_15lanes_R2-paired.8.fastq > KAP-00074_PA42_with_mt-paired.8.sam &
wait

# 4. Combine the SAM files using Picard.

$picard MergeSamFiles I=KAP-00074_PA42_with_mt-paired.1.sam I=KAP-00074_PA42_with_mt-paired.2.sam I=KAP-00074_PA42_with_mt-paired.3.sam I=KAP-00074_PA42_with_mt-paired.4.sam I=KAP-00074_PA42_with_mt-paired.5.sam I=KAP-00074_PA42_with_mt-paired.6.sam I=KAP-00074_PA42_with_mt-paired.7.sam I=KAP-00074_PA42_with_mt-paired.8.sam O=KAP-00074_PA42_with_mt.sam

# 5. Convert the SAM file to the BAM file.
$samtools view -bS KAP-00074_PA42_with_mt.sam > KAP-00074_PA42_with_mt.bam

# 6. Sort the BAM file using Picard.
echo "Sorting the bam file using Picard."
$picard SortSam INPUT=KAP-00074_PA42_with_mt.bam OUTPUT=Sorted_KAP-00074_PA42_with_mt.bam SORT_ORDER=coordinate

# 7. Add read groups to the sorted BAM file.
echo "Adding read groups to the sorted bam file."
$picard AddOrReplaceReadGroups INPUT=Sorted_KAP-00074_PA42_with_mt.bam OUTPUT=RG_Sorted_KAP-00074_PA42_with_mt.bam RGID=Daphnia RGLB=bar RGPL=illumina RGSM=KAP-00074 RGPU=6

# 8. Mark duplicate reads.
echo "Marking duplicates using Picard."
$picard MarkDuplicates INPUT=RG_Sorted_KAP-00074_PA42_with_mt.bam OUTPUT=dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam METRICS_FILE=KAP-00074_PA42_with_mt_metrics.txt

# 9. Index the BAM file using Picard.
echo "Indexing the bam file using Picard."
$picard BuildBamIndex INPUT=dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam

# 10. Define intervals to target for the local realignment.
echo "Defining intervals to target for local realignment using Picard."
$GATK -T RealignerTargetCreator -R ../assembly/$assemblyName -I dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam -o KAP-00074_PA42_with_mt.intervals

# 11. Locally realign reads around indels.
echo "Performing the local realignment using Picard."
$GATK -T IndelRealigner -R ../assembly/$assemblyName -I dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam -targetIntervals KAP-00074_PA42_with_mt.intervals -o realigned_dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam

# 12. Clip overlapping read pairs.
echo "Clipping the overlapping read pairs using bamUtil."
$bamUtil clipOverlap --in realigned_dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam --out Clipped_realigned_dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam

# 13. Index the clipped BAM file using Samtools
echo "Indexing the clipped BAM file using Samtools."
$samtools index Clipped_realigned_dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam

# 14. Make the mpileup file from the BAM file.
echo "Creating the mpileup file from the BAM file."
$samtools mpileup -f ../assembly/$assemblyName Clipped_realigned_dedup_RG_Sorted_KAP-00074_PA42_with_mt.bam > KAP-00074_PA42_with_mt.mpileup

# Remaining issues:
# 1. Some Java in the IU computing system (e.g., java/1.7.0_51 on Mason) fail to create the virtual machine.
# 2. Need to explore how to specify optimal memory options for using Java to avoid the memory issues.
# 3. Should understand the meanings of the read groups better.
# 4. They might change the version of Novoalign in the directory in the future.
# 5. Samtools 1.3.1 in the directory cannot be used on Big Red II.  Samtools in general needs to be used by using "module load" on Big Red II.
