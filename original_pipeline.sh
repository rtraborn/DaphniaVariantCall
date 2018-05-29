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

######### Paths to reads #########
fastqBase=/N/dc2/projects/daphpops/Population_samples/KAP2013/140501
SampleDir=Sample_KAP-00074
CloneID=KAP-00074

# please install ngsutils (see INSTALL.txt for instructions)
# then, provide the path to the fastqutils binary 
fastqutils=/N/u/rtraborn/Carbonate/scratch/DaphniaVariantCall/software/ngsutils/bin/fastqutils

####### Path to adapters #######
adapterTrim=../adapters/Bioo_Adapters.fa

###### Misc names ######
assemblyName=PA42_with_mt.fasta
novoIndexName=PA42_with_mt.ndx
sequenceDict=PA42_with_mt.dict

cd $WD

echo "Creating symbolic links to files (from clone ${CloneID}) to a new directory called fastq/."
mkdir fastq
cd fastq
   for fq in `ls $fastqBase/$SampleDir/*.fastq`; do
       ln -s $fq $(basename $fq) 
   done

echo "Creating a single, combined fastq file (from clone ${CloneID} from the R1 and R2 fastq reads, respectively."

R1_fqs=`ls *_R1_001.fastq`
R2_fqs=`ls *_R2_001.fastq`
cat $R1_fqs > ${CloneID}_merged_R1.fastq
cat $R2_fqs > ${CloneID}_merged_R2.fastq

cd ..

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

#echo "Performing fastqc on read pairs for this sample."
#cd ../fastq
#fastqc $clone_R1 $clone_R2 

# 1. After preparing the FASTA file of adapter sequences, trim adapter sequences from sequence reads.

echo "Trimming adapter sequences from sequence reads."
cd ../fastq
$Trimmomatic PE ${CloneID}_merged_R1.fastq ${CloneID}_merged_R2.fastq ${CloneID}_R1-paired.fastq ${CloneID}_R1-unpaired.fastq ${CloneID}_R2-paired.fastq ${CloneID}_R2-unpaired.fastq HEADCROP:3 ILLUMINACLIP:$adapterTrim:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:30
cd ..

# The additional information on Trimmomatic arguments was omitted because it is available in the Trimmomatic documentation.

echo "Splitting the fastq files into 8 parts. We are doing this because we are using the free version of Novoalign, which does not permit multi-threading."
$fastqutils split ${CloneID}_R1-paired.fastq ${CloneID}_R1-paired.fastq 8 &
$fastqutils split ${CloneID}_R2-paired.fastq ${CloneID}_R2-paired.fastq 8 &
#$fastqutils split ${CloneID}_merged_R1_unpaired.fastq ${CloneID}_merged_R1_unpaired.fastq 8 &
#$fastqutils split ${CloneID}_merged_R2_unpaired.fastq ${CloneID}_merged_R2_unpaired.fastq 8 &

wait

# 3. Map reads to the reference sequence.

for i in {1..8}; do
$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f ${CloneID}_R1-paired.${i}.fastq ${CloneID}_R2-paired.${i}.fastq > ${CloneID}_${assemblyID}.${i}.sam &
#$novoalign -d ../assembly/$novoIndexName -r None -o Sam -f ${CloneID}_R1-unpaired.${i}.fastq ${CloneID}_R2-unpaired.${i}.fastq > ${CloneID}_${assemblyID}_unpaired.${i}.sam &
wait

# 4. Combine the SAM files using Picard.
echo "Combining the SAM files using Picard."

$picard MergeSamFiles  I=${CloneID}_${assemblyID}.1.sam I=${CloneID}_${assemblyID}.2.sam I=${CloneID}_${assemblyID}.3.sam I=${CloneID}_${assemblyID}.4.sam ${CloneID}_${assemblyID}.5.sam ${CloneID}_${assemblyID}.6.sam ${CloneID}_${assemblyID}.7.sam ${CloneID}_${assemblyID}.8.sam O=${CloneID}_${assemblyID}.sam

# 5. Convert the SAM file to the BAM file.
echo "Converting the file from SAM to BAM format."
$samtools view -bS ${CloneID}_${assemblyID}.sam > ${CloneID}_${assemblyID}.bam

# 6. Sort the BAM file using Picard.
echo "Sorting the bam file using Picard."
$picard SortSam INPUT=${CloneID}_${assemblyID}.bam OUTPUT=${CloneID}_${assemblyID}_sorted.bam SORT_ORDER=coordinate

# 7. Add read groups to the sorted BAM file.
echo "Adding read groups to the sorted bam file."
$picard AddOrReplaceReadGroups INPUT=${CloneID}_${assemblyID}_sorted.bam OUTPUT=${CloneID}_${assemblyID}_sorted_rg.bam RGID=Daphnia RGLB=bar RGPL=illumina RGSM=$CloneID RGPU=6

# 8. Mark duplicate reads.
echo "Marking duplicates using Picard."
$picard MarkDuplicates INPUT=${CloneID}_${assemblyID}_sorted_rg.bam OUTPUT=${CloneID}_${assemblyID}_sorted_rg_dedup.bam METRICS_FILE=${CloneID}_${assemblyID}_metrics.txt

# 9. Index the BAM file using Picard.
echo "Indexing the bam file using Picard."
$picard BuildBamIndex INPUT=${CloneID}_${assemblyID}_sorted_rg_dedup.bam

# 10. Define intervals to target for the local realignment.
echo "Defining intervals to target for local realignment using Picard."
$GATK -T RealignerTargetCreator -R ../assembly/$assemblyName -I ${CloneID}_${assemblyID}_sorted_rg_dedup.bam -o ${CloneID}_${assemblyID}.intervals

# 11. Locally realign reads around indels.
echo "Performing the local realignment using Picard."
$GATK -T IndelRealigner -R ../assembly/$assemblyName -I ${CloneID}_${assemblyID}_sorted_rg_dedup.bam -targetIntervals ${CloneID}_${assemblyID}.intervals -o ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned.bam

# 12. Clip overlapping read pairs.
echo "Clipping the overlapping read pairs using bamUtil."
$bamUtil clipOverlap --in ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned.bam --out ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned_clip.bam 

# 13. Index the clipped BAM file using Samtools
echo "Indexing the clipped BAM file using Samtools."
$samtools index ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned_clip.bam

# 14. Make the mpileup file from the BAM file.
echo "Creating the mpileup file from the BAM file."
$samtools mpileup -f ../assembly/$assemblyName ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned_clip.bam -o ${CloneID}_${assemblyID}.mpileup

# Remaining issues (Ed: these are notes from from the original pipeline.)
# 1. Some Java in the IU computing system (e.g., java/1.7.0_51 on Mason) fail to create the virtual machine.
# 2. Need to explore how to specify optimal memory options for using Java to avoid the memory issues.
# 3. Should understand the meanings of the read groups better.
# 4. They might change the version of Novoalign in the directory in the future.
# 5. Samtools 1.3.1 in the directory cannot be used on Big Red II.  Samtools in general needs to be used by using "module load" on Big Red II.
