#!/usr/bin/bash

######### Local paths
WD=/N/u/rtraborn/Carbonate/scratch/DaphniaVariantCall

######### Loading carbonate-specific resources (change if not using this cluster) #######
module load java
#module load trimmomatic #this loads version 0.36
module load fastqc
module load hisat2
module load samtools # this will load Samtools v 1.5 on Carbonate

######### Paths to binaries #########
picard='java -jar /N/soft/rhel6/picard/2.8.1/picard.jar'
Trimmomatic='java -jar /N/dc2/projects/daphpops/Software/Trimmomatic-0.36/trimmomatic-0.36.jar'
GATK='java -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar'
bamUtil=/N/soft/rhel6/bamUtil/1.0.13/bam

######### Paths to assembly ######
assemblyDir=/N/u/rtraborn/Carbonate/scratch/projects/PA42_4_0
assemblyName=PA42.4.0.fasta

######### Paths to reads #########
fastqBase=/N/dc2/projects/daphpops/Population_samples/KAP2013/140501
SampleDir_1=Sample_KAP-00074
CloneID_1=KAP-00074
SampleDir_2=Sample_KAP-00075
CloneID_2=KAP-00075
SampleDir_3=Sample_KAP-00079
CloneID_3=KAP-00079

# please install ngsutils (see INSTALL.txt for instructions)
# then, provide the path to the fastqutils binary 
# please note that fastqutils is only required for the original pipeline (i.e. the one using novoalign)
# fastqutils=/N/u/rtraborn/Carbonate/scratch/DaphniaVariantCall/software/ngsutils/bin/fastqutils

####### Path to adapters #######
adapterTrim=../adapters/Bioo_Adapters.fa

###### Misc names ######
assemblyID=PA42_4_0
sequenceDict=$(basename $assemblyName .fasta).dict

##### Number of threads #####
nThreads=8

cd $WD

echo "Creating symbolic links to files (from clone ${CloneID_1}) to a new directory called fastq/."
mkdir fastq
cd fastq
   for fq in `ls $fastqBase/$SampleDir_1/*.fastq`; do
       ln -s $fq $(basename $fq) 
   done

echo "Creating symbolic links to files (from clone ${CloneID_2}) to a new directory called fastq/."
   for fq in `ls $fastqBase/$SampleDir_2/*.fastq`; do
       ln -s $fq $(basename $fq) 
   done

echo "Creating symbolic links to files (from clone ${CloneID_3}) to a new directory called fastq/."
   for fq in `ls $fastqBase/$SampleDir_3/*.fastq`; do
       ln -s $fq $(basename $fq) 
   done

echo "Creating a single, combined fastq file (from clones ${CloneID_1}, ${CloneID_2} and ${CloneID_3} from the R1 and R2 fastq reads, respectively."

Clone1_R1_fqs=`ls ${CloneID_1}_*_R1_001.fastq`
Clone1_R2_fqs=`ls ${CloneID_1}_*_R2_001.fastq`
Clone2_R1_fqs=`ls ${CloneID_2}_*_R1_001.fastq`
Clone2_R2_fqs=`ls ${CloneID_2}_*_R2_001.fastq`
Clone3_R1_fqs=`ls ${CloneID_3}_*_R1_001.fastq`
Clone3_R2_fqs=`ls ${CloneID_3}_*_R2_001.fastq`
cat $Clone1_R1_fqs > ${CloneID_1}_merged_R1.fastq
cat $Clone1_R2_fqs > ${CloneID_1}_merged_R2.fastq
cat $Clone2_R1_fqs > ${CloneID_2}_merged_R1.fastq
cat $Clone2_R2_fqs > ${CloneID_2}_merged_R2.fastq
cat $Clone3_R1_fqs > ${CloneID_3}_merged_R1.fastq
cat $Clone3_R2_fqs > ${CloneID_3}_merged_R2.fastq

cd ..

echo "Creating a symbolic link to the PA42 assembly"
mkdir assembly
cd assembly
ln -s ${assemblyDir}/${assemblyName} $assemblyName

echo "Making a samtools index file"
samtools faidx $assemblyName

echo "Making a dictionary file of a reference."
# The Java version should not matter as long as it works. (ed: ?)
$picard CreateSequenceDictionary R=$assemblyName O=$sequenceDict

# 0. Performing fastqc on sample (uncomment the following lines to run this step)

#echo "Performing fastqc on read pairs for this sample."
#cd ../fastq
#fastqc $clone_R1 $clone_R2 
#cd ..

# 1. After preparing the FASTA file of adapter sequences, trim adapter sequences from sequence reads.

echo "Trimming adapter sequences from sequence reads."
cd ../fastq

for CloneID in KAP-00074 KAP-00075 KAP-00079; do

echo $CloneID

$Trimmomatic PE ${CloneID}_merged_R1.fastq ${CloneID}_merged_R2.fastq ${CloneID}_R1-paired.fastq ${CloneID}_R1-unpaired.fastq ${CloneID}_R2-paired.fastq ${CloneID}_R2-unpaired.fastq HEADCROP:3 ILLUMINACLIP:$adapterTrim:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:30

done

cd ..

# The additional information originally provided on Trimmomatic arguments was omitted because it is available in the Trimmomatic documentation.

# 2. Building the Hisat2 index using the assembly

cd assembly
hisat2-build $assemblyName $assemblyID

# 3. Map reads to the reference sequence using Hisat2.
echo "Mapping reads to the reference genome using Hisat2."
for CloneID in KAP-00074 KAP-00075 KAP-00079; do

if [ -s ../fastq/${CloneID}_R1-unpaired.fastq ] || [ -s ../fastq/${CloneID}_R2-unpaired.fastq ];
then echo "One of the unpaired read files are empty"
hisat2 --no-spliced-alignment -p $nThreads -q -x $assemblyID -1 ../fastq/${CloneID}_R1-paired.fastq -2 ../fastq/${CloneID}_R2-paired.fastq -S ../fastq/${CloneID}_${assemblyID}.sam
else echo "All of the read files are fine."
hisat2 --no-spliced-alignment -p $nThreads -q -x $assemblyID -1 ../fastq/${CloneID}_R1-paired.fastq -2 ../fastq/${CloneID}_R2-paired.fastq -U ../fastq/${CloneID}_R1-unpaired.fastq ../fastq/${CloneID}_R2-unpaired.fastq  -S ../fastq/${CloneID}_${assemblyID}.sam
fi

# 4. Convert the SAM file to the BAM file.
echo "Converting the sam file to bam and removing non-primary alignments."
cd ../fastq

samtools view -b -F 256 ${CloneID}_${assemblyID}.sam > ${CloneID}_${assemblyID}.bam 

# 5. Sort the BAM file using Picard.
echo "Sorting the bam file using Picard."
$picard SortSam INPUT=${CloneID}_${assemblyID}.bam OUTPUT=${CloneID}_${assemblyID}_sorted.bam SORT_ORDER=coordinate

# 6. Add read groups to the sorted BAM file.
echo "Adding read groups to the sorted bam file."
$picard AddOrReplaceReadGroups INPUT=${CloneID}_${assemblyID}_sorted.bam OUTPUT=${CloneID}_${assemblyID}_sorted_rg.bam RGID=Daphnia RGLB=bar RGPL=illumina RGSM=$CloneID RGPU=6

# 7. Mark duplicate reads.
echo "Marking duplicates using Picard."
$picard MarkDuplicates INPUT=${CloneID}_${assemblyID}_sorted_rg.bam OUTPUT=${CloneID}_${assemblyID}_sorted_rg_dedup.bam METRICS_FILE=${CloneID}_${assemblyID}_metrics.txt

# 8. Index the BAM file using Picard.
echo "Indexing the bam file using Picard."
$picard BuildBamIndex INPUT=${CloneID}_${assemblyID}_sorted_rg_dedup.bam

# 9. Define intervals to target for the local realignment.
echo "Defining intervals to target for local realignment using Picard."
$GATK -T RealignerTargetCreator -R ${WD}/assembly/${assemblyName} -I ${CloneID}_${assemblyID}_sorted_rg_dedup.bam -o ${CloneID}_${assemblyID}.intervals

# 10. Locally realign reads around indels.
echo "Performing the local realignment using Picard."
$GATK -T IndelRealigner -R ${WD}/assembly/${assemblyName} -I ${CloneID}_${assemblyID}_sorted_rg_dedup.bam -targetIntervals ${CloneID}_${assemblyID}.intervals -o ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned.bam

# 11. Clip overlapping read pairs.
echo "Clipping the overlapping read pairs using bamUtil."
$bamUtil clipOverlap --in ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned.bam --out ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned_clip.bam 

# 12. Index the clipped BAM file using Samtools
echo "Indexing the clipped BAM file using Samtools."
samtools index ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned_clip.bam

# 13. Make the mpileup file from the BAM file.
echo "Creating the mpileup file from the BAM file."
samtools mpileup -f ../assembly/$assemblyName ${CloneID}_${assemblyID}_sorted_rg_dedup_realigned_clip.bam -o ${CloneID}_${assemblyID}.mpileup

done

echo "Job complete!"
