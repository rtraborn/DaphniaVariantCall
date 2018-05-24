# Updated on 02/08/17
# Use trimmomatic-0.36.jar in /N/dc2/projects/daphpops/Software/Trimmomatic-0.36
# Use fastqutils in /N/dc2/projects/daphpops/Software/ngsutils-ngsutils-0.5.9/bin/
# Use novoalign (V3.02.11) in /N/soft/rhel6/novoalign/novocraft
# Use samtools (Version: 1.3.1) in /N/soft/rhel6/samtools/1.3.1/bin
# Use picard.jar (Version 2.8.1) in /N/soft/rhel6/picard/2.8.1
# Use GenomeAnalysisTK.jar (Version 3.4-0) in /N/soft/rhel6/gatk/3.4-0
# Use bam (Version: 1.0.13) in /N/soft/rhel6/bamUtil/1.0.13

# 0. Make two index files (one is a Novoalign index file and the other is a Samtools index file) and a dictionary file of the reference FASTA file.
# We need to use the same files in the shared directory (/N/dc2/projects/daphpops/PA42_with_mt).  
# This step is written just for showing the entire process.
# Do not make and use different files for our Daphnia project.
# Make a Novoalign index file
/N/soft/rhel6/novoalign/novocraft/novoindex PA42_with_mt.ndx PA42_with_mt.fasta
# Make a Samtools index file
/N/soft/rhel6/samtools/1.3.1/bin/samtools faidx PA42_with_mt.fasta
# Make a dictionary file of a reference file
module load java
# The Java version should not matter as long as it works.
java -jar /N/soft/rhel6/picard/2.8.1/picard.jar CreateSequenceDictionary R=PA42_with_mt.fasta O=PA42_with_mt.dict

# 1. After preparing the FASTA file of adapter sequences, trim adapter sequences from sequence reads.
module load java
java -jar /N/dc2/projects/daphpops/Software/Trimmomatic-0.36/trimmomatic-0.36.jar PE /N/dc2/projects/daphpops/Population_samples/KAP2013/KAP-00030_15lanes_R1.fastq /N/dc2/projects/daphpops/Population_samples/KAP2013/KAP-00030_15lanes_R2.fastq KAP-00030_15lanes_R1-paired.fastq KAP-00030_15lanes_R1-unpaired.fastq KAP-00030_15lanes_R2-paired.fastq KAP-00030_15lanes_R2-unpaired.fastq HEADCROP:3 ILLUMINACLIP:KAP-00030_Adapters.fa:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:30
# The order of the options matters here, as the trimming in Trimmomatic is done in the specified order.
# HEADCROP: Removes the specified number of bases (3 here), regardless of the quality, from the beginning of the read.
# ILLUMINACRIP: Cuts adapter and other Illumina-specific sequences from the read.  The meanings of the parameters are as follows:  
# The first parameter (KAP-00001_Adapters.fa here) specifies the path to a FASTA file containing all adapters.  Make sure to prepare the FASTA file of adapter sequences before running the command.  
# The second parameter (2 here) specifies the maximum mismatch count which will still allow a full match to be performed.
# The third parameter (30 here) specifies how accurate the match between the 'adapter ligated' reads must be for the PE palindrome read alignment.
# The fourth parameter (10 here) specifies how accurate the match between any adapter sequence must be against a read.
# The fifth parameter (2 here) specifies the minimum length of the adapters detected.
# SLIDINGWINDOW: Performs a sliding-window approach.  It starts scanning at the 5' end and clips the read once the average quality within the window falls below a threshold.  The meanings of the parameters are as follows:
# The first parameter (4 here) specifies the number of bases to average across.
# The second parameter (15 here) specifies the average quality required.
# MINLEN: Drops the read if it is below a specified length (30 here).

# 2. Split the FASTQ files into eight pieces.
/N/dc2/projects/daphpops/Software/ngsutils-ngsutils-0.5.9/bin/fastqutils split KAP-00030_15lanes_R1-paired.fastq KAP-00030_15lanes_R1-paired 8 &
/N/dc2/projects/daphpops/Software/ngsutils-ngsutils-0.5.9/bin/fastqutils split KAP-00030_15lanes_R2-paired.fastq KAP-00030_15lanes_R2-paired 8 &
/N/dc2/projects/daphpops/Software/ngsutils-ngsutils-0.5.9/bin/fastqutils split KAP-00030_15lanes_R1-unpaired.fastq KAP-00030_15lanes_R1-unpaired 8 &
/N/dc2/projects/daphpops/Software/ngsutils-ngsutils-0.5.9/bin/fastqutils split KAP-00030_15lanes_R2-unpaired.fastq KAP-00030_15lanes_R2-unpaired 8 &

wait

# 3. Map reads to the reference sequence.
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-paired.1.fastq KAP-00030_15lanes_R2-paired.1.fastq > KAP-00030_PA42_with_mt-paired.1.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-paired.2.fastq KAP-00030_15lanes_R2-paired.2.fastq > KAP-00030_PA42_with_mt-paired.2.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-paired.3.fastq KAP-00030_15lanes_R2-paired.3.fastq > KAP-00030_PA42_with_mt-paired.3.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-paired.4.fastq KAP-00030_15lanes_R2-paired.4.fastq > KAP-00030_PA42_with_mt-paired.4.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-paired.5.fastq KAP-00030_15lanes_R2-paired.5.fastq > KAP-00030_PA42_with_mt-paired.5.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-paired.6.fastq KAP-00030_15lanes_R2-paired.6.fastq > KAP-00030_PA42_with_mt-paired.6.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-paired.7.fastq KAP-00030_15lanes_R2-paired.7.fastq > KAP-00030_PA42_with_mt-paired.7.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-paired.8.fastq KAP-00030_15lanes_R2-paired.8.fastq > KAP-00030_PA42_with_mt-paired.8.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-unpaired.1.fastq > KAP-00030_PA42_with_mt_R1-unpaired.1.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-unpaired.2.fastq > KAP-00030_PA42_with_mt_R1-unpaired.2.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-unpaired.3.fastq > KAP-00030_PA42_with_mt_R1-unpaired.3.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-unpaired.4.fastq > KAP-00030_PA42_with_mt_R1-unpaired.4.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-unpaired.5.fastq > KAP-00030_PA42_with_mt_R1-unpaired.5.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-unpaired.6.fastq > KAP-00030_PA42_with_mt_R1-unpaired.6.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-unpaired.7.fastq > KAP-00030_PA42_with_mt_R1-unpaired.7.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R1-unpaired.8.fastq > KAP-00030_PA42_with_mt_R1-unpaired.8.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R2-unpaired.1.fastq > KAP-00030_PA42_with_mt_R2-unpaired.1.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R2-unpaired.2.fastq > KAP-00030_PA42_with_mt_R2-unpaired.2.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R2-unpaired.3.fastq > KAP-00030_PA42_with_mt_R2-unpaired.3.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R2-unpaired.4.fastq > KAP-00030_PA42_with_mt_R2-unpaired.4.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R2-unpaired.5.fastq > KAP-00030_PA42_with_mt_R2-unpaired.5.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R2-unpaired.6.fastq > KAP-00030_PA42_with_mt_R2-unpaired.6.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R2-unpaired.7.fastq > KAP-00030_PA42_with_mt_R2-unpaired.7.sam &
/N/soft/rhel6/novoalign/novocraft/novoalign -d /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.ndx -r None -o Sam -f KAP-00030_15lanes_R2-unpaired.8.fastq > KAP-00030_PA42_with_mt_R2-unpaired.8.sam &

wait

# 4. Combine the SAM files using Picard.
module load java
java -jar /N/soft/rhel6/picard/2.8.1/picard.jar MergeSamFiles I=KAP-00030_PA42_with_mt-paired.1.sam I=KAP-00030_PA42_with_mt-paired.2.sam I=KAP-00030_PA42_with_mt-paired.3.sam I=KAP-00030_PA42_with_mt-paired.4.sam I=KAP-00030_PA42_with_mt-paired.5.sam I=KAP-00030_PA42_with_mt-paired.6.sam I=KAP-00030_PA42_with_mt-paired.7.sam I=KAP-00030_PA42_with_mt-paired.8.sam I=KAP-00030_PA42_with_mt_R1-unpaired.1.sam I=KAP-00030_PA42_with_mt_R1-unpaired.2.sam I=KAP-00030_PA42_with_mt_R1-unpaired.3.sam I=KAP-00030_PA42_with_mt_R1-unpaired.4.sam I=KAP-00030_PA42_with_mt_R1-unpaired.5.sam I=KAP-00030_PA42_with_mt_R1-unpaired.6.sam I=KAP-00030_PA42_with_mt_R1-unpaired.7.sam I=KAP-00030_PA42_with_mt_R1-unpaired.8.sam I=KAP-00030_PA42_with_mt_R2-unpaired.1.sam I=KAP-00030_PA42_with_mt_R2-unpaired.2.sam I=KAP-00030_PA42_with_mt_R2-unpaired.3.sam I=KAP-00030_PA42_with_mt_R2-unpaired.4.sam I=KAP-00030_PA42_with_mt_R2-unpaired.5.sam I=KAP-00030_PA42_with_mt_R2-unpaired.6.sam I=KAP-00030_PA42_with_mt_R2-unpaired.7.sam I=KAP-00030_PA42_with_mt_R2-unpaired.8.sam O=KAP-00030_PA42_with_mt.sam

# 5. Convert the SAM file to the BAM file.
/N/soft/rhel6/samtools/1.3.1/bin/samtools view -bS KAP-00030_PA42_with_mt.sam > KAP-00030_PA42_with_mt.bam

# 6. Sort the BAM file using Picard.
module load java
java -jar /N/soft/rhel6/picard/2.8.1/picard.jar SortSam INPUT=KAP-00030_PA42_with_mt.bam OUTPUT=Sorted_KAP-00030_PA42_with_mt.bam SORT_ORDER=coordinate

# 7. Add read groups to the sorted BAM file.
module load java
java -jar /N/soft/rhel6/picard/2.8.1/picard.jar AddOrReplaceReadGroups INPUT=Sorted_KAP-00030_PA42_with_mt.bam OUTPUT=RG_Sorted_KAP-00030_PA42_with_mt.bam RGID=Daphnia RGLB=bar RGPL=illumina RGSM=KAP-00030 RGPU=6
# The following five read-group fields are required for using GATK.
# ID: globally unique string identifying the sequencing run.
# LB: an identifier of the library from which the DNA was sequenced.
# PL: the platform used.
# SM: the name associated with the DNA sample in the file.
# PU: the platform unit identifier for the sequencing run.

# 8. Mark duplicate reads.
module load java
java -jar /N/soft/rhel6/picard/2.8.1/picard.jar MarkDuplicates INPUT=RG_Sorted_KAP-00030_PA42_with_mt.bam OUTPUT=dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam METRICS_FILE=KAP-00030_PA42_with_mt_metrics.txt

# 9. Index the BAM file using Picard.
module load java
java -jar /N/soft/rhel6/picard/2.8.1/picard.jar BuildBamIndex INPUT=dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam

# 10. Define intervals to target for the local realignment.
module load java
java -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.fasta -I dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam -o KAP-00030_PA42_with_mt.intervals

# 11. Locally realign reads around indels.
module load java
java -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -T IndelRealigner -R /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.fasta -I dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam -targetIntervals KAP-00030_PA42_with_mt.intervals -o realigned_dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam

# 12. Clip overlapping read pairs.
/N/soft/rhel6/bamUtil/1.0.13/bam clipOverlap --in realigned_dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam --out Clipped_realigned_dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam

# 13. Index the clipped BAM file using Samtools
/N/soft/rhel6/samtools/1.3.1/bin/samtools index Clipped_realigned_dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam

# 14. Make the mpileup file from the BAM file.
/N/soft/rhel6/samtools/1.3.1/bin/samtools mpileup -f /N/dc2/projects/daphpops/PA42_with_mt/PA42_with_mt.fasta Clipped_realigned_dedup_RG_Sorted_KAP-00030_PA42_with_mt.bam > KAP-00030_PA42_with_mt.mpileup

# Remaining issues:
# 1. Some Java in the IU computing system (e.g., java/1.7.0_51 on Mason) fail to create the virtual machine.
# 2. Need to explore how to specify optimal memory options for using Java to avoid the memory issues.
# 3. Should understand the meanings of the read groups better.
# 4. They might change the version of Novoalign in the directory in the future.
# 5. Samtools 1.3.1 in the directory cannot be used on Big Red II.  Samtools in general needs to be used by using "module load" on Big Red II.