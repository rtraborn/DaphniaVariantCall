#! /bin/usr/perl -w

#The reference genome must be first indexed using the following command:
#bwa index $ref_genome.fasta 

#$ref_genome="/PATH/TO/Reference/Genome";
$ref_genome="PA42.4.0";

# The adapter file: an example (Bioo_Adapters.fa) can be found in the same directory
#$Adapters="/PATH/TO/Adapters.fa";
$Adapters="/N/u/xw63/Carbonate/daphnia/Bioo_Adapters.fa";

#Save all your raw reads in the DATA_DIR in a sub dir named as SampleID
#Name you files like: 
#SampleID-001-R1.fastq
#SampleID-001-R2.fastq
#SampleID-002-R1.fastq
#SampleID-002-R2.fastq
#......
#SampleID-100-R1.fastq
#SampleID-100-R2.fastq

$SampleID="PA2014";  
$DATA_DIR="/N/dc2/scratch/xw63/".$SampleID;
$MaxNumberofSamples=100;
$emailaddress='ouqd@hotmail.com';

# The paths to the software used in this pipeline
# You must first make sure you have these software installed and they are all functional
 
$BWA="~/bwa-0.7.17/bwa";
$novoalign="/N/soft/rhel6/novoalign/novocraft/novoalign";
$novoindex="/N/u/xw63/Carbonate/daphnia/$ref_genome.ndx";
$PICARD="/N/soft/rhel6/picard/2.8.1/picard.jar";
$samtools="/N/soft/rhel6/samtools/1.3.1/bin/samtools";
$Trimmomatic="~/Trimmomatic-0.36/trimmomatic-0.36.jar";
$fastqutils="/N/dc2/projects/daphpops/Software/ngsutils-ngsutils-0.5.9/bin/fastqutils";
$GATK="/N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar";

#Now we find the raw reads and produce pbs files for them
$n=0;
$n1=0;
while ($n<=$MaxNumberofSamples) {
	$n=$n+1;
	$nstr=chr($n);
	$nlen=length($nstr);
	$nstr001= sprintf ("%03d", $n);
	print $nstr001.": ";
		$Sample=$DATA_DIR."/".$SampleID."-".$nstr001;
		$Sample_R1=$DATA_DIR."/fastq/".$SampleID."-".$nstr001."-R1.fastq";
		$Sample_R2=$DATA_DIR."/fastq/".$SampleID."-".$nstr001."-R2.fastq";
		$OUTPUT_DIR=$DATA_DIR."/Bwa";
		$OUTPUT=$OUTPUT_DIR."/".$SampleID."-".$nstr001;
	print $Sample_R1."/R2.fastq";
	if(-e $Sample_R1 && -e $Sample_R2){  
		print ", Okay, this pair-end reads fastq file is found! lets make a pbs file:";  
		$n1=$n1+1;	
		
		$pbsfile=$DATA_DIR."/bwa-".$SampleID."-".$nstr001.".pbs";
		print $pbsfile."\n";
			
		open OUT, ">$pbsfile" or die "cannot open file: $!";
		print OUT 
		"	#!/bin/bash	
			#PBS -N Bwa-$SampleID-$nstr001
			#PBS -l nodes=1:ppn=8
			#PBS -l vmem=100gb
			#PBS -l walltime=12:00:00
			#PBS -M $emailaddress
			#PBS -m abe
			#PBS -j oe

			#module load samtools
			#module load java
			
			mkdir $OUTPUT_DIR
			
			ulimit -s
			
			#Notice: The following line should be executed before submitting the jobs
			#bwa index $ref_genome.fasta 

			# 1. After preparing the FASTA file of adapter sequences, trim adapter sequences from sequence reads.
			module load java
			java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=~/tmp -jar $Trimmomatic PE $Sample_R1 $Sample_R2 $Sample_R1-paired.fq $Sample_R1-unpaired.fq $Sample_R2-paired.fq $Sample_R2-unpaired.fq  HEADCROP:3 ILLUMINACLIP:$Adapters:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:30

			# 3. Map reads to the reference sequence and output bam file.

			$BWA mem -t 8 -M -k 30 $ref_genome.fasta $Sample_R1-paired.fq $Sample_R2-paired.fq > $OUTPUT-paired.sam &
			$BWA mem -t 8 -M -k 30 $ref_genome.fasta $Sample_R1-unpaired.fq  > $OUTPUT-R1-unpaired.sam & 
			$BWA mem -t 8 -M -k 30 $ref_genome.fasta $Sample_R2-unpaired.fq  > $OUTPUT-R2-unpaired.sam &

			wait 

			# 4. Combine the SAM files using Picard.
			module load java
			java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=~/tmp -jar $PICARD MergeSamFiles I=$OUTPUT-paired.sam I=$OUTPUT-R1-unpaired.sam I=$OUTPUT-R2-unpaired.sam O=$OUTPUT.sam

			# 5. Convert the SAM file to the BAM file.
			$samtools view -bS $OUTPUT.sam > $OUTPUT.bam

			# 6. Sort the BAM file using Picard.
			module load java
			java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=/N/dc2/scratch/xw63/tmp -jar $PICARD SortSam INPUT=$OUTPUT.bam OUTPUT=$OUTPUT-Sorted.bam SORT_ORDER=coordinate

			# 7. Add read groups to the sorted BAM file.
			module load java
			java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=/N/dc2/scratch/xw63/tmp -jar $PICARD AddOrReplaceReadGroups INPUT=$OUTPUT-Sorted.bam OUTPUT=$OUTPUT-RG_Sorted.bam RGID=Daphnia RGLB=bar RGPL=illumina RGSM=$Sample RGPU=6

			# 8. Mark duplicate reads.
			module load java
			java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=/N/dc2/scratch/xw63/tmp -jar $PICARD MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 MAX_RECORDS_IN_RAM=5000000 VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true INPUT=$OUTPUT-RG_Sorted.bam OUTPUT=$OUTPUT-RG_Sorted_dedup.bam METRICS_FILE=$OUTPUT-metrics.txt

			# 9. Index the BAM file using Picard.
			module load java
			java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=/N/dc2/scratch/xw63/tmp -jar $PICARD BuildBamIndex INPUT=$OUTPUT-RG_Sorted_dedup.bam

			# 10. Define intervals to target for the local realignment.
			module load java
			java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=~/tmp -jar $GATK -T RealignerTargetCreator -R $ref_genome.fasta -I $OUTPUT-RG_Sorted_dedup.bam -o $OUTPUT.intervals

			# 11. Locally realign reads around indels.
			module load java
			java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=~/tmp -jar $GATK -T IndelRealigner -R $ref_genome.fasta -I $OUTPUT-RG_Sorted_dedup.bam -targetIntervals $OUTPUT.intervals -o $OUTPUT-RG_Sorted_dedup_realigned.bam

			# 12. Clip overlapping read pairs.
			/N/soft/rhel6/bamUtil/1.0.13/bam clipOverlap --in $OUTPUT-RG_Sorted_dedup_realigned.bam --out $OUTPUT-RG_Sorted_dedup_realigned_Clipped.bam

			# 13. Index the clipped BAM file using Samtools
			$samtools index $OUTPUT-RG_Sorted_dedup_realigned_Clipped.bam

			# 14. Make the mpileup file from the BAM file.
			$samtools mpileup -f $ref_genome.fasta $OUTPUT-RG_Sorted_dedup_realigned_Clipped.bam > $OUTPUT.mpileup
		"
}else{  
		print ", Ops, this file is not found! \n";  
	}  
}
print "\n\nTotal number pbs files produced:".$n1."\n\n\n";

