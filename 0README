======================================
	How to run:
======================================
1. Save all Your raw reads in DATA_DIR in a sub dir named as SampleID, name you files in this way: 
	SampleID-001-R1.fastq
	SampleID-001-R2.fastq
	SampleID-002-R1.fastq
	SampleID-002-R2.fastq
			...
	SampleID-100-R1.fastq
	SampleID-100-R2.fastq

2. Make index: the reference genome must be first indexed using the following command:
	======================================
		export ref_genome="PA42.4.0"
		bwa index $ref_genome.fasta 
	======================================
	
3. Edit Make_pipeline_pbs-bwa.pl, change the settings:

	Your reference genome: 
		$ref_genome="PA42.4.0"
	Your adapter sequence file name: 
		$Adapters="/N/u/xw63/Carbonate/daphnia/Bioo_Adapters.fa";
	Your data directory:
		$DATA_DIR="/N/dc2/scratch/xw63/".$SampleID;
	Your sample ID: 
		$SampleID="PA2014";
	Your number of files: 
		$MaxNumberofSamples=100;
	Your output dir: 
		$OUTPUT_DIR=$DATA_DIR."/Bwa";
	Your output filename:
		$OUTPUT=$OUTPUT_DIR."/".$SampleID."-".$nstr001;
	Your email address: $emailaddress='ouqd@hotmail.com'
	
4. Make pipeline pbs files:
	======================================
		./Make_pipeline_pbs.sh
	======================================
	
	pipeline .pbs files will be produced for each pair of reads, 
	and the pbs files are saved in the current directory:
	./SampleID-001.pbs
	./SampleID-002.pbs
			...
	./SampleID-100.pbs
	
5. Submit the pipeline pbs to the system for computing
	======================================
		chmod 755 qsub_all_pbs.sh
		./qsub_all_pbs.sh
	======================================

That's all, thank you!
