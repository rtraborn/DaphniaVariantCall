## A repository containing the Daphnia populations variant calling pipeline
### Lynch Lab Biodesign CME
#### Curated and updated by R. Taylor Raborn <rtraborn@asu.edu>
#### Initialized May 23, 2018

### Running this code:
#### To run this code on your system, please do the following:
1. Clone the repository on your system as follows:
   ```
   git clone https://github.com/LynchLab/DaphniaVariantCall
   cd DaphniaVariantCall
   ```
2. Install `NGSUtils` on your system. Please consult `INSTALL.txt` for instructions.
3. Update the links in the pipeline (e.g. `original_pipeline.sh`) that you would like to test using your favorite text editor.
4. Execute the pipeline using the appropriate batch script on your cluster as follows:
   ```
   # To run the original pipeline (i.e. with Novoalign)
   # Please make changes to the pipeline (original_pipeline_novoalign.sh) and batch files as needed for your system.
   # The batch scripts are set up for the PBS job scheduler on IU Carbonate. 
   qsub original_batch_novoalign.sh
   # To run the Hisat2 pipeline:
   # Please make changes to the pipeline (hisat2_pipeline.sh) and batch files as needed for your system.
   qsub hisat2_batch.sh                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
   
  
  
   
   
   ```
