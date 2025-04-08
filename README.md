# A pipeline for calling SVs from long-read sequencing data (deprecated)

This pipeline is designed for calling somatic structural and small variants from long-read sequencing data (e.g. Pacbio-HiFi and Ont). 
It works on LSF job scheduler and can run multiple jobs in parallel.
Basically, this pipeline is composed of two steps: `LR_0_StrucVar.pl` for structral variants, and `LR_1_SmallVar.pl` for small variants.

Usage: a. perl LR_0_StrucVar.pl config4StrucVar.tsv    b. perl LR_1_SmallVar.pl config4SmallVar.tsv

## Prerequisite ##
**Option 1:** set a environment for LSF job on compute1 by adding the following to ~/.bashrc file:

    export PATH=/rdcw/fs1/dinglab/Active/Projects/yuweiz/anaconda/envs/longread_cv/bin/:$PATH
    
    export STORAGE1=/rdcw/fs1/dinglab/Active export SCRATCH1=/scratch1/fs1/dinglab/Active
    
    export LSF_DOCKER_VOLUMES="$STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"
    
  do NOT forget to run `source ~/.bashrc`
  
**Option 2:** modify the configure file (config.tsv) to specify the location of the softwares (minimap2 & Sniffles2) and the minimap2 index (if available) in the configure file.
    for example, MINIMAP2	${your_path}/minimap2

## Step 1 ## 
  prepare your sample list file (sample.lst)
  
      format: #id    sequencing_platform    path_of_fastq
  
  Note: specify the sequencing platform as the parameter -x used in minimap2 (https://lh3.github.io/minimap2/minimap2.html). 
  
## Step 2 ## 
  modify the parameters configure file (config.tsv), inclduing the paths of output "OUTDIR", sample list "SAMPLE", softwares and index and the bsub setting. 

## Step 3 ##
  run **perl LR_0_StrucVar.pl config4StrucVar.tsv**
  
  Please be sure the output dir is writtable and all softwares can be invoked. 

  Take care of the Sniffles mode!!! This pipeline defaultly runs both basic and mosaic (for low-frequency/non-germline SVs) modes.

## Step 4 ##
  Small variant calling is based on ClairS (https://github.com/HKU-BAL/ClairS). 
  This step requires the output bam file from Minimap2!
  
  a. load new environment: 
  
  `export PATH=/rdcw/fs1/dinglab/Active/Projects/yuweiz/anaconda3/envs/clair3/bin:$PATH`

  b. set up a configure file and a sample list. 

      sample list format: #id    tumor_bam_path    normal_bam_path (if none, will run in the tumor only mode)
  
  c. run **perl LR_1_SmallVar.pl config4SmallVar.tsv**

  
## Contact ##
Yuwei ZHANG ywhang0713@gmail.com
