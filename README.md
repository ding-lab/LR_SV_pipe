# A pipeline for calling SVs from long-read sequencing data

This pipeline is designed for calling structural variants from 3rd generation sequencing data (e.g. Pacbio-HiFi and Ont). 
It works on LSF job scheduler and can run multiple jobs in parallel.
Basically, this pipeline is composed of two steps: mapping by minimap2 and SV calling by Sniffles2.  

Usage: perl LR_SV_pipe.pl config.tsv

## Prerequisite ##
**Option 1:** set a environment for LSF job on compute1 by adding the following to ~/.bashrc file:

    export PATH=/rdcw/fs1/dinglab/Active/Projects/yuweiz/anaconda/envs/longread_cv/bin/:$PATH
    
    export STORAGE1=/rdcw/fs1/dinglab/Active export SCRATCH1=/rdcw/fs1/dinglab/
    
    export LSF_DOCKER_VOLUMES="$STORAGE1:$STORAGE1"
    
  do NOT forget to run `source ~/.bashrc`
  
**Option 2:** modify the configure file (config.tsv)
    Specify the location of the softwares (minimap2 & Sniffles2) and the minimap2 index (if available) in the configure file.
    for example, MINIMAP2	${your_path}/minimap2

## Step 1 ## 
  prepare your sample list file (sample.lst)
  
  format: #id sequencing_platform path_of_fastq
  
  Note: specify the sequencing platform as the parameter -x used in minimap2 (https://lh3.github.io/minimap2/minimap2.html). 
  
## Step 2 ## 
  modify the parameters configure file (config.tsv), inclduing the paths of output "OUTDIR", sample list "SAMPLE", softwares and index and the bsub setting. 

## Step 3 ##
  run **perl LR_SV_pipe.pl config.tsv**
  
  Please be sure the output dir is writtable and all softwares can be invoked. 

## Contact ##
Yuwei ZHANG ywhang0713@gmail.com
