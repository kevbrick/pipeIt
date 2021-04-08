# RDCO lab: pipeIt code and nextflow DSL2 pipelines :

This repo contins nextflow DSL2 coded pipelines used by the RDCO lab @ NIDDK, NIH. 
Little effort has been put into portability, so using these in another setting may require extensive customization.

pipeIt2 is the control script to allow different pipelines to be run in a standardized manner. 
The individual pipelines are highly portable and written in nextflow DSL2

## Requirements:
- nextflow	20.10.0+
- singularity	3.7.3+

## Global variables required:
$NXF_PIPEDIR   : Path to folder containing SSDSPipeline_1.8.nf

$NXF_GENOMES   : Path to folder containing reference genomes for alignment
                 ** This folder requires a very specific structure (see below) **

$SLURM_JOBID   : Specifies the temporary subfolder to use  (see Temp folder requirements below)

### NXF_GENOMES Folder structure
The Illumina igenomes folder structure is not yet fully integrated. 

#### Alternative genomes structure
Each reference genome should be contained in a separate folder (i.e. $NXF_GENOMES/mouse_mm10). The sub-structure within this folder should be as follows:

$NXF_GENOMES/\<genome\>/genome.fa                : Genome fasta file

$NXF_GENOMES/\<genome\>/genome.fa.fai            : Index of genome fasta file (samtools faidx)

$NXF_GENOMES/\<genome\>/genome.dict              : Sequence dictionary for genome fasta file (use picard CreateSequenceDictionary)

$NXF_GENOMES/\<genome\>/BWAIndex/version0.7.10/  : BWA 0.7 index files (should also contain soft links to the three files above)

** NOTE: The genome files MUST be named genome.XXX - other names will cause errors (i.e. mm10.fa / hg19_genome.fa / etc ...)

### Temp folder requirements
The pipeline requires a high-level temporary folder called /lscratch. On a SLURM-based HPC, each job is assigned a global id ($SLURM_JOBID) and this is appended to the temp folder name for each process. This can be modified in the config.nf file. Thus, there is a requirement for :

/lscratch folder for temporary files
SLURM_JOBID global variable for each HPC job.

## USAGE:

```
pipeIt2 --g <genome name>
        --f1 <fastq read1>
        --f2 <fastq read2 (if applicable)>
        --pipe <pipeline>
```

## Direct pipeline usage

For example, the simple BWA alignment pipeline: 
```
nextflow run pipelines/align.nf
             -c config/nextflow.config.nf
             -profile singularity 
             --fq1 data.R1.fastq.gz 
             --fq2 data.R2.fastq.gz
             --pe true 
             --aligner bwa
             --genome mm10
```
