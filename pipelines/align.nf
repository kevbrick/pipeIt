nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "Alignment PIPELINE (Version 2.0)                                "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "The pipeline is run using a parent perl script: \\"
  log.info "/data/RDCO/code/nxfDSL2/pipeIt2 --h for help \\ "
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run \$DSL2DIR/align.nf \\"
  log.info " --bam           <bam file>"
  log.info " --sra           <sra names (or comma separated list: align as one file)>"
  log.info " --fq1           <Read1 fastq file>"
  log.info " --fq2           <Read2 fastq file>"
  log.info " --obj           <Sample name regex for RDCO object store>"
  log.info " --pe            pass arg if paired-end (required for --sra or --obj)"
  log.info " --genome        <string> \\"
  log.info " --name          <string default=bam file stem> \\"
  log.info " --sample_name   <string: default = outName> \\"
  log.info " --outdir        <string: default = outName> \\"
  log.info " "
  log.info "HELP: nextflow run \$DSL2DIR/align.nf --help"
  log.info " "
  log.info "=========================================================================="
  log.info "Required Arguments:"
  log.info " "
  log.info "          --aligner		 STRING     bwa / bwaaln / bowtie / mm2 (default = bwa (uses mem) )" OR
  log.info "          --bam				 STRING     BAM file" OR
  log.info "          --sra				 STRING     SRA ACCESSION(s) ; for multiple, separate with commas" OR
  log.info "          --obj  			 STRING     Sample name regex to search for in object storage" OR
  log.info "          --fq1				 STRING     Read1 FASTQ / FASTQ.GZ"
  log.info "          --fq2				 STRING     Read2 FASTQ / FASTQ.GZ "
  log.info "          --genome     STRING     reference genome name"
  log.info " "
  log.info "Output and Tempory directory Arguments"
  log.info "          --name       STRING     Output file prefix"
  log.info "          --outdir     STRING     Output directory"
  log.info " "
  log.info "====================================================================================="
  log.info "Optional Pipeline Arguments:"
  log.info "====================================================================================="
  log.info " "
  log.info "Optional MINIMAP/BWA Alignment Arguments:"
  log.info " "
  log.info "          --splitSz         NUMERIC   split size for fastqs (default = 20000000)\""
  log.info "          --extra_args      STRING     Optional arguments eg: \"-I 250,50\""
  log.info "          --genomes2screen  STRING    genomes to screen with fastqScreen"
  log.info "                                      DEFAULT = mm10,hg38,rn6,sacCer3,canFam3,"
  log.info "                                                monDom5,ecoli,bsub214,phiX,UniVec"
  log.info ""
  log.info "====================================================================================="
  log.info "Genomes & indices:"
  log.info "====================================================================================="
  log.info " "
  log.info "          --genome_fastq    PATH      reference genome fasta file \""
  log.info "          --genome_fai      PATH      reference genome fasta.fai index file \""
  log.info "          --genome_bwaidx   PATH      path to bwa 0.7 index \""
  log.info "          --genome_mm2idx   PATH      path to minimap2 index \""
  log.info "          --genome_bt2idx   PATH      path to bowtie2 index \""
  log.info " "
  log.info "====================================================================================="
  exit 1
  }

// Params:
params.name           = 'ssds_test_output'
params.outdir         = "${launchDir}/output"
params.genome         = ''
params.genomes2screen = 'mm10,hg38,rn6,sacCer3,canFam3,monDom5,ecoli,bsub214,phiX,UniVec'
params.splitSz        = 20000000
params.aligner        = "bwa"
params.pe             = true
params.outdir         = "out"
params.sra            = ''
params.obj            = ''
params.bam            = ''
params.fq1            = ''
params.fq2            = ''
params.extra_args     = ''
params.gzipoutput     = false
params.sortFQ         = false

def isPE              = params.pe

if (params.fq1 && params.fq2){
  if (!isPE){println("** WARNING ** FQ1 and FQ2 provided !! assuming PE, despite --pe false ....")}
  isPE             = true
}

if (params.fq1 && !params.fq2){
  if (isPE){println("** WARNING ** FQ2 NOT provided !! assuming SR, despite --pe true ....")}
  isPE             = false
}

params.genome_fasta   = params.genome?"\$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa":""
params.genome_fai     = params.genome?"\$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa.fai":""
params.genome_bwaidx  = params.genome?"\$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa":""
params.genome_mm2idx  = params.genome?"\$NXF_GENOMES/${params.genome}/Minimap2Index/genome.fa":""
params.genome_bt2idx  = params.genome?"\$NXF_GENOMES/${params.genome}/Bowtie2Index/genome":""

def inputType
if (params.sra){inputType = 'sra'}
if (params.obj){inputType = 'obj'}
if (params.bam){inputType = 'bam'}
if (params.fq1){inputType = 'fqsr'}
if (params.fq2){inputType = 'fqpe'}

//log.info
log.info "===================================================================="
log.info "ALIGNMENT PIPELINE 2.0 : Map, mark duplicates, sort and index       "
log.info "===================================================================="
log.info "name               : ${params.name}"
log.info "outdir             : ${params.outdir}"
log.info "aligner            : ${params.aligner} "
log.info "ref genome         : ${params.genome}"
log.info "genome fasta       : ${params.genome_fasta}"
log.info "genome idx         : ${params.genome_fai}"
if (params.aligner =~ "bwa")   {log.info "bwa index          : ${params.genome_bwaidx}"}
if (params.aligner == "mm2")   {log.info "mm2 index          : ${params.genome_mm2idx}"}
if (params.aligner == "bowtie"){log.info "bowtie2 index      : ${params.genome_bt2idx}"}
log.info "source type        : ${inputType}"
if (params.bam){log.info "bam                : ${params.bam}"}
if (params.sra){log.info "sra                : ${params.sra}"}
if (params.obj){log.info "obj                : ${params.obj}"}
if (params.fq1){log.info "fq1                : ${params.fq1}"}
if (params.fq2){log.info "fq2                : ${params.fq2}"}
log.info isPE?"sequencing type    : paired-end":"sequencing type    : single-end"
if (params.extra_args){log.info "extra args         : ${params.extra_args}"}
log.info "FQ split size      : ${params.splitSz}"
log.info "FQ screen genomes  : ${params.genomes2screen}"
log.info "--------------------------------------------------------------------"
log.info "Work folder        : ${workflow.workDir}"
log.info "Config files       : ${workflow.configFiles}"
log.info "===================================================================="
log.info ""

// import modules for pipeline
include { getFQs; fastqC } from "./modules/getFQ.modules.nf" \
  params(inputType: inputType, \
         outdir: params.outdir, \
         genome: params.genome, \
         bam: params.bam , \
         sra: params.sra, \
         obj: params.obj, \
         fq1: params.fq1 , \
         fq2: params.fq2, \
         name: params.name, \
         sortFQ: params.sortFQ, \
         gzipoutput: params.gzipoutput, \
         genomes2screen: params.genomes2screen, \
         aligner: params.aligner)

include { alignFromFQ; makeDeeptoolsBW; makeFRBW; samStats } from "./modules/align.modules.nf" \
   params(name: params.name, \
          outdir: params.outdir, \
          aligner: params.aligner, \
          isPE: isPE, \
          splitSz: params.splitSz, \
          genome: params.genome, \
          genome_fasta: params.genome_fasta, \
          genome_fai: params.genome_fai, \
          genome_bwaidx: params.genome_bwaidx, \
          genome_mm2idx: params.genome_mm2idx, \
          genome_bt2idx: params.genome_bt2idx)

include { multiQC } \
  from "./modules/generalRDCO.modules.nf"

// OK ... let's start
workflow {
  fastqs = getFQs()
  aln    = alignFromFQ(fastqs.fq)

  mqcRep = multiQC(fastqs.fqc.mix(fastqs.fqscr,
                                    aln.repMD,
                                    aln.repST,
                                    aln.repDT,
                                    aln.repAln).collect())

  // publish:
  // aln.bam      to: "${params.outdir}/bam", mode: 'copy', overwrite: true
  // aln.repMD    to: "${params.outdir}/reports", mode: 'copy', overwrite: true
  // aln.repST    to: "${params.outdir}/reports", mode: 'copy', overwrite: true
  // aln.repDT    to: "${params.outdir}/reports", mode: 'copy', overwrite: true
  // aln.repAln   to: "${params.outdir}/reports", mode: 'copy', overwrite: true
  // fastqs.fqc   to: "${params.outdir}/reports", mode: 'copy', overwrite: true
  // fastqs.fqscr to: "${params.outdir}/reports", mode: 'copy', overwrite: true
  // mqcRep       to: "${params.outdir}/reports", mode: 'copy', overwrite: true
  // aln.bwBT     to: "${params.outdir}/bigwig", mode: 'copy', overwrite: true
  // aln.bwDT     to: "${params.outdir}/bigwig", mode: 'copy', overwrite: true
  }
