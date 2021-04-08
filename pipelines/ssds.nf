nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================================================="
  log.info "SSDS PIPELINE (Version 2.0)                                "
  log.info "=========================================================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "The pipeline is run using a parent perl script: \\"
  log.info "/data/RDCO/code/nxfDSL2/pipeIt2 --h for help \\ "
  log.info " "
  log.info "---------------------------------------------------------------------------------------------------------- "
  log.info "nextflow run \$DSL2DIR/ssds.nf \\"
  log.info " --bam           <bam file>"
  log.info " --sra           <sra names (or comma separated list: align as one file)>"
  log.info " --fq1           <Read1 fastq file>"
  log.info " --fq2           <Read2 fastq file>"
  log.info " --obj           <Sample name regex for RDCO object store>"
  log.info " --genome        <string> \\"
  log.info " --name          <string default=bam file stem> \\"
  log.info " --outdir        <string: default = outName> \\"
  log.info " "
  log.info "HELP: nextflow run \$DSL2DIR/align.nf --help"
  log.info " "
  log.info "=========================================================================================================="
  log.info "Required Arguments:"
  log.info " "
  log.info "   --aligner		 STRING     bwa / bwaaln / bowtie / mm2 (default = bwa (uses mem) )" OR
  log.info "   --bam				 STRING     BAM file" OR
  log.info "   --sra				 STRING     SRA ACCESSION(s) ; for multiple, separate with commas" OR
  log.info "   --obj  			 STRING     Sample name regex to search for in object storage" OR
  log.info "   --fq1				 STRING     Read1 FASTQ / FASTQ.GZ"
  log.info "   --fq2				 STRING     Read2 FASTQ / FASTQ.GZ "
  log.info "   --genome      STRING     reference genome name"
  log.info " "
  log.info "Output and Tempory directory Arguments"
  log.info "   --name       STRING     Output file prefix"
  log.info "   --outdir     STRING     Output directory"
  log.info " "
  log.info "=========================================================================================================="
  log.info "Optional Pipeline Arguments:"
  log.info "=========================================================================================================="
  log.info " "
  log.info "   --R1len           STRING    trim length for read 1 (default = 36 )"
  log.info "   --R2len           STRING    trim length for read 2 (default = 40 )"
  log.info "   --original        STRING    use the original SSDS pipeline with BWA-RA"
  log.info "                               ** NEW PIPELINE IS RECOMMENDED ** "
  log.info "                               Pros: much faster, more complete capture of ssDNA"
  log.info "                                   : works for long reads OR short reads"
  log.info "                               Cons: not completely tested (but it appears good)"
  log.info "   --splitSz         NUMERIC   split size for fastqs (default = 10000000)\""
  log.info "   --genomes2screen  STRING    genomes to screen with fastqScreen"
  log.info "                               DEFAULT = mm10,hg38,rn6,sacCer3,canFam3,monDom5,ecoli,bsub214,phiX,UniVec_June15,"
  log.info ""
  log.info "=========================================================================================================="
  log.info "Genomes & indices:"
  log.info "=========================================================================================================="
  log.info " "
  log.info "   --genome_fastq    PATH      reference genome fasta file \""
  log.info "   --genome_fai      PATH      reference genome fasta.fai index file \""
  log.info "   --genome_bwaidx   PATH      path to bwa 0.7 index \""
  log.info "   --genome_mm2idx   PATH      path to minimap2 index \""
  log.info " "
  log.info "=========================================================================================================="
  exit 1

  }

// Parse args'
params.name           = 'ssds_test_output'
params.outdir         = "${launchDir}/output"
params.genome         = ''
params.genomes2screen = 'mm10,hg38,rn6,sacCer3,canFam3,monDom5,ecoli,bsub214,phiX,UniVec_June15'
params.r1Len          = 36
params.r2Len          = 40
params.splitSz        = 10000000
params.original       = false
params.gzipoutput     = false
params.sortFQ         = false

params.sra            = ''
params.obj            = ''
params.bam            = ''
params.fq1            = ''
params.fq2            = ''

params.genome_fasta   = "\$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa"
params.genome_fai     = "\$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa.fai"
params.genome_bwaidx  = "\$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa"
params.genome_mm2idx  = "\$NXF_GENOMES/${params.genome}/Minimap2Index/genome.fa"

def inputType
if (params.sra){inputType = 'sra'}
if (params.obj){inputType = 'obj'}
if (params.bam){inputType = 'bam'}
if (params.fq2){inputType = 'fqpe'}

def idx = params.original ? params.genome_bwaidx : params.genome_mm2idx

//log.info
log.info "===================================================================="
log.info "SSDS PIPELINE 2.0 : "
log.info "===================================================================="
log.info "name               : ${params.name}"
log.info "outdir             : ${params.outdir}"
log.info "ref genome         : ${params.genome}"
log.info "genome fasta       : ${params.genome_fasta}"
log.info "genome idx         : ${params.genome_fai}"
if (params.original) {
  log.info "pipeline version   : Orignal"
  log.info "bwa index          : ${params.genome_bwaidx}"
  log.info "Read 1 trim        : ${params.r1Len} bp"
  log.info "Read 2 trim        : ${params.r2Len} bp"
}else{
  log.info "pipeline version   : 2.0"
  log.info "mm2 index          : ${params.genome_mm2idx}"
}
log.info "source type        : ${inputType}"
if (params.bam){log.info "bam                : ${params.bam}"}
if (params.sra){log.info "sra                : ${params.sra}"}
if (params.obj){log.info "obj                : ${params.obj}"}
if (params.fq1){log.info "fq1                : ${params.fq1}"}
if (params.fq2){log.info "fq2                : ${params.fq2}"}
log.info "FQ split size      : ${params.splitSz}"
log.info "FQ screen genomes  : ${params.genomes2screen}"
log.info "--------------------------------------------------------------------"
log.info "Work folder        : ${workflow.workDir}"
log.info "Config files       : ${workflow.configFiles}"
log.info "===================================================================="
log.info ""

// import modules for pipeline
// import MUST be AFTER parameter definitions
include { getFQs } from "${projectDir}/modules/getFQ.modules.nf" \
  params(inputType: inputType, \
         genome: params.genome, \
         outdir: params.outdir , \
         bam: params.bam , \
         sra: params.sra, \
         obj: params.obj, \
         fq1: params.fq1 , \
         fq2: params.fq2, \
         name: params.name, \
         sortFQ: params.sortFQ, \
         gzipoutput: params.gzipoutput, \
         genomes2screen: params.genomes2screen)

//include { getFQs }                                                                                        from "${projectDir}/modules/getFQ.modules.nf"
//include { getFQs; SRAtoFQ; BAMtoFQ; FQtoFQpe; FQtoFQsr; OBJtoFQ; mergeFQ; fastqC; fastqScreen}            from "${projectDir}/modules/getFQ.modules.nf"
include { ssdsAlign; mergeBAMssds; parseITRs; gatherITROutputs; makeFRBWssds; makeSSreport; multiQCssds } from "${projectDir}/modules/SSDS.modules.nf"
include { trimFASTQpe; getPicardMetrics; makeDeeptoolsBW; samStats; }                                     from "${projectDir}/modules/align.modules.nf"
include { multiQC }                                                                                       from "${projectDir}/modules/generalRDCO.modules.nf"

// OK ... let's start
workflow {
  def msg = params.original ?
           "NOTE: Using original BWA-RA SSDS pipeline for short reads (trim to 36-40 bp PE) ..." : \
           "Using v2.0 SSDS pipeline ... to use the original pipeline, please use --original"

  println("${msg}")

  fastq = getFQs()

  // Trim FASTQs ... align ... and merge alignment BAMs
  // scatter/gather is to assure constant memory footprint
  trimFASTQpe(fastq.fq) | splitFastq ( by: params.splitSz, pe: true, file:true) |ssdsAlign |collect |mergeBAMssds

  //Parse by type
  Channel.from(["ssDNA_type1","ssDNA_type2","dsDNA","dsDNA_strict","unclassified"]).set {pType}
  parseITRs(ssdsAlign.out)
  gatherITROutputs(parseITRs.out.pBam.collect(),parseITRs.out.pBed.collect(),pType)

  getPicardMetrics(mergeBAMssds.out.bam,"PE")

  gatherITROutputs.out.bambai | (makeDeeptoolsBW & samStats & makeFRBWssds)

  makeSSreport(mergeBAMssds.out.bam, gatherITROutputs.out.bed.collect()) | multiQCssds

  multiQC(mergeBAMssds.out.mdmetrics.mix(fastq.fqc, \
                                           fastq.fqscr, \
                                           getPicardMetrics.out.report, \
                                           makeDeeptoolsBW.out.tab, \
                                           samStats.out.report, \
                                           gatherITROutputs.out.report).collect())
  }
