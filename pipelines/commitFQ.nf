nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "Commit to Object Store PIPELINE (Version 2.0)                                "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "The pipeline is run using a parent perl script: \\"
  log.info "/data/RDCO/code/nxfDSL2/pipeIt2 --h for help \\ "
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run \$DSL2DIR/commitFQ.nf \\"
  log.info " --bam           <bam file>"
  log.info " --sra           <sra names (or comma separated list: align as one file)>"
  log.info " --fq1           <Read1 fastq file>"
  log.info " --fq2           <Read2 fastq file>"
  log.info " --obj           <Sample name regex for RDCO object store>"
  log.info " --pe            pass arg if paired-end (required for --sra or --obj)"
  log.info " --genome        <string> \\"
  log.info " --name          <string default=bam file stem> \\"
  log.info " "
  log.info "HELP: nextflow run \$DSL2DIR/commitFQ.nf --help"
  exit 1

}

// Params:
params.name           = 'commited'
params.genome         = ''
params.pe             = true
params.outdir         = "out"
params.sra            = ''
params.obj            = ''
params.bam            = ''
params.fq1            = ''
params.fq2            = ''
params.extra_args     = ''
params.gzipoutput     = true
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

def inputType
if (params.sra){inputType = 'sra'}
if (params.obj){inputType = 'obj'}
if (params.bam){inputType = 'bam'}
if (params.fq1){inputType = 'fqsr'}
if (params.fq2){inputType = 'fqpe'}

//log.info
log.info "===================================================================="
log.info "COMMIT FQ PIPELINE 2.0 : Map, mark duplicates, sort and index       "
log.info "===================================================================="
log.info "name               : ${params.name}"
log.info "outdir             : ${params.outdir}"
log.info "ref genome         : ${params.genome}"
log.info "source type        : ${inputType}"
if (params.bam){log.info "bam                : ${params.bam}"}
if (params.sra){log.info "sra                : ${params.sra}"}
if (params.obj){log.info "obj                : ${params.obj}"}
if (params.fq1){log.info "fq1                : ${params.fq1}"}
if (params.fq2){log.info "fq2                : ${params.fq2}"}
log.info isPE?"sequencing type    : paired-end":"sequencing type    : single-end"
log.info "--------------------------------------------------------------------"
log.info "Work folder        : ${workflow.workDir}"
log.info "Config files       : ${workflow.configFiles}"
log.info "===================================================================="
log.info ""

// import modules for pipeline
include { getFQs; commitToObj} from "${projectDir}/modules/getFQ.modules.nf" \
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

// OK ... let's start
workflow {
  fastq = getFQs()

  commitToObj(fastq.fq)

  }
