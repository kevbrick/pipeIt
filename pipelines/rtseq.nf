nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "Replication timing Seq (RT-Seq) PIPELINE (Version 2.0)                                "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "The pipeline is run using a parent perl script: \\"
  log.info "/data/RDCO/code/nxfDSL2/pipeIt2 --h for help \\ "
  log.info " "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run \$DSL2DIR/rtseq.nf \\"
  log.info " --bam               <bam file>"
  log.info " --sra               <sra names (or comma separated list: align as one file)>"
  log.info " --fq1               <Read1 fastq file>"
  log.info " --fq2               <Read2 fastq file>"
  log.info " --obj               <Sample name regex for RDCO object store>"
  log.info " --pe                pass arg if paired-end (required for --sra or --obj)"
  log.info " --genome            Name of genome reference (mm10|hg38|canFam3 are the ONLY options) "
  log.info " --gcCorrection      GC calibration file (choose one from accessoryFiles/gcCalibration) "
  log.info " -                   If omitted, pipeline will generate a calibration file from the file "
  log.info " -                   provided. This is a time-consuming step. "
  log.info " --name              <string default=bam file stem> "
  log.info " --yeast             adapt pipeline for tiiiiny yeast genome "
  log.info " --outdir            <string: default = outName> \\"
  log.info " "
  log.info "HELP: nextflow run \$DSL2DIR/rtseq.nf --help"
  log.info " "
  exit 1

}

// Params:
params.name             = 'ssds_test_output'
params.genome           = ''
params.genomes2screen   = 'mm10,hg38,rn6,sacCer3,canFam3,monDom5,ecoli,bsub214,phiX,UniVec'
params.gcCorrection     = ''
params.splitSz          = 20000000
params.aligner          = "bwa"
params.pe               = true
params.outdir           = "out"
params.covPC            = "66" // Each bin in final bedgraph must have at least x% of bases with coverage
params.sra              = ''
params.obj              = ''
params.bam              = ''
params.bai              = ''
params.fq1              = ''
params.fq2              = ''
params.extra_args       = ''
params.gzipoutput       = false
params.test             = false
params.sortFQ           = false
params.skipAlignment    = false

//params.accessoryDir     = "$NXF_PIPEDIR/accessoryFiles/rtSeq"
params.accessoryDir     = "/usr/local/rtseq" // location in docker container
params.rtData           = "${params.accessoryDir}"
params.genome_mask_fa   = "${params.rtData}/genomeMASK/genome.fa"
params.pseudoReadBase   = "${params.rtData}/mapability/"

params.genome_fasta   = params.genome?"$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa":""
params.genome_fai     = params.genome?"$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa.fai":""
params.genome_bwaidx  = params.genome?"$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa":""

def isPE                = params.pe

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

// Make output name
def oname
if (params.name){
  oname           = "${params.name}.RTseq"
}else{
  if (params.bam){oname = params.bam.replaceFirst(".bam",".RTseq")}
}

def outname = oname.replaceFirst(".bwaMem(SR|PE).+",".RTseq")
def outName = "${outname}.${params.genome}"

log.info "bam               : ${params.bam}"
log.info "sa               : ${params.skipAlignment}"

def skipAlignment = params.bam ? true : false

log.info "sa               : ${skipAlignment}"

//log.info
log.info "===================================================================="
log.info "RT-Seq PIPELINE 2.0 : Align and infer RT       "
log.info "===================================================================="
log.info "name               : ${params.name}"
log.info "outdir             : ${params.outdir}"
log.info "GC-correction file : ${params.gcCorrection} "
log.info "ref genome         : ${params.genome}"
log.info "genome fasta       : ${params.genome_fasta}"
log.info "genome idx         : ${params.genome_fai}"
log.info "source type        : ${inputType}"
if (params.bam){log.info "bam                : ${params.bam}"}
if (params.sra){log.info "sra                : ${params.sra}"}
if (params.obj){log.info "obj                : ${params.obj}"}
if (params.fq1){log.info "fq1                : ${params.fq1}"}
if (params.fq2){log.info "fq2                : ${params.fq2}"}

if (!skipAlignment){
  if (params.aligner =~ "bwa")   {log.info "bwa index          : ${params.genome_bwaidx}"}
  log.info isPE?"sequencing type    : paired-end":"sequencing type    : single-end"
  log.info "FQ split size      : ${params.splitSz}"
  log.info "FQ screen genomes  : ${params.genomes2screen}"
  log.info "Skip alignment     : ${params.skipAlignment}"
}else{
  log.info "** BAM provided: Skipping alignment step "
}

log.info "--------------------------------------------------------------------"
log.info "Work folder        : ${workflow.workDir}"
log.info "Config files       : ${workflow.configFiles}"
log.info "===================================================================="
log.info ""

// import modules for pipeline
include { getFQs; fastqC } from "${baseDir}/modules/getFQ.modules.nf" \
  params(inputType: inputType, \
         outdir: params.outdir, \
         genome: params.genome, \
         bam: params.bam , \
         sra: params.sra, \
         obj: params.obj, \
         fq1: params.fq1 , \
         fq2: params.fq2, \
         name: params.name, \
         gzipoutput: params.gzipoutput, \
         genomes2screen: params.genomes2screen, \
         aligner: params.aligner)

include { alignFromFQ; makeDeeptoolsBW; makeFRBW; samStats } from "${baseDir}/modules/align.modules.nf" \
   params(name: params.name, \
          outdir: params.outdir, \
          isPE: isPE, \
          splitSz: params.splitSz, \
          genome: params.genome, \
          genome_fasta: params.genome_fasta, \
          genome_fai: params.genome_fai, \
          genome_bwaidx: params.genome_bwaidx,
          aligner: "bwa")

include { multiQC } from "${baseDir}/modules/generalRDCO.modules.nf"

include { getCoverage; mergeCoverageBedgraphs; generateMpileup; getGCstats ; getGCcorrFile} from "${baseDir}/modules/rtseq.modules.nf" \
    params(name: params.name, \
           outdir: params.outdir, \
           accessoryDir: params.accessoryDir, \
           isPE: isPE, \
           rtData: params.rtData, \
           genome: params.genome, \
           genome_fasta: params.genome_fasta, \
           genome_fai: params.genome_fai, \
           test: params.test, \
           fullChromTest: params.test, \
           sortFQ: params.sortFQ, \
           pseudoReadBase: params.pseudoReadBase,
           gcCorrection: params.gcCorrection,
           genome_mask_fa: params.genome_mask_fa,
           covPC: params.covPC)

// Get chromosomes //////////////////////////////////
if( params.test ){
    log.info "======== Using TEST chromosomes : 17-19 & X ============================="
    if (params.genome == 'hg38'){
      chromNames = Channel.from('chr20')
    }else{
      if (params.genome == 'mm10'){
        chromNames = Channel.from('chr18','chr19')
      }else{
        chromNames = Channel.from('chr10','chr11')
      }
    }
  }else{
    // Get CS names
    def fai= new File("${params.genome_fai}");
    chromosomeNames = [];
    fai.splitEachLine('\t') {
      def cs = it[0]
      if (cs =~ /(rand|Rand|Un|un|conti|EBV|chrM)/){
        println("Skipping ... ")
      }else{
        chromosomeNames << cs;
        //log.info "${cs}";
      }
    }

    chromNames = Channel.from(chromosomeNames)
  }

// OK ... let's start
workflow {
  //Do alignment if necessary
  if (skipAlignment){
    bam = Channel.value(params.bam)
           .map{ b -> [file("$b"), file("${b}.bai") ] }

    bam.view()

    printf('NOT Aligning !!!')

  }else{
    fastqs = getFQs()
    aln    = alignFromFQ(fastqs.fq)

    mqcRep = multiQC(fastqs.fqc.merge(fastqs.fqscr,
                                      aln.repMD,
                                      aln.repST,
                                      aln.repDT,
                                      aln.repAln))
    bam = aln.bam
  }

  //Get correction file if necessary
  if (params.gcCorrection){
    corrFile = getGCcorrFile()
  }else{
    corrFile = generateMpileup(bam, chromNames) | collect | getGCstats
  }
  
  bam.view()
  
  cov     = getCoverage(bam, corrFile, chromNames)
  
  mergeBG = mergeCoverageBedgraphs(outName,
                                   cov.bg.collect(),
                                   cov.gcData.collectFile(name: 'allGCData.txt', newLine: true))
  }
