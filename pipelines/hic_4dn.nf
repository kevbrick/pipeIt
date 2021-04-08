nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "======================================================================================\n"
  log.info "Hi-C 4D-Nucleome PIPELINE (Version 0.1.0)                                "
  log.info "======================================================================================\n"
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "The pipeline is run using a parent perl script: \\"
  log.info "\$DSL2DIR/pipeIt2 --h for help \\ "
  log.info " "
  log.info "--------------------------------------------------------------------------------------- "
  log.info "nextflow run \$DSL2DIR/hic_4dn.nf \\"
  log.info " --ss            Sample sheet"
  log.info " --reads         <fastq files>"
  log.info " --genome        <string> \\"
  log.info " --re            <string> \\"
  log.info " --refile        <string> \\"
  log.info " --dnase         <string> \\"
  log.info " --name          <string default=bam file stem> \\"
  log.info " --sample_name   <string: default = outName> \\"
  log.info " --outdir        <string: default = outName> \\"
  log.info " "
  log.info "HELP: nextflow run \$DSL2DIR/hic_4dn.nf --help"
  log.info " "
  log.info "======================================================================================\n";
  log.info "Required Arguments:"
  log.info " "
  log.info "          --re	     	 STRING     restriction enzyme default = MBOI" OR
  log.info "          --ss				 STRING     Sample sheet "
  log.info "          --genome     STRING     reference genome name"
  log.info " "
  log.info "          --dnase	 	   STRING     Micro-C "
  log.info "======================================================================================\n"
  log.info "Output and Tempory directory Arguments"
  log.info "          --name       STRING     Output file prefix"
  log.info "          --outdir     STRING     Output directory"
  log.info " "
  log.info "= REFERENCE ==========================================================================\n"
  log.info " This pipeline is based on the 4D nucleome protocols. Details are at :\n"
  log.info ' https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline#recheck'."\n"
  log.info "======================================================================================\n"
  exit 1
  }

// Params:
params.name           = ''
params.outdir         = "${launchDir}/output"
params.genome         = ''
params.ss             = ''
params.dnase          = false
params.re             = 'MboI'
params.ligationsite   = 'GATCGATC'
params.refile         = ''
params.splitSize      = 50000000
params.aligner        = "bwa"
params.vcf            = '/data/RDCO/annotation/mm10/SNV/mgp.v5.merged.snps_all.dbSNP142.B6vCAST.vcf'

//def file = new File(params.genomes["${params.genome}"])

// KB 03-14-2021: Use igenomes by default; use NXF_GENOMES if igenomes is not available
def bwaidx  = params.genomes["${params.genome}"]?params.genomes["${params.genome}"]["bwa"]:"\$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa"
def bt2idx  = params.genomes["${params.genome}"]?params.genomes["${params.genome}"]["bowtie"]:"\$NXF_GENOMES/${params.genome}/Bowtie2Index/genome"
def fasta   = params.genomes["${params.genome}"]?params.genomes["${params.genome}"]["fasta"]:"\$NXF_GENOMES/${params.genome}/genome.fa"
def fai     = "${fasta}.fai"

log.info "= GENOME DATA ===========================================================================\n"
log.info "FASTA:  ${fasta}"
log.info "FAI:    ${fai}"
log.info "BWAIDX: ${bwaidx}"
log.info "=========================================================================================\n"

// import modules for pipeline
include { bwa4D; pairsam; mergesampairs ; pairsQC; merge_biological_replicates;
          get_RE_file; addfrag2pairs; run_cooler;
          run_juicebox_pre; cool2multirescool; hicnormvector_to_mcool ;
          bowtie2_end2end; trim_hic_reads; bowtie2_on_trimmed_reads;
          bowtie2_mergeR1R2; bowtie2_make_paired_bam;
          fixMDtags; markAlleleOfOrigin; splitByPhase} from "./modules/4dnucleome.modules.nf" \
  params(bwaidx: bwaidx, \
         bt2idx: bt2idx, \
         re: params.re, \
         fasta: fasta, \
         fai: fai, \
         outdir: params.outdir, \
         dnase: params.dnase, \
         vcf: params.vcf, \
         ligation_site: params.ligationsite)

// Sample sheet is required for this pipeline
samplesheet = Channel.fromPath( file(params.ss) )
       .splitCsv(header: true, sep: '\t')
       .map{row ->
           def id = "${row['Sample']}_${row['SeqRep']}_${row['BioRep']}"
           def sample = "${row['Sample']}"
           def sequencing_replicate = row['SeqRep']
           def biological_replicate = row['BioRep']
           def reads1 = file("${row['R1']}")
           def reads2 = file("${row['R2']}")
           return [ id, sample, sequencing_replicate, biological_replicate, reads1, reads2 ]}
       .splitFastq(by: params.splitSize, pe:true, file:true,compress:true)

// OK ... let's start
workflow {

  if (params.dnase){
    re = Channel.fromPath('none.txt')
  }else{
    if (params.refile){
      re = Channel.fromPath(params.refile)
    }else{
      re = get_RE_file()
    }
  }

  if (params.aligner == 'bowtie'){
    if (params.dnase){
      bt2_e2e = bowtie2_end2end(samplesheet,50,0)
      alignedR1R2 = bt2_e2e.bamR1.mix(bt2_e2e.bamR2).groupTuple(by: [0,1,2,3], size: 2)
                  .map {row ->
                    def read1 = (row[5][0] =~ /.R1./) ? row[5][0] : row[5][1]
                    def read2 = (row[5][0] =~ /.R2./) ? row[5][0] : row[5][1]
                    return [row[0], row[1], row[2], row[3], read1, read2]}
    }else{
      bt2_e2e = bowtie2_end2end(samplesheet,50,1)
      bt2trim   = trim_hic_reads(bt2_e2e.unmappedR1.mix(bt2_e2e.unmappedR2)) |bowtie2_on_trimmed_reads

      // TUPLE GROUPED BY READ1/READ2
      bt2All    = bt2_e2e.bamR1.mix(bt2_e2e.bamR2.mix(bt2trim)).groupTuple(by: [0,1,2,3,4], size: 2)

      // MERGE e2e and trim for R1 & R2
      alignedR1R2 = bowtie2_mergeR1R2(bt2All).groupTuple(by: [0,1,2,3])
                                         .map {row ->
                                           def read1 = (row[4][0] =~ /.R1./) ? row[4][0] : row[4][1]
                                           def read2 = (row[4][0] =~ /.R2./) ? row[4][0] : row[4][1]
                                           return [row[0], row[1], row[2], row[3], read1, read2]}
    }

    // Join Read 1 and Read 2
    aln = bowtie2_make_paired_bam(alignedR1R2)
  }else{
    aln = bwa4D(samplesheet)
  }

  if (params.phased){
    println("PHASING")
    fixMD = fixMDtags(aln)
    dedup = markAlleleOfOrigin(fixMD.bam)
  }else{
    dedup = fixMDtags(aln)
  }

  bamPair     = pairsam(dedup.bam)

  if (params.phased){
    samPEpair   = mergesampairs(bamPair.pairs.groupTuple(by: [0,1,2])) | splitByPhase
  }else{
    samPEpair   = mergesampairs(bamPair.pairs.groupTuple(by: [0,1,2]))
  }

  se = samPEpair.transpose().map {row ->
      def type = row[2].name.replaceFirst("^.+\\.(all|hom1|hom2|hom|het)\\.pairs.+\$",'$1')
      return [ "${row[0]}.${type}", row[2], row[3]]}

  se.groupTuple(by: 0).view()
  if (params.dnase){
    validpairs   = merge_biological_replicates(se.groupTuple(by: 0))
  }else{
    validpairs   = merge_biological_replicates(se.groupTuple(by: 0)) | combine(re) | addfrag2pairs
  }

  cool    = run_cooler(validpairs)
  hic     = run_juicebox_pre(validpairs)

  coolM   = cool2multirescool(cool)

  hicool  = hicnormvector_to_mcool(cool.join(hic))
}
