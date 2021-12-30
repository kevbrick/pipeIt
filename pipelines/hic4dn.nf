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
params.saveBAM        = false
params.splitSz        = 50000000
params.aligner        = ""
params.vcf            = '/data/RDCO/annotation/mm10/SNV/mgp.v5.merged.snps_all.dbSNP142.B6vCAST.vcf'
params.blacklist      = ''
params.genomes2screen = 'mm10,hg38,univec'
params.gzipoutput     = true
params.sortFQ         = false
params.justPairs      = false
params.saveAllPairs   = true

// KB 03-14-2021: Use igenomes by default; use NXF_GENOMES if igenomes is not available
def bwaidx  = params.genomes["${params.genome}"]?params.genomes["${params.genome}"]["bwa"]:"\$NXF_GENOMES/${params.genome}/BWAIndex/version0.7.10/genome.fa"
def bt2idx  = params.genomes["${params.genome}"]?params.genomes["${params.genome}"]["bowtie2"]:"\$NXF_GENOMES/${params.genome}/Bowtie2Index/genome"
def fasta   = params.genomes["${params.genome}"]?params.genomes["${params.genome}"]["fasta"]:"\$NXF_GENOMES/${params.genome}/genome.fa"
def fai     = "${fasta}.fai"

def inputType
if (params.sra){inputType = 'sra'}
if (params.obj){inputType = 'obj'}
if (params.bam){inputType = 'bam'}
if (params.fq1){inputType = 'fqsr'}
if (params.fq2){inputType = 'fqpe'}

log.info "= GENOME DATA ===========================================================================\n"
log.info "FASTA:  ${fasta}"
log.info "FAI:    ${fai}"
if (params.aligner == 'bwa'){
  log.info "BWAIDX: ${bwaidx}"
}else{
  log.info "BOWTIE IDX: ${bt2idx}"
}
log.info "\n"
log.info "FASTQ split size = ${params.splitSz}"
log.info "=========================================================================================\n"

// import modules for pipeline
include { bwa4D; pairsam; mergepairs_dedup ; pairsQC; mergepairs_nodedup;
          get_RE_file; addfrag2pairs; run_cooler;
          run_juicebox_pre; cool2multirescool; hicnormvector_to_mcool ;
          bowtie2_end2end; trim_hic_reads; bowtie2_on_trimmed_reads;
          bowtie2_mergeR1R2; bowtie2_make_paired_bam;
          fixMDtags; markAlleleOfOrigin; splitByPhase;
          pairs_to_bam; clean_dnase_bam; make_mnase_bigwig ;
          pairtools_stats; multiqc ; filter_pairs; concatenate_phasing_reports;
          concatenate_allelicstatus_reports; balance_cool_matrix;
          call_compartments_from_cool } from "./modules/4dnucleome.modules.nf" \
  params(bwaidx: bwaidx, \
         bt2idx: bt2idx, \
         re: params.re, \
         fasta: fasta, \
         fai: fai, \
         outdir: params.outdir, \
         dnase: params.dnase, \
         phased: params.phased, \
         saveBAM: params.saveBAM, \
         saveAllPairs: params.saveAllPairs, \
         vcf: params.vcf, \
         blacklist: params.blacklist, \
         ligation_site: params.ligationsite)

include { getFQs; fastqC } from "./modules/getFQ.modules.nf" \
 params(inputType: inputType, \
       outdir: params.outdir, \
       genome: params.genome, \
       bam: params.bam, \
       sra: params.sra, \
       obj: params.obj, \
       fq1: params.fq1 , \
       fq2: params.fq2, \
       name: params.name, \
       sortFQ: params.sortFQ, \
       gzipoutput: params.gzipoutput, \
       genomes2screen: params.genomes2screen, \
       aligner: params.aligner)

// OK ... let's start
workflow {
  if (params.ss){
    // Sample sheet is required for this pipeline if you wish to merge samples and replicates
    samplesheet = Channel.fromPath( file(params.ss) )
           .splitCsv(header: true, sep: '\t')
           .map{row ->
               def seqrep = row['SeqRep'] == 'NA' ? '' : "_${row['SeqRep']}"
               def biorep = row['BioRep'] == 'NA' ? '' : "_${row['BioRep']}"
               def id = "${row['Sample']}${seqrep}${biorep}"
               def sample = "${row['Sample']}"
               def sequencing_replicate = row['SeqRep']
               def biological_replicate = row['BioRep']
               def reads1 = file("${row['R1']}")
               def reads2 = file("${row['R2']}")
               return [ id, sample, sequencing_replicate, biological_replicate, reads1, reads2 ]}
           .splitFastq(by: params.splitSz, pe:true, file:true,compress:true)
  }else{
    // Otherwise, just use usual pipeIt stuff
    fastqs = getFQs()

    samplesheet = fastqs.fq.map{row ->
        def id = row[0].name.replaceFirst('.merged.R1.fastq.+$','')
        def sample = id
        def sequencing_replicate = 'NA'
        def biological_replicate = 'NA'
        def reads1 = row[0]
        def reads2 = row[1]
        return [ id, sample, sequencing_replicate, biological_replicate, reads1, reads2 ]}
    .splitFastq(by: params.splitSz, pe:true, file:true, compress:true)
  }

  // Get RE file (if available) ... otherwise, make one
  if (params.dnase){
    re = Channel.fromPath('none.txt')
  }else{
    if (params.refile){
      re = Channel.fromPath(params.refile)
    }else{
      re = get_RE_file()
    }
  }

  // ALIGNMENT STEP
  // BOWTIE USES HiCPro-style alignment
  // BWA uses 4D nucleome
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

  // MARK DUPLICATES
  // IF PHASED, MARK ALLELE OF ORIGIN TOO
  if (params.phased){
    println("PHASING")
    fixMD = fixMDtags(aln)
    bamMDOK = markAlleleOfOrigin(fixMD.bam)
  }else{
    bamMDOK = fixMDtags(aln)
  }

  // MAKE INITIAL PAIRS FILE
  initPairs = pairsam(bamMDOK.bam)

  if (params.phased){
    samPEall    = mergepairs_dedup(initPairs.pairs.groupTuple(by: 0))
    samPEphased = splitByPhase(samPEall)
    samPEpair   = samPEall.mix(samPEphased.pairs)
  }else{
    samPEpair   = mergepairs_dedup(initPairs.pairs.groupTuple(by: 0))
  }

  // COLLAPSE SEQUENCING REPLICATES & SPLITFQs INTO A TUPLE
  allPairFilesForStats = samPEpair.transpose().map {row ->
                  def type = row[1].name.replaceFirst("^.+\\.(all|hom1|hom2|hom|het)\\.pairs.+\$",'$1')
                  return [ "${row[0]}.${type}",  row[1]]}
               .groupTuple(by: 0)

  // GET SOME PAIRS STATS
  samPEreport    = pairtools_stats(allPairFilesForStats)

  // COLLAPSE SEQUENCING REPLICATES & SPLITFQs INTO A TUPLE
  allPairFiles = samPEpair.transpose().map {row ->
                   def type = row[1].name.replaceFirst("^.+\\.(all|hom1|hom2|hom|het)\\.pairs.+\$",'$1')
                   return [ "${row[0]}.${type}",  row[1], row[2]]}
                .groupTuple(by: 0)

  //Add frags to hic files (not dnase)
  if (params.dnase){
    validpairs  = mergepairs_nodedup(allPairFiles) | filter_pairs
  }else{
    validpairs  = mergepairs_nodedup(allPairFiles)| combine(re) | addfrag2pairs | filter_pairs
  }

  if (!params.justPairs){
    cool    = run_cooler(validpairs)
    balcool = balance_cool_matrix(cool)
    coolM   = cool2multirescool(cool)

    hic     = run_juicebox_pre(validpairs)
    hicool  = hicnormvector_to_mcool(cool.join(hic))
  }

  // Reports
  if (params.phased){
    ch_rep = bamMDOK.report.mix(samPEphased.report).map{ ["${it[0]}",it[1] ] }.groupTuple(by: 0).view()
    phasingReport = concatenate_phasing_reports(ch_rep)
    allelicReport = concatenate_allelicstatus_reports(bamMDOK.rep.groupTuple(by: 0))
  }

  if (params.ss){
    if (params.phased){
      rep = samPEreport.mix(phasingReport)
    }else{
      rep = samPEreport
    }
  }else{
    if (params.phased){
      rep = fastqs.fqc.mix(fastqs.fqscr).mix(samPEreport).mix(phasingReport)
    }else{
      rep = fastqs.fqc.mix(fastqs.fqscr).mix(samPEreport)
    }
  }

  multiqc(rep.collect())

  //
  // // MAKE INITIAL PAIRS FILE
  // bamPair     = pairsam(dedup.bam)
  //
  // if (params.phased){
  //   samPEall    = mergepairs_dedup(bamPair.pairs.groupTuple(by: 0))
  //   samPEphased = splitByPhase(samPEall)
  //   samPEpair   = samPEall.mix(samPEphased.pairs)
  // }else{
  //   samPEpair   = mergepairs_dedup(bamPair.pairs.groupTuple(by: 0))
  // }
  //
  // // COLLAPSE SEQUENCING REPLICATES & SPLITFQs INTO A TUPLE
  // allPairFilesForStats = samPEpair.transpose().map {row ->
  //                  def type = row[2].name.replaceFirst("^.+\\.(all|hom1|hom2|hom|het)\\.pairs.+\$",'$1')
  //                  return [ "${row[0]}.${type}",  "${row[1]}.${type}", row[2]]}
  //               .groupTuple(by: [0,1])
  //
  // // GET SOME PAIRS STATS
  // samPEreport    = pairtools_stats(allPairFilesForStats)
  //
  // // GENERATE DNase BAM & WIG files for DNase-HiC
  // if (params.dnase){
  //   pairs_to_bam(allPairFilesForStats) | clean_dnase_bam | make_mnase_bigwig
  // }
  //
  // allPairFiles = samPEpair.transpose().map {row ->
  //                  def type = row[2].name.replaceFirst("^.+\\.(all|hom1|hom2|hom|het)\\.pairs.+\$",'$1')
  //                  return [ "${row[0]}.${type}",  "${row[1]}.${type}", row[2], row[3]]}
  //               .groupTuple(by: [0,1])
  //
  // filtered_pairs = filter_pairs(allPairFiles)
  //
  // //se.groupTuple(by: 0).view()
  // if (params.dnase){
  //   validpairs   = mergepairs_nodedup(filtered_pairs.groupTuple(by: 0))
  // }else{
  //   validpairs   = mergepairs_nodedup(filtered_pairs.groupTuple(by: 0)) | combine(re) | addfrag2pairs
  // }
  //
  // if (!params.justPairs){
  //   cool    = run_cooler(validpairs)
  //   balcool = balance_cool_matrix(cool)
  //   coolM   = cool2multirescool(cool)
  //
  //   hic     = run_juicebox_pre(validpairs)
  //
  //   hicool  = hicnormvector_to_mcool(cool.join(hic))
  // }

}
