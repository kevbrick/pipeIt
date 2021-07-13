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
params.dnase          = false
params.re             = 'MboI'
params.ligationsite   = 'GATCGATC'
params.refile         = ''
params.blacklist      = ''
params.saveAllPairs   = false
params.ss             = ''

def fasta   = params.genomes["${params.genome}"]?params.genomes["${params.genome}"]["fasta"]:"\$NXF_GENOMES/${params.genome}/genome.fa"
def fai     = "${fasta}.fai"

log.info "= GENOME DATA ===========================================================================\n"
log.info "FASTA:  ${fasta}"
log.info "FAI:    ${fai}"
log.info "=========================================================================================\n"

// import modules for pipeline
include { mergepairs_nodedup as mergepairs_tech_replicates;
          mergepairs_nodedup as mergepairs_bio_replicates;
          mergepairs_dedup;
          get_RE_file; addfrag2pairs; run_cooler;
          run_juicebox_pre; cool2multirescool;
          hicnormvector_to_mcool ;
          pairtools_stats; multiqc ;
          filter_pairs; concatenate_phasing_reports;
          balance_cool_matrix; call_compartments_from_cool} from "./modules/4dnucleome.modules.nf" \
  params(re: params.re, \
         fasta: fasta, \
         fai: fai, \
         outdir: params.outdir, \
         dnase: params.dnase, \
         saveAllPairs: params.saveAllPairs, \
         phased: params.phased, \
         blacklist: params.blacklist, \
         ligation_site: params.ligationsite)

// OK ... let's start
workflow {
  Channel.fromPath( file(params.ss) )
         .splitCsv(header: true, sep: '\t')

  // Sample sheet is required for this pipeline if you wish to merge samples and replicates
  samplesheet = Channel.fromPath( file(params.ss) )
         .splitCsv(header: true, sep: '\t')
         .map{row ->
             def sample               = "${row['Sample']}"
             def seqrep               = row['SeqRep'] == 'NA' ? '' : "${row['SeqRep']}"
             def techrep              = row['TechRep'] == 'NA' ? '' : "${row['TechRep']}"
             def biorep               = row['BioRep'] == 'NA' ? '' : "${row['BioRep']}"
             def pairs                = file("${row['Pairs']}")
             return [ sample, techrep, biorep, pairs ]}
//             return [ id, sample, biological_replicate, pairs ]}

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

  //println("${params.ss}")
  //samplesheet.view()

  // MAKE INITIAL PAIRS FILE
  dedup = mergepairs_dedup(samplesheet.groupTuple(by: [0,1,2]).map {row ->
                                        def name  = "${row[0]}_BioRep${row[2]}_TechRep${row[1]}"
                                        def file    = row[3]
                                        return [name, file]})

  // COLLAPSE SEQUENCING REPLICATES & SPLITFQs INTO A TUPLE
  allPairFilesForStats = dedup.map {row -> return [row[0], row[1]]}

  //allPairFiles   = dedup.map {row -> return [row[0], row[1], row[2]]}

  // GET SOME PAIRS STATS
  samPEreport    = pairtools_stats(allPairFilesForStats)

  // Pare off the last piece of name (seq Rep)
  // merge technical reps
  biological_reps = mergepairs_tech_replicates(dedup.map {row ->
                                                        def name    = row[0]
                                                        def newname = name.replaceAll("_TechRep[0-9]","")
                                                        return [newname, row[1], row[2]]}.groupTuple(by: 0))

  // Pare off the last piece of name (tech Rep)
  // merge biological reps
  full_dataset    = mergepairs_bio_replicates(biological_reps.map {row ->
                                                        def n1      = row[0]
                                                        def newname = n1.replaceAll("_BioRep[0-9]\$","")
                                                        return [newname, row[1], row[2]]}.groupTuple(by: 0))

  //Add frags to hic files (not dnase)
  if (params.dnase){
    validpairs  = full_dataset.mix(biological_reps) | filter_pairs
  }else{
    // perSamplePairs = merge_technical_replicates(filtered_pairs.groupTuple(by: 0)) | combine(re) | addfrag2pairs
    // validpairs     = merge_biological_replicates(filtered_pairs.groupTuple(by: 0)) | combine(re) | addfrag2pairs
    validpairs  = full_dataset.mix(biological_reps) | combine(re) | addfrag2pairs | filter_pairs
  }

  cool    = run_cooler(validpairs)
  balcool = balance_cool_matrix(cool)
  coolM   = cool2multirescool(cool)

  hic     = run_juicebox_pre(validpairs)
  hicool  = hicnormvector_to_mcool(cool.join(hic))

}
