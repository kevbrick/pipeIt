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
params.in             = ''
params.dnase          = true
params.outdir         = "${launchDir}/output"

// import modules for pipeline
include { pairs_to_bam; clean_dnase_bam; make_mnase_bigwig } from "./modules/4dnucleome.modules.nf" \
  params(outdir: params.outdir, \
         dnase: params.dnase)

// OK ... let's start
workflow {
  def allPairFilesForStats = Channel.fromPath(params.in)
                                    .map { row ->
                                      def name=row.name.replaceFirst(~/\.pairs.gz$/, '')
                                      return [name,name,file(row)]}

  allPairFilesForStats.view()

  pairs_to_bam(allPairFilesForStats) | clean_dnase_bam | make_mnase_bigwig
}
