nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "======================================================================================\n"
  log.info "merge DNAse bams (Version 0.1.0)                                "
  log.info "======================================================================================\n"
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "--------------------------------------------------------------------------------------- "
  log.info "nextflow run \$DSL2DIR/mergeDNAseBAMs.nf \\"
  log.info " --ss            Sample sheet"
  log.info " --genome        <string> \\"
  log.info " --name          <string default=bam file stem> \\"
  log.info " --outdir        <string: default = outName> \\"
  log.info " "
  log.info "======================================================================================\n";
  log.info "Required Arguments:"
  log.info "          --ss				 STRING     Sample sheet "
  log.info "          --genome     STRING     reference genome name"
  log.info "          --name       STRING     Output file prefix"
  log.info "          --outdir     STRING     Output directory"
  log.info "======================================================================================\n"
  exit 1
  }

// Params:
params.ss             = ''
params.outdir         = "${launchDir}/output"
//params.bedlist        = "/data/RDCO/hotspots/mm10/B6fXCASTm_hotspots.raw.bedgraph"
//params.TSS            = "/fdb/GENCODE/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf"

// import modules for pipeline
include { make_mnase_bigwig } from "${baseDir}/modules/4dnucleome.modules.nf" \
  params(outdir: params.outdir, \
         dnase: params.dnase)

include { mergeBAMv2 as mergeBAMs; mergeBAMv2 as mergeBAMsv2; markBAMduplicates } from "${baseDir}/modules/align.modules.nf" \
 params(outdir: params.outdir)

 include { computeMatrixRefPoint; computeMatrixRefPointMultiple; bamCoverage } from "${baseDir}/modules/deeptools.modules.nf" \
  params(outdir: params.outdir)

// OK ... let's start
workflow {
  def sample_sheet = Channel.fromPath(params.ss)
                            .splitCsv(header: true, sep: '\t')
                            .map { row ->
                                  def name="${row['Sample']}.Rep${row['BioRep']}"
                                  def bam=row["bam"]
                                  return [name, bam]}

  def annot = Channel.from(["TSS",file("/data/RDCO/hotspots/mm10/B6xCASTF1/TSSpos.bed")],
                           ["HS_B6",file("/data/RDCO/hotspots/mm10/B6xCASTF1/B6.oneMotif.500bp.bed")],
                           ["HS_F1_CAST",file("/data/RDCO/hotspots/mm10/B6xCASTF1/B6xCST_CSTHS.oneMotif.500bp.bed")],
                           ["HS_F1_B6",file("/data/RDCO/hotspots/mm10/B6xCASTF1/B6xCST_B6HS.oneMotif.500bp.bed")],
                           ["HS_CAST",file("/data/RDCO/hotspots/mm10/B6xCASTF1/CST.oneMotif.500bp.bed")])


  def annotBed = Channel.from([file("/data/RDCO/hotspots/mm10/B6xCASTF1/TSSpos.bed"),
                            file("/data/RDCO/hotspots/mm10/B6xCASTF1/B6.oneMotif.500bp.bed"),
                            file("/data/RDCO/hotspots/mm10/B6xCASTF1/B6xCST_CSTHS.oneMotif.500bp.bed"),
                            file("/data/RDCO/hotspots/mm10/B6xCASTF1/B6xCST_B6HS.oneMotif.500bp.bed"),
                            file("/data/RDCO/hotspots/mm10/B6xCASTF1/CST.oneMotif.500bp.bed")])

  // def annotHS  = Channel.fromPath(params.hotspots)
  //                        .map{row -> return ["Hotspots",file(row)]}
  // def annotTSS = Channel.fromPath(params.TSS)
  //                        .map{row -> return ["TSS",file(row)]}
  // def annot    = annotHS.mix(annotTSS)

  //annot.view()

  //sample_sheet.groupTuple(by: 0).view()

  replicateBAMs = mergeBAMsv2(sample_sheet.groupTuple(by: 0))
  mdBAM         = markBAMduplicates(replicateBAMs)

  // def sampleMerge = mdBAM.bam.groupTuple(by: 0)
  // sampleMerge.view()

  def sampleMerge = mdBAM.bam
                         .map { row -> [row[0].replaceFirst('.Rep[0-9]+',''), row[1]] }
                         .groupTuple(by: 0)

  //sampleMerge.view()

  full    = mergeBAMs(sampleMerge)
  allBAMs = mdBAM.bam.mix(full).map{row -> return [row[0],row[0],file(row[1]),file(row[2]) ] }

  bw      = make_mnase_bigwig(allBAMs)

  bw.view()
  annot.view()

  //img     = computeMatrixRefPoint(bw.combine(annot).collect())
  imgs    = computeMatrixRefPointMultiple(bw,annotBed.collect())
  //imgTSS  = computeMatrixRefPoint(bw,Channel.fromPath(params.TSS).map{row -> return ["TSS",file(row[0])]})
}
