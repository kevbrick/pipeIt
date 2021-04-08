#!/usr/bin/env nextflow

nextflow.preview.dsl=2

// set some default params
params.help=""

if (params.help) {
  log.info " "
  log.info "=========================================================================="
  log.info "MSH4 Pipeline: Kevin Brick : May 20 2020       "
  log.info "=========================================================================="
  log.info " "
  log.info "USAGE: "
  log.info "------------------------------------------------------------------------- "
  log.info "nextflow run MSH4_pipe.nf \\"
  log.info " -with-trace -with-timeline -with-report"
  log.info " "
  log.info "HELP: nextflow run MSH4_pipe.nf --help"
  log.info " "
  log.info "================================================================================================================="
  log.info "Required Arguments:"
  log.info " --projectdir    default = ."
  log.info " "
  log.info "Output Arguments"
  log.info " --outdir       Output dir (default = output)"
  log.info " "
  log.info "================================================================================================================"
  exit 1
  }

//output and tempory directories
params.projectdir      	      = "./"

params.accessorydir	          = "${params.projectdir}/accessoryFiles"
params.bamdir     	          = "${params.accessorydir}/bam/"
params.codedir 	              = "${params.accessorydir}/scripts/"
params.confdir 	              = "${params.accessorydir}/nxfConfig/"
params.genomedir    	        = "${params.accessorydir}/genomeFiles/"
params.prdm9dir               = "${params.accessorydir}/prdm9CHIPSeq/"
params.motifdir               = "${params.accessorydir}/PRDM9motif/"
params.k4m3dir                = "${params.accessorydir}/publishedH3K4me3"
params.k4m312dpp              = "${params.k4m3dir}/H3K4me3_12dpp_B6_Paigen_rep1_SRR1035576.q30.bam"

params.outdir      		        = "./output"
params.outdirAnnot 		        = "${params.outdir}/annotation"
params.outdirBW   	          = "${params.outdir}/bigwig"
params.outdirPeaks	          = "${params.outdir}/peaks"
params.outdirFigs  	          = "${params.outdir}/figs"
params.outdirStrengthCompFigs = "${params.outdir}/figs/strengthComparison"

params.hmWidth     	          = 1000
params.wcWidth         	      = 500

params.getData                = false

//output and tempory directories
log.info "===================================================================="
log.info "MSH4 Pipeline: Kevin Brick : May 2020 "
log.info "===================================================================="
log.info " "
log.info "- pipeline args ----------------------------------------------------"
log.info "Accessory dir    : ${params.accessorydir}"
log.info "Output dir       : ${params.outdir}"
log.info "BAM dir          : ${params.bamdir}"
log.info "heatmap width    : ${params.hmWidth} bp"
log.info " "
log.info "--------------------------------------------------------------------"

Channel.fromPath("${params.genomedir}/mm10_genome.fa")
  .set {mm10FA}

Channel.fromPath("${params.genomedir}/mm10_genome.fa.fai")
  .set {mm10IDX}

Channel.fromPath("${params.genomedir}/mm10_w1k_s100.bed")
  .set {mm10w1ks100}

// import modules for pipeline
include { getMouseMeiosisChromatinMods; getGENCODE; getB6xCASThotspots } from "./modules/getannotation.nf"

// OK ... let's start
workflow {
  hs      = getB6xCASThotspots(mm10FA,mm10IDX,mm10w1ks100)
  //gencode = getGENCODE(mm10FA,mm10IDX,mm10w1ks100)
  //mmChrom = getMouseMeiosisChromatinMods(mm10FA,mm10IDX,mm10w1ks100)

  publish:
  hs      to: "${params.outdir}/annotation", mode: 'copy', overwrite: true
  //gencode to: "${params.outdir}/annotation", mode: 'copy', overwrite: true
  //mmChrom to: "${params.outdir}/annotation", mode: 'copy', overwrite: true

  }

//
// process makeBigWigs {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 4
//   memory '16g'
//
//   time { 1.hour }
//
// 	tag {bam}
//
// 	module 'samtools/1.9'
// 	module 'picard/2.9.2'
// 	module 'bedtools/2.27.1'
// 	module 'deeptools/3.0.1'
// 	module 'ucsc/388'
//
//   publishDir params.outdirBW , mode: 'copy', overwrite: true, pattern: '*bigwig'
//
//   input:
//   each file(bam) from bamFiles
//
//   output:
//   file("*ALL.bigwig") into bwAll
// 	file("*MN.bigwig")  into (bwMN, bwMN_a)
//
//   script:
//   //def s1   = bam.name.replaceAll(/.bwaMemPE.mm10.q30.bam/,'_q30.bwaMemPE.mm10.bam')
// 	//def stem = s1.replaceAll(/.bwaMemPE.mm10.bam/,'')
//
//   """
//   s1=`basename ${bam}`
//   s2=\${s1/.bwaMemPE.mm10.q30.bam/_q30.bwaMemPE.mm10.bam}
//   stem=\${s2/.bwaMemPE.mm10.bam/}
//
// 	bai=`readlink ${bam}`
// 	ln -s \$bai".bai" .
//
// 	bamCoverage -b ${bam} -o \$stem".ALL.bigwig" --centerReads -p 12
// 	bamCoverage -b ${bam} -o \$stem".MN.bigwig"  --centerReads -p 12 --MNase --centerReads
//   """
//   }
//
// process callZcwpw1Hotspots {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 4
//   memory '16g'
//
//   time { 3.hour }
//
// 	tag {strain}
//
// 	module 'macs/2.1.2'
// 	module 'picard/2.9.2'
// 	module 'bedtools/2.27.1'
// 	module 'ucsc/388'
//
//   publishDir params.outdirPeaks , mode: 'copy', overwrite: true
//
//   input:
// 	set file(bamGFP), file(baiGFP), file(bamK4m3), file(baiK4m3), val(strain), file(baiZCW), file(bamZCW) from bamTC
//   file(mm10FA)      from mm10FA
//   file(mm10IDX)     from mm10IDX
//   file(mm10w1ks100) from mm10w1ks100
//
//   output:
// 	file('Zcwpw1*peaks.bed')              into (zcwPk, zcwPk_spot)
// 	file('Zcwpw1*peaks*bedgraph')         into (zcwBG, zcwBG_a, zcwBG_b)
// 	file('Zcwpw1*SPoT.txt')               into txtSPOTvals
//
//   script:
//   //def bamGFP  = file(getGFPbam(sample))
//   //def naiGFP  = file(getGFPbai(sample))
//   //def bamK4m3 = file(getK4m3bam(sample))
//   //def baiK4m3 = file(getK4m3bai(sample))
//   //def strain  = getStrainName(sample)
//   //def bamZCW  = file(getZcwbam(sample))
//   //def baiZCW  = file(getZcwbai(sample))
// 	"""
//   #nm="Zcwpw1_${strain}_OLD"
// 	#macs2 callpeak -p 0.001 -t ${bamZCW} -c ${bamGFP} -g mm --name \$nm
// 	#cut -f1-3 \$nm"_peaks.narrowPeak"  >\$nm"_peaks.bed"
// 	#sortBed -g ${mm10IDX}	-i \$nm"_peaks.bed" >\$nm"_peaks.gSort.bed"
// 	#intersectBed -a \$nm"_peaks.gSort.bed" -b ${bamZCW} -c -sorted -g ${mm10IDX} |sort -k1,1 -k2n,2n >\$nm"_peaks.bedgraph"
//
//   nm="Zcwpw1_${strain}"
// 	macs2 callpeak -p 0.001 --broad-cutoff 0.001 -t ${bamZCW} -c ${bamGFP} -g mm --broad --name \$nm
// 	cut -f1-3 \$nm"_peaks.broadPeak"  >\$nm"_peaks.bed"
// 	sortBed -g ${mm10IDX}	-i \$nm"_peaks.bed" >\$nm"_peaks.gSort.bed"
//
// 	slopBed -i \$nm"_peaks.bed" -g ${mm10IDX} -l -0.5 -r -0.5 -pct |
//   #slopBed -i - -g ${mm10IDX} -l 250 -r 250 |mergeBed -i - -c 3 -o count >\$nm"_peaks.250bpMerge.bed"
//
// 	intersectBed -a \$nm"_peaks.gSort.bed" -b ${bamZCW} -c -sorted -g ${mm10IDX} |sort -k1,1 -k2n,2n >\$nm"_peaks.bedgraph"
// 	#intersectBed -a \$nm"_peaks.gSort.bed" -b ${bamZCW} -c -sorted -g ${mm10IDX} |sort -k1,1 -k2n,2n >\$nm"_peaks.250bpMerge.bedgraph"
//
//   bash ${params.codedir}/getSignalPortionOfTagsFromBAM.sh ${bamZCW} \$nm"_peaks.bedgraph" ${mm10IDX} >\$nm".SPoT.txt"
//
//   """
//   }
//
// process callH3K4me3Peaks {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 4
//   memory '16g'
//
//   time { 3.hour }
//
// 	tag {strain}
//
// 	module 'macs/2.1.2'
// 	module 'picard/2.9.2'
// 	module 'bedtools/2.27.1'
// 	module 'ucsc/388'
//
//   publishDir params.outdirPeaks , mode: 'copy', overwrite: true
//
//   input:
// 	set file(bamGFP), file(baiGFP), file(bamK4m3), file(baiK4m3), val(strain), file(baiZCW), file(bamZCW) from bamTC_b
//   file(mm10FA)      from mm10FA
//   file(mm10IDX)     from mm10IDX
//   file(mm10w1ks100) from mm10w1ks100
//
//   output:
// 	file('H3K4me3*peaks.bed')      into (k4me3Pk, k4me3Pk_spot)
// 	file('H3K4me3*peaks.bedgraph') into (k4me3BG, k4me3BG_a)
// 	file('*SPoT.txt')              into txtK4m3SPOTvals
// 	val(strain)                    into strains
//
//   script:
//   //def bamGFP  = file(getGFPbam(sample))
//   //def naiGFP  = file(getGFPbai(sample))
//   //def bamK4m3 = file(getK4m3bam(sample))
//   //def baiK4m3 = file(getK4m3bai(sample))
//   //def strain  = getStrainName(sample)
//   //def bamZCW  = file(getZcwbam(sample))
//   //def baiZCW  = file(getZcwbai(sample))
//   """
//     macs2 callpeak -p 0.001 -t ${bamK4m3} -c ${bamGFP} -g mm --name H3K4me3_${strain}
//     cut -f1-3 H3K4me3_${strain}_peaks.narrowPeak >H3K4me3_${strain}_peaks.bed
//     sortBed -g ${mm10IDX}	-i H3K4me3_${strain}_peaks.bed >H3K4me3_${strain}_peaks.gSort.bed
//     intersectBed -a H3K4me3_${strain}_peaks.gSort.bed -b ${bamK4m3} -c -sorted -g ${mm10IDX} |sort -k1,1 -k2n,2n >H3K4me3_${strain}_peaks.bedgraph
//
//     bash ${params.codedir}/getSignalPortionOfTagsFromBAM.sh ${bamK4m3} H3K4me3_${strain}_peaks.bedgraph ${mm10IDX} >H3K4me3_${strain}.SPoT.txt
//   """
//   }
//
// process checkPeakOverlaps {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 4
//   memory '16g'
//
//   time { 6.hour }
//
// 	tag {strain}
//
// 	module 'picard/2.9.2'
// 	module 'bedtools/2.27.1'
// 	module 'R/3.5.2'
//
//   publishDir params.outdirPeaks , mode: 'copy', overwrite: true
// 	publishDir params.outdirFigs,   mode: 'copy', overwrite: true, pattern: '*svg'
// 	publishDir params.outdirFigs,   mode: 'copy', overwrite: true, pattern: '*png'
//
//   input:
// 	val(strain)    from strains
// 	file(zcw)      from zcwPk.collect()
// 	file(zcwBG)    from zcwBG.collect()
// 	file(h3k4Bed)  from k4me3Pk.collect()
// 	file(h3k4BG)   from k4me3BG.collect()
// 	file(hsOne)    from hotspotsOneMotif500bp.collect()
// 	file(hs)       from hotspotsBG3Kbp.collect()
// 	file (tss)     from tssBEDc
// 	file (tes)     from gencodeTESBED
// 	file (cgi)     from mouseCGIBED
// 	file (h3k4GL)  from h3k4m3GL
// 	file (spo11BG) from spo11BG
//   file(q30BAM)   from bamQ30Files.collect()
//
//   output:
// 	val(strain)                                         into strainz2use
// 	file("singleMotifhotspots_with_Zcwpw1_peak*bed")    into hsNoZcw
// 	file("singleMotifhotspots_without_Zcwpw1_peak*bed") into hsAtZcw
// 	file('zcwpw1_overlaps*tab')                         into overlapTable
// 	file('*svg')                                        into olSVG
// 	file('*png')                                        into olPNG
//
//   script:
// 	"""
// 	zcwPk="Zcwpw1_${strain}_peaks.bed"
// 	zcwBG="Zcwpw1_${strain}_peaks.bedgraph"
// 	k4me3Pk="H3K4me3_${strain}_peaks.bed"
// 	hs="${strain}_maleHS.3Kb.bedgraph"
//   cgi="CGI.mm10.bed"
//
// 	intersectBed -a ${strain}.oneMotif.500bp.bed -b \$zcwPk -u |cut -f1-3 >singleMotifhotspots_with_Zcwpw1_peak_${strain}.bed
// 	intersectBed -a ${strain}.oneMotif.500bp.bed -b \$zcwPk -v |cut -f1-3 >singleMotifhotspots_without_Zcwpw1_peak_${strain}.bed
//
// 	echo -e "cs\tfrom\tto\tstrength\tspo11\ths\ttss\ttes\tcgi\tk4m3\tzyk4m3" >zcwpw1_overlaps_${strain}.tab
//
// 	sort -k1,1 -k2n,2n ${spo11BG} >spo11.bg
//   mapBed       -a \$zcwPk   -b spo11.bg   -c 4 -o sum |cut -f4 |perl -lane 'print (\$F[0] eq "."?0:\$F[0])' >spo11.ol
//
// 	intersectBed -a \$zcwPk   -b \$hs       -c          |cut -f4 |perl -lane 'print (\$F[0]>0?TRUE:FALSE)'    >hotspot.ol
// 	intersectBed -a \$zcwPk   -b ${tss}     -c          |cut -f4 |perl -lane 'print (\$F[0]>0?TRUE:FALSE)'    >tss.ol
// 	intersectBed -a \$zcwPk   -b ${tes}     -c          |cut -f4 |perl -lane 'print (\$F[0]>0?TRUE:FALSE)'    >tes.ol
// 	intersectBed -a \$zcwPk   -b ${cgi}     -c          |cut -f4 |perl -lane 'print (\$F[0]>0?TRUE:FALSE)'    >cgi.ol
// 	intersectBed -a \$zcwPk   -b \$k4me3Pk  -c          |cut -f4 |perl -lane 'print (\$F[0]>0?TRUE:FALSE)'    >h3k4.ol
// 	intersectBed -a \$zcwPk   -b ${h3k4GL}  -c          |cut -f4 |perl -lane 'print (\$F[0]>0?TRUE:FALSE)'    >zygoK4m3.ol
//
// 	paste \$zcwBG spo11.ol hotspot.ol tss.ol tes.ol cgi.ol h3k4.ol zygoK4m3.ol >>zcwpw1_overlaps_${strain}.tab
//
//   cp zcwpw1_overlaps_${strain}.tab zcwpw1_overlaps.tab
//   cp ${params.codedir}/drawOverlapFigs.R .
//   cp ${params.codedir}/genericFunctions.R .
//
// 	R --vanilla <drawOverlapFigs.R
//
// 	mv Zcwpw1_UpSet.png Zcwpw1_UpSet_${strain}.png
// 	mv Zcwpw1_UpSet.svg Zcwpw1_UpSet_${strain}.svg
// 	mv Zcwpw1_Venn.svg  Zcwpw1_Venn_${strain}.svg
//
// 	mv Zcwpw1_Strength_HS_v_others.png Zcwpw1_Strength_HS_v_others_${strain}.png
// 	mv Zcwpw1_Strength_HS_v_others.svg Zcwpw1_Strength_HS_v_others_${strain}.svg
//
//   """
//   }
//
// process makeHeatmaps{
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 16
//   memory '16g'
//
//   time { 5.hour }
//
//   tag {strain}
//
//   module 'deeptools/3.0.1'
//   module 'bedtools/2.27.1'
//
//   publishDir params.outdirAnnot, mode: 'copy', overwrite: true, pattern: "*bed"
//   publishDir params.outdirFigs , mode: 'copy', overwrite: true, pattern: "*png"
// 	publishDir params.outdirFigs , mode: 'copy', overwrite: true, pattern: "*svg"
//
//   input:
//   file (bwMN)    from bwMN.collect()
//   file (hsAtZcw) from hsAtZcw.collect()
//   file (hsNoZcw) from hsNoZcw.collect()
//   val(strain)    from strainz2use
//   file (tss)     from tssBEDd
//   file (tes)     from gencodeTESBED_a
//   file (cgi)     from mouseCGIBED_a
//   file (h3k4m312dppBIGWIG)
//   //file (tss)     from refseqTSSb
//
//   output:
//   file "*matrix.gz"    into hmCompOutMatrix
//   file "*png"          into hmCompOutPNG
//   file "*svg"          into hmCompOutSVG
//
//   script:
//
//   """
//   computeMatrix reference-point -R singleMotifhotspots_with_Zcwpw1_peak_${strain}.bed singleMotifhotspots_without_Zcwpw1_peak_${strain}.bed \
//   -a ${params.hmWidth} -b ${params.hmWidth} \
//   -S Zcwpw1_${strain}.MN.bigwig H3K4me3_${strain}.MN.bigwig GFP_Control_${strain}.MN.bigwig \
//   -o Zcwpw1_${strain}_MN_singleMotifHS_1Kb_vOL.matrix.gz --referencePoint center -p max
//
//   computeMatrix reference-point -R ${tss} -a ${params.hmWidth} -b ${params.hmWidth} \
//   -S Zcwpw1_${strain}.MN.bigwig H3K4me3_${strain}.MN.bigwig GFP_Control_${strain}.MN.bigwig \
//   -o Zcwpw1_${strain}_MN_vTSS.matrix.gz --referencePoint center -p max \
//   -bs 50 --skipZeros --missingDataAsZero
//
//   computeMatrix reference-point -R ${tss} ${tes} ${cgi} -a ${params.hmWidth} -b ${params.hmWidth} \
//   -S Zcwpw1_${strain}.MN.bigwig H3K4me3_${strain}.MN.bigwig GFP_Control_${strain}.MN.bigwig \
//   -o Zcwpw1_${strain}_MN_vTSSTESCGI.matrix.gz --referencePoint center -p max \
//   -bs 50 --skipZeros --missingDataAsZero
//
//   computeMatrix reference-point -R singleMotifhotspots_with_Zcwpw1_peak_${strain}.bed singleMotifhotspots_without_Zcwpw1_peak_${strain}.bed ${tss} ${tes} ${cgi} \
//   -a ${params.hmWidth} -b ${params.hmWidth} \
//   -S Zcwpw1_${strain}.MN.bigwig H3K4me3_${strain}.MN.bigwig GFP_Control_${strain}.MN.bigwig \
//   -o Zcwpw1_${strain}_MN_vALL.matrix.gz --referencePoint center -p max \
//   -bs 50 --skipZeros --missingDataAsZero
//
//   computeMatrix reference-point -R singleMotifhotspots_with_Zcwpw1_peak_${strain}.bed singleMotifhotspots_without_Zcwpw1_peak_${strain}.bed ${tss} ${tes} ${cgi} \
//   -a ${params.hmWidth} -b ${params.hmWidth} \
//   -S H3K4me3_12dpp_B6_Paigen_rep1_SRR1035576.q30.bigwig H3K4me3_${strain}.MN.bigwig \
//   -o H3K4me3_${strain}_v_12dpp_MN_vALL.matrix.gz --referencePoint center -p max \
//   -bs 50 --skipZeros --missingDataAsZero
//
//   ## MAKE PNGs
//   plotHeatmap --matrixFile Zcwpw1_${strain}_MN_singleMotifHS_1Kb_vOL.matrix.gz --colorMap PuOr \
//   -x "Dist. to center (Kb)" --refPointLabel "HS" \
//   --regionsLabel "HS: with ZCWPW1 Peak" "HS: no ZCWPW1 Peak" --outFileName Zcwpw1_${strain}_MN_singleMotifHS_1Kb_vOL.png \
//   --legendLocation upper-left --heatmapWidth 6
//
//   plotHeatmap --matrixFile Zcwpw1_${strain}_MN_vTSS.matrix.gz --colorMap PuOr \
//   -x "Dist. to TSS (Kb)" --refPointLabel "TSS" --regionsLabel "GENCODE TSS" \
//   --outFileName Zcwpw1_${strain}_MN_vTSS.png --legendLocation upper-left \
//   --heatmapWidth 6
//
//   plotHeatmap --matrixFile Zcwpw1_${strain}_MN_vTSSTESCGI.matrix.gz --colorMap PuOr \
//   -x "Dist. to TSS (Kb)" --refPointLabel "Mid" --regionsLabel "TSS" "TES" "CGI"\
//   --outFileName Zcwpw1_${strain}_MN_vTSSTESCGI.png --legendLocation upper-left \
//   --heatmapWidth 6
//
//   plotHeatmap --matrixFile Zcwpw1_${strain}_MN_vALL.matrix.gz --colorMap PuOr \
//   -x "Dist. to TSS (Kb)" --refPointLabel "Mid" --regionsLabel "HS: with ZCWPW1 Peak" "HS: no ZCWPW1 Peak" "TSS" "TES" "CGI"\
//   --outFileName Zcwpw1_${strain}_MN_vALL.png --legendLocation upper-left \
//   --heatmapWidth 6
//
//   plotHeatmap --matrixFile H3K4me3_${strain}_v_12dpp_MN_vALL.matrix.gz --colorMap PuOr \
//   -x "Dist. to site (Kb)" --refPointLabel "Mid" --regionsLabel "HS: with ZCWPW1 Peak" "HS: no ZCWPW1 Peak" "TSS" "TES" "CGI"\
//   --outFileName H3K4me3_${strain}_v_12dpp_MN_vALL.png --legendLocation upper-left \
//   --heatmapWidth 6
//
//   ## MAKE SVGs
//   plotHeatmap --matrixFile Zcwpw1_${strain}_MN_singleMotifHS_1Kb_vOL.matrix.gz --colorMap PuOr \
//   -x "Dist. to center (Kb)" --refPointLabel "HS" \
//   --regionsLabel "HS: with ZCWPW1 Peak" "HS: no ZCWPW1 Peak" \
//   --outFileName Zcwpw1_${strain}_MN_singleMotifHS_1Kb_vOL.svg \
//   --legendLocation upper-left --heatmapWidth 6
//
//   plotHeatmap --matrixFile Zcwpw1_${strain}_MN_vTSS.matrix.gz --colorMap PuOr \
//   -x "Dist. to TSS (Kb)" --refPointLabel "TSS" --regionsLabel "GENCODE TSS" \
//   --outFileName Zcwpw1_${strain}_MN_vTSS.svg --legendLocation upper-left \
//   --heatmapWidth 6
//
//   plotHeatmap --matrixFile Zcwpw1_${strain}_MN_vTSSTESCGI.matrix.gz --colorMap PuOr \
//   -x "Dist. to TSS (Kb)" --refPointLabel "Mid" --regionsLabel "TSS" "TES" "CGI"\
//   --outFileName Zcwpw1_${strain}_MN_vTSSTESCGI.svg --legendLocation upper-left \
//   --heatmapWidth 6
//
//   plotHeatmap --matrixFile Zcwpw1_${strain}_MN_vALL.matrix.gz --colorMap PuOr \
//   -x "Dist. to TSS (Kb)" --refPointLabel "Mid" --regionsLabel "HS: with ZCWPW1 Peak" "HS: no ZCWPW1 Peak" "TSS" "TES" "CGI"\
//   --outFileName Zcwpw1_${strain}_MN_vALL.svg --legendLocation upper-left \
//   --heatmapWidth 6
//
//   plotHeatmap --matrixFile H3K4me3_${strain}_v_12dpp_MN_vALL.matrix.gz --colorMap PuOr \
//   -x "Dist. to site (Kb)" --refPointLabel "Mid" --regionsLabel "HS: with ZCWPW1 Peak" "HS: no ZCWPW1 Peak" "TSS" "TES" "CGI"\
//   --outFileName H3K4me3_${strain}_v_12dpp_MN_vALL.svg --legendLocation upper-left \
//   --heatmapWidth 6
//   """
//   }
//
// process plotRecombAtZcwpw1Hotspots{
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 16
//   memory '16g'
//
//   time { 3.hour }
//
// 	tag {strain}
//
//   module 'deeptools/3.0.1'
//   module 'bedtools/2.27.1'
//
//   publishDir params.outdirFigs , mode: 'copy', overwrite: true, pattern: "*png"
// 	publishDir params.outdirFigs , mode: 'copy', overwrite: true, pattern: "*svg"
//
//   input:
// 	file (zcw)    from zcwBG_b.collect()
// 	file spo11BW
// 	file (bw)     from bwMN_a.collect()
//
//   output:
//   file "*matrix.gz"    into zcwOutMatrix
//   file "*png"          into zcwVOthersOutPNG
//   file "*svg"          into zcwVOthersSVG
//
//   script:
//
//   """
//   #!/bin/bash
// 	wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664275/suppl/GSM2664275_Testis_SSDS_T1.ssDNA_type1.bam
//
// 	bamCoverage -b GSM2664275_Testis_SSDS_T1.ssDNA_type1.bam -o B6_SSDS.bigwig --centerReads -p max
//
// 	computeMatrix reference-point -R Zcwpw1_B6_peaks.bedgraph \
// 	-a ${params.hmWidth} -b ${params.hmWidth} \
// 	-S Zcwpw1_B6.MN.bigwig H3K4me3_B6.MN.bigwig GFP_Control_B6.MN.bigwig B6_SSDS.bigwig ${spo11BW} \
// 	-o Zcwpw1_B6_vOthers.matrix.gz --referencePoint center -p max
//
// 	## MAKE PNGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vOthers.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "Zcwpw1" "H3K4m3" "GFP" "DMC1 SSDS" "Spo11 oligos" \
// 	--outFileName Zcwpw1_B6_vOthers.png \
// 	--legendLocation upper-left
//
// 	## MAKE SVGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vOthers.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "Zcwpw1" "H3K4m3" "GFP" "DMC1 SSDS" "Spo11 oligos" \
// 	--outFileName Zcwpw1_B6_vOthers.svg \
// 	--legendLocation upper-left
//
//   #########################################################################################
// 	computeMatrix reference-point -R Zcwpw1_B6_peaks.bedgraph \
// 	-a ${params.hmWidth} -b ${params.hmWidth} \
// 	-S Zcwpw1_B6.MN.bigwig \
// 	-o Zcwpw1_B6_vSelf.matrix.gz --referencePoint center -p max
//
// 	## MAKE PNGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vSelf.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "Zcwpw1" \
// 	--outFileName Zcwpw1_B6_vSelf.png \
// 	--legendLocation upper-left
//
// 	## MAKE SVGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vSelf.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "Zcwpw1" \
// 	--outFileName Zcwpw1_B6_vSelf.svg \
// 	--legendLocation upper-left
//
//   #########################################################################################
// 	computeMatrix reference-point -R Zcwpw1_B6_peaks.bedgraph \
// 	-a ${params.hmWidth} -b ${params.hmWidth} \
// 	-S H3K4me3_B6.MN.bigwig \
// 	-o Zcwpw1_B6_vH3K4.matrix.gz --referencePoint center -p max
//
// 	## MAKE PNGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vH3K4.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "H3K4m3" \
// 	--outFileName Zcwpw1_B6_vH3K4.png \
// 	--legendLocation upper-left
//
// 	## MAKE SVGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vH3K4.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "H3K4m3" \
// 	--outFileName Zcwpw1_B6_vH3K4.svg \
// 	--legendLocation upper-left
//
//   #########################################################################################
// 	computeMatrix reference-point -R Zcwpw1_B6_peaks.bedgraph \
// 	-a ${params.hmWidth} -b ${params.hmWidth} \
// 	-S GFP_Control_B6.MN.bigwig \
// 	-o Zcwpw1_B6_vGFP.matrix.gz --referencePoint center -p max
//
// 	## MAKE PNGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vGFP.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "GFP" \
// 	--outFileName Zcwpw1_B6_vGFP.png \
// 	--legendLocation upper-left
//
// 	## MAKE SVGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vGFP.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "GFP" \
// 	--outFileName Zcwpw1_B6_vGFP.svg \
// 	--legendLocation upper-left
//
//   #########################################################################################
// 	computeMatrix reference-point -R Zcwpw1_B6_peaks.bedgraph \
// 	-a ${params.hmWidth} -b ${params.hmWidth} \
// 	-S B6_SSDS.bigwig  \
// 	-o Zcwpw1_B6_vDMC1.matrix.gz --referencePoint center -p max
//
// 	## MAKE PNGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vDMC1.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "DMC1 SSDS" \
// 	--outFileName Zcwpw1_B6_vDMC1.png \
// 	--legendLocation upper-left
//
// 	## MAKE SVGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vDMC1.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "DMC1 SSDS" \
// 	--outFileName Zcwpw1_B6_vDMC1.svg \
// 	--legendLocation upper-left
//
//   #########################################################################################
// 	computeMatrix reference-point -R Zcwpw1_B6_peaks.bedgraph \
// 	-a ${params.hmWidth} -b ${params.hmWidth} \
// 	-S ${spo11BW} \
// 	-o Zcwpw1_B6_vSPO11.matrix.gz --referencePoint center -p max
//
// 	## MAKE PNGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vSPO11.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "Spo11 oligos" \
// 	--outFileName Zcwpw1_B6_vSPO11.png \
// 	--legendLocation upper-left
//
// 	## MAKE SVGs
// 	plotHeatmap --matrixFile Zcwpw1_B6_vSPO11.matrix.gz --colorMap PuOr \
// 	-x "Dist. to center (Kb)" --refPointLabel "Center" \
// 	--regionsLabel "ZCWPW1 Peak" \
// 	--samplesLabel "Spo11 oligos" \
// 	--outFileName Zcwpw1_B6_vSPO11.svg \
// 	--legendLocation upper-left
//   """
//   }
//
// process checkZcwpw1Vstrength {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 4
//   memory '16g'
//
//   time { 3.hour }
//
//   tag {strain}
//
//   module 'picard/2.9.2'
//   module 'bedtools/2.27.1'
//   module 'ucsc/388'
//   module 'R/3.5.2'
//
//   publishDir params.outdirPeaks , mode: 'copy', overwrite: true, pattern: '*tab'
//   publishDir params.outdirStrengthCompFigs , mode: 'copy', overwrite: true, pattern: '*pdf'
//   publishDir params.outdirStrengthCompFigs , mode: 'copy', overwrite: true, pattern: '*png'
//
//   input:
// 	file(zcwBG)      from zcwBG_a.collect()
//
//   file(bam)        from bamFiles_a.collect()
//   file(bQ30)       from bamQ30Files_b.collect()
//   file(idx30)      from bamQ30IDX.collect()
// 	file(k4me3BG)    from k4me3BG_a.collect()
// 	file(k36m3B6)    from h3k36m3B6
// 	file(k36m3B6S)   from h3k36m3B6Spo11
// 	file(prdm9BG)    from prdm9BG
//
//   file(mm10FA)      from mm10FA
//   file(mm10IDX)     from mm10IDX
//   file(mm10w1ks100) from mm10w1ks100
//
//   output:
// 	file('allData.tab')  into strengthTable
// 	file('*png')         into strengthVZcwPNG
// 	file('*pdf')         into strengthVZcwPDF
//
//   script:
//   """
// 	wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0492-5/MediaObjects/41586_2018_492_MOESM3_ESM.zip
// 	unzip 41586_2018_492_MOESM3_ESM.zip
// 	mv 2017-12-16654C-s3.txt BrickEtAlTable_2018.tabinit
// 	grep -vP '^(\\#|cs)' BrickEtAlTable_2018.tabinit |perl -lane '\$F[1] += 1000; \$F[2] -= 1000; print join("\\t",@F) unless(join(":",@F[55..57]) =~ /TRUE/)' |grep -P '^(cs\\s|chr[0-9]+\\s)' >BrickEtAlTable_2018.bedgraph
// 	grep -vP '^(\\#)'    BrickEtAlTable_2018.tabinit |perl -lane '\$F[1] += 1000 if (\$F[1] !~ /start/); \$F[2] -= 1000 if (\$F[1] !~ /start/); print join("\\t",@F) unless(join(":",@F[55..57]) =~ /TRUE/)' |grep -P '^(cs\\s|chr[0-9]+\\s)' >BrickEtAlTable_2018.tab
//
// 	## LEPTO H3K4me3 Peaks
//
// 	wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734408/suppl/GSM3734408%5FLE%2ER1%2EH3K4me3%2Epeaks%2Ebed%2Egz
// 	gunzip -c GSM3734408_LE.R1.H3K4me3.peaks.bed.gz |cut -f1-3 >Leptotene_H3K4me3.peaks.bed
//
// 	wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734408/suppl/GSM3734408%5FLE%2ER1%2EH3K4me3%2EmonoCorrected%2Ews25bp%2Ebigwig
//   bigWigToBedGraph GSM3734408_LE.R1.H3K4me3.monoCorrected.ws25bp.bigwig GSM3734408_LE.R1.H3K4me3.monoCorrected.ws25bp.bg
// 	mapBed -sorted -a Leptotene_H3K4me3.peaks.bed -b GSM3734408_LE.R1.H3K4me3.monoCorrected.ws25bp.bg -c 4 -o sum >h3k4m3Lep.recombMetrics.bedgraph
//
// 	## ZYGO H3K4me3 Peaks
//   wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734414/suppl/GSM3734414%5FZY%2ER1%2EH3K4me3%2Epeaks%2Ebed%2Egz
//   gunzip -c GSM3734414_ZY.R1.H3K4me3.peaks.bed.gz |cut -f1-3 >Zygotene_H3K4me3.peaks.bed
//
// 	wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734414/suppl/GSM3734414%5FZY%2ER1%2EH3K4me3%2EmonoCorrected%2Ews25bp%2Ebigwig
// 	bigWigToBedGraph GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bigwig GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bg
// 	mapBed -sorted -a Zygotene_H3K4me3.peaks.bed -b GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bg -c 4 -o sum >h3k4m3Zyg.recombMetrics.bedgraph
//
//   cut -f1-3 BrickEtAlTable_2018.bedgraph >BrickEtAlTable_2018.bed
//
//   echo 'h3k4m3Lep'     >h3k4m3Lep.OL
//   echo 'h3k4m3Zyg'     >h3k4m3Zyg.OL
//   echo 'zcwpw1'        >zcwpw1.OL
//   echo 'zcwpw1HS'      >zcwpw1HS.OL
// 	#echo 'zcwpw1HS250'   >zcwpw1HS250.OL
// 	echo 'h3k4m3CUTnRUN' >h3k4m3CUTnRUN.OL
// 	echo 'h3k4m3MM'      >h3k4m3MM.OL
// 	echo 'h3k36m3'       >h3k36m3WT.OL
// 	echo 'h3k36m3SPO'    >h3k36m3SPO11.OL
// 	echo 'prdm9'         >prdm9CHIPSeq.OL
//
//   mapBed -a BrickEtAlTable_2018.bed -b H3K4me3_B6_peaks.bedgraph           -c 4 -o sum |sort -k1,1 -k2n,2n |cut -f4 |perl -pi -e 's/^\\.\$/0/' >>h3k4m3CUTnRUN.OL
//   mapBed -a BrickEtAlTable_2018.bed -b h3k4m3Lep.recombMetrics.bedgraph    -c 4 -o sum |sort -k1,1 -k2n,2n |cut -f4 |perl -pi -e 's/^\\.\$/0/' >>h3k4m3Lep.OL
//   mapBed -a BrickEtAlTable_2018.bed -b h3k4m3Zyg.recombMetrics.bedgraph    -c 4 -o sum |sort -k1,1 -k2n,2n |cut -f4 |perl -pi -e 's/^\\.\$/0/' >>h3k4m3Zyg.OL
//   mapBed -a BrickEtAlTable_2018.bed -b Zcwpw1_B6_peaks.bedgraph            -c 4 -o sum |sort -k1,1 -k2n,2n |cut -f4 |perl -pi -e 's/^\\.\$/0/' >>zcwpw1HS.OL
// 	#mapBed -a BrickEtAlTable_2018.bed -b Zcwpw1_B6_peaks.250bpMerge.bedgraph -c 4 -o sum |sort -k1,1 -k2n,2n |cut -f4 |perl -pi -e 's/^\\.\$/0/' >>zcwpw1HS250.OL
//
// 	sort -k1,1 -k2n,2n ${k36m3B6}  >k36m3.bedgraph
// 	sort -k1,1 -k2n,2n ${k36m3B6S} >k36m3S11.bedgraph
// 	sort -k1,1 -k2n,2n ${prdm9BG}  >prdm9.bedgraph
//
//
// 	mapBed -a BrickEtAlTable_2018.bed -b k36m3.bedgraph        -c 4 -o sum |sort -k1,1 -k2n,2n |cut -f4 |perl -pi -e 's/^\\.\$/0/' >>h3k36m3WT.OL
// 	mapBed -a BrickEtAlTable_2018.bed -b k36m3S11.bedgraph     -c 4 -o sum |sort -k1,1 -k2n,2n |cut -f4 |perl -pi -e 's/^\\.\$/0/' >>h3k36m3SPO11.OL
// 	mapBed -a BrickEtAlTable_2018.bed -b prdm9.bedgraph        -c 4 -o sum |sort -k1,1 -k2n,2n |cut -f4 |perl -pi -e 's/^\\.\$/0/' >>prdm9CHIPSeq.OL
//
//   sortBed -i BrickEtAlTable_2018.bed -g ${mm10IDX} >BrickEtAlTable_2018.gSort.bed
// 	intersectBed -a BrickEtAlTable_2018.gSort.bed -b Zcwpw1_B6.bwaMemPE.mm10.q30.bam  -c -g ${mm10IDX} -sorted |sort -k1,1 -k2n,2n |cut -f4 >>zcwpw1.OL
// 	intersectBed -a BrickEtAlTable_2018.gSort.bed -b H3K4me3_B6.bwaMemPE.mm10.q30.bam -c -g ${mm10IDX} -sorted |sort -k1,1 -k2n,2n |cut -f4 >>h3k4m3MM.OL
//
//   grep -vP '^#' BrickEtAlTable_2018.tab >merge.tab
//   paste merge.tab *.OL >allData.tab
//
//   #cp ${params.codedir}/drawStrengthCCs.R .
//   cp ${params.codedir}/genericFunctions.R .
// 	cp ${params.codedir}/drawScatterPlots.R .
//
//   R --vanilla <drawScatterPlots.R
//
//   """
//   }
//
// process getSPOTvals {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 4
//   memory '16g'
//
//   time { 4.hour }
//
// 	tag {bam}
//
// 	module 'picard/2.9.2'
// 	module 'bedtools/2.27.1'
// 	module 'R/3.5.2'
//
//   publishDir params.outdirPeaks , mode: 'copy', overwrite: true
//
//   input:
// 	file (zcw)     from zcwPk_spot.collect()
// 	file (h3k4Bed) from k4me3Pk_spot.collect()
// 	file (hs)      from hotspotsBED3Kbp
// 	file (hsbxc)   from hotspotsBEDBxC3Kbp
// 	file (tss)     from tssBED_spot
// 	file (tes)     from gencodeTESBED_spot
// 	file (cgi)     from mouseCGIBED_spot
//   file (bam)     from spotBAMs
//
//   file(mm10FA)      from mm10FA
//   file(mm10IDX)     from mm10IDX
//
//   output:
// 	file('*SPoT.txt') into allSPOTtxts
//
//   script:
// 	"""
// 	## GET SPOTs
// 	b="${bam}"
// 	sNm=\${b/.bam/}
//
// 	rm -f *250* *gSort*
//
// 	for bed in *bed; do
// 	  bash ${params.codedir}/getSignalPortionOfTagsFromBAM.sh ${bam} \$bed  ${mm10IDX}  >>\$sNm".SPoT.txt"
// 	done
//   """
//   }
//
// process gatherSPOTvals {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:40'
//   echo true
//   cpus 1
//   memory '16g'
//
//   time { 1.hour }
//
//   module 'R/3.5.2'
//
//   publishDir params.outdirFigs , mode: 'copy', overwrite: true, pattern: '*svg'
// 	publishDir params.outdirFigs , mode: 'copy', overwrite: true, pattern: '*pdf'
// 	publishDir params.outdirFigs , mode: 'copy', overwrite: true, pattern: '*png'
// 	publishDir params.outdirPeaks, mode: 'copy', overwrite: true, pattern: '*tab'
//
//   input:
//   file (spot) from allSPOTtxts.collect()
//
//   output:
// 	file ('allSPoTVals.tab') into spotTable
//   file ('S*png')            into spotPNG
// 	file ('S*pdf')            into spotPDF
// 	file ('S*svg')            into spotSVG
//
//   script:
// 	"""
// 	## GET SPOTs
// 	echo -e "bam\tinterval\treadsInside\treadsOutside\tpcInside\tpcOutside" >allSPoTVals.tab
//
// 	cat *txt |sort -k2,2 |grep -v IDX >>allSPoTVals.tab
//
// 	perl -pi -e 's/(.bwaMemPE.mm10.q30.bam|.3Kb.gS.bed|.mm10.gS.bed|_peaks.gS.bed|_rep1_SRR1035576.q30.bam|.1Kb.gS.bed|graph|.gS.bed)//g' allSPoTVals.tab
//
//   cp ${params.codedir}/genericFunctions.R .
// 	cp ${params.codedir}/drawSPoTTables.R .
//
// 	R --vanilla <drawSPoTTables.R
//   """
//   }
