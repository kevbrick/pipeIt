process getB6xCASThotspots {

  cpus 2
  memory '4 GB'

  time { 1.hour * task.attempt }

  errorStrategy { 'retry' }
  maxRetries 1

  input:
  path(mm10FA)
  path(mm10IDX)
  path(mm10w1ks100)

  output:
	path('*.oneMotif.500bp.bed', emit: hotspotOneMotif500bp)
	path('*_maleHS.3Kb.bedgraph', emit: hotspotsBG3Kbp)
	path('B6_maleHS.3Kb.bed', emit: hotspotsBED3Kbp)
	path('*oneMotif.500bp.bed', emit: hotspotsOneMotif500bp)
	path('*oneMotif.3Kb.bed', emit: hotspotsOneMotif3Kb)
  path('B6_maleHS.1bp.bedgraph', emit: hotspotBG1bp)
  path('B6_maleHS.500bp.bedgraph', emit: hotspotBG500bp)
  path('B6_maleHS.1Kb.bed', emit: hotspotBED)
	path('CST_maleHS.1bp.bedgraph', emit: cstHotspotBG1bp)
  path('CST_maleHS.500bp.bedgraph', emit: cstHotspotBG500bp)
	path('CST_maleHS.3Kb.bedgraph', emit: cstHotspotBG3Kbp)
  path('CST_maleHS.1Kb.bed', emit: cstHotspotBED)
	path('B6xCST.heat.bedgraph', emit: b6xcst_HeatBG)
	path('B6xCST.bias.bedgraph', emit: b6xcst_BiasBG)
  path('B6xCST.details.tab', emit: b6xcst_details)
	path('B6xCST_maleHS.3Kb.bed', emit: hotspotsBEDBxC3Kbp)

  script:
  """
	## get B6 hotspots
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2664nnn/GSM2664275/suppl/GSM2664275_Testis_SSDS_T1.DSBhotspots.bedgraph.gz
  gunzip -c GSM2664275_Testis_SSDS_T1.DSBhotspots.bedgraph.gz |cut -f1-3,6 |grep -P \'^chr[0-9]+\' >B6_maleHS.bedgraph

 	bedtools slop -l -0.5 -r -0.5 -pct -i B6_maleHS.bedgraph -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6_maleHS.1bp.bedgraph
  cut -f1-3 B6_maleHS.1bp.bedgraph                                                                                       >B6_maleHS.1bp.bed

 	bedtools slop -l 250  -r 250       -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6_maleHS.500bp.bedgraph
	bedtools slop -l 1500 -r 1500      -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6_maleHS.3Kb.bedgraph
	bedtools slop -l 1500 -r 1500      -i B6_maleHS.1bp.bedgraph          -g ${mm10IDX} |cut -f1-3 >B6_maleHS.3Kb.bed

	bedtools slop -l 500 -r 500      -i B6_maleHS.1bp.bedgraph            -g ${mm10IDX} |cut -f1-3 >B6_maleHS.1Kb.bed

  perl -lane \'\$nm=join("_",@F[0..2]); print join("\\t",@F[0..2],\$nm,\$nm,\$F[3])' B6_maleHS.500bp.bedgraph >B6HS500forFIMO.bed
  bedtools getfasta -fi ${mm10FA} -bed B6HS500forFIMO.bed -name -fo B6_maleHS.500bp.fa

	rm -rf fimo
  fimo --max-stored-scores 1000000 --thresh 1e-3 --o fimo1 ${params.accessorydir}/PRDM9motif/PRBS_B6.MEMEv4.pwm B6_maleHS.500bp.fa

  perl ${params.codedir}/getHotspotsWithSingleMotif.pl --fimo  ./fimo1/fimo.tsv --w 250 --out B6.oneMotif.500bp.bed
	bedtools slop -l -0.5 -r -0.5 -pct -i B6.oneMotif.500bp.bed -g ${mm10IDX} >B6.oneMotif.1bp.bed
	bedtools slop -l 1500 -r 1500 -pct -i B6.oneMotif.1bp.bed   -g ${mm10IDX} >B6.oneMotif.3Kb.bed

	## get CST hotspots
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1954nnn/GSM1954846/suppl/GSM1954846%5FCAST%5Fhotspots%2Etab%2Egz
  gunzip -c GSM1954846_CAST_hotspots.tab.gz |cut -f1-3,4 |grep -P \'^chr[0-9]+\' >CST_maleHS.bedgraph

 	bedtools slop -l -0.5 -r -0.5 -pct -i CST_maleHS.bedgraph -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >CST_maleHS.1bp.bedgraph
  cut -f1-3 CST_maleHS.1bp.bedgraph                                                                                       >CST_maleHS.1bp.bed

 	bedtools slop -l 250  -r 250       -i CST_maleHS.1bp.bedgraph          -g ${mm10IDX}            >CST_maleHS.500bp.bedgraph
	bedtools slop -l 1500 -r 1500      -i CST_maleHS.1bp.bedgraph          -g ${mm10IDX}            >CST_maleHS.3Kb.bedgraph

	bedtools slop -l 500 -r 500      -i CST_maleHS.1bp.bedgraph            -g ${mm10IDX} |cut -f1-3 >CST_maleHS.1Kb.bed

  perl -lane \'\$nm=join("_",@F[0..2]); print join("\\t",@F[0..2],\$nm,\$nm,\$F[3])' CST_maleHS.500bp.bedgraph >CSTHS500forFIMO.bed
  bedtools getfasta -fi ${mm10FA} -bed CSTHS500forFIMO.bed -name -fo CST_maleHS.500bp.fa

	rm -rf fimo
  fimo --max-stored-scores 1000000 --thresh 1e-3 --o fimo2 ${params.accessorydir}/PRDM9motif/PRBS_CST.MEMEv4.pwm CST_maleHS.500bp.fa

  perl ${params.codedir}/getHotspotsWithSingleMotif.pl --fimo  ./fimo2/fimo.tsv --w 250 --out CST.oneMotif.500bp.bed
	bedtools slop -l -0.5 -r -0.5 -pct -i CST.oneMotif.500bp.bed -g ${mm10IDX} >CST.oneMotif.1bp.bed
	bedtools slop -l 1500 -r 1500 -pct -i CST.oneMotif.1bp.bed   -g ${mm10IDX} >CST.oneMotif.3Kb.bed

	## get CSTXb6 hotspots
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2049nnn/GSM2049312/suppl/GSM2049312%5Fdmc1hotspots%5FB6CASTF1%2EPRDM9bc%2Etxt%2Egz

	gunzip -c GSM2049312_dmc1hotspots_B6CASTF1.PRDM9bc.txt.gz |grep -v heat |perl -lane 'print "chr".join("\\t",\$F[0],\$F[1]-500,\$F[1]+500,\$F[2])' >B6xCST.heat.bedgraph
	gunzip -c GSM2049312_dmc1hotspots_B6CASTF1.PRDM9bc.txt.gz |grep -v heat |perl -lane 'print "chr".join("\\t",\$F[0],\$F[1]-500,\$F[1]+500,(\$F[3] == NA?"0.5":\$F[3]))' >B6xCST.bias.bedgraph
  gunzip -c GSM2049312_dmc1hotspots_B6CASTF1.PRDM9bc.txt.gz |grep -v heat |perl -lane 'print "chr".join("\\t",\$F[0],\$F[1]-500,\$F[1]+500,\$F[2],\$F[3])' >B6xCST.bg

	bedtools slop -l -0.5 -r -0.5 -pct -i B6xCST.bg -g ${mm10IDX} |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6xCST_maleHS.1bp.bedgraph

 	bedtools slop -l 250  -r 250       -i B6xCST_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6xCST_maleHS.500bp.bedgraph
	bedtools slop -l 1500 -r 1500      -i B6xCST_maleHS.1bp.bedgraph          -g ${mm10IDX}            >B6xCST_maleHS.3Kb.bedgraph
	bedtools slop -l 1500 -r 1500      -i B6xCST_maleHS.1bp.bedgraph          -g ${mm10IDX} |cut -f1-3 >B6xCST_maleHS.3Kb.bed

  intersectBed -a B6xCST.bg -b CST_maleHS.1Kb.bed -c >C1.bg
  intersectBed -a C1.bg     -b B6_maleHS.1Kb.bed -c   >C2.bg

  cat C2.bg |perl -lane '\$type = ""; \$type = "CST" if ((\$F[5] > 0 && \$F[6] == 0) || (\$F[5] == 0 && \$F[6] == 0 && (\$F[4] > 0.75))); \$type = "B6" if ((\$F[5] == 0 && \$F[6] > 0) || (\$F[5] == 0 && \$F[6] == 0 && (\$F[4] < 0.25))); \$type = "Ambiguous" if (not \$type); print join("\\t",@F,\$type)' |sort -k1,1 -k2n,2n -k3n,3n >B6CST_all.tab

  grep CST B6CST_all.tab |cut -f1-5 >B6xCST.CSTHS.tab
  grep B6 B6CST_all.tab  |cut -f1-5 >B6xCST.B6HS.tab

	bedtools slop -l -0.5 -r -0.5 -pct -i B6xCST.CSTHS.tab -g ${mm10IDX} |cut -f1-3 |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6xCST.CSTHS.1bp.bedgraph
	bedtools slop -l -0.5 -r -0.5 -pct -i B6xCST.B6HS.tab   -g ${mm10IDX} |cut -f1-3 |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6xCST.B6HS.1bp.bedgraph

	bedtools slop -l 250  -r 250       -i B6xCST.CSTHS.1bp.bedgraph          -g ${mm10IDX}            >B6xCST.CSTHS.500bp.bedgraph
	bedtools slop -l 250  -r 250       -i B6xCST.B6HS.1bp.bedgraph            -g ${mm10IDX}            >B6xCST.B6HS.500bp.bedgraph

  paste B6xCST.bias.bedgraph B6xCST.heat.bedgraph |cut -f1-4,8 >B6xCST.heatbias.tmp
  intersectBed -a B6xCST.bias.bedgraph -b B6xCST.B6HS.1bp.bedgraph  -c |cut -f5 >likely_b6_defined.tab
  intersectBed -a B6xCST.bias.bedgraph -b B6xCST.CSTHS.1bp.bedgraph -c |cut -f5 >likely_cast_defined.tab

  echo -e "cs\\tfrom\\tto\\tbias\\theat\\tB6\\tCAST"                       >B6xCST.details.tab
  paste B6xCST.heatbias.tmp likely_b6_defined.tab likely_cast_defined.tab >>B6xCST.details.tab

  perl -lane \'\$nm=join("_",@F[0..2]); print join("\\t",@F[0..2],\$nm,\$nm,\$F[3])' B6xCST.CSTHS.500bp.bedgraph >BxC_CSTHS500forFIMO.bed
	perl -lane \'\$nm=join("_",@F[0..2]); print join("\\t",@F[0..2],\$nm,\$nm,\$F[3])' B6xCST.B6HS.500bp.bedgraph   >BxC_B6HS500forFIMO.bed

	bedtools getfasta -fi ${mm10FA} -bed BxC_CSTHS500forFIMO.bed -name -fo BxC_CST_maleHS.500bp.fa
	bedtools getfasta -fi ${mm10FA} -bed BxC_B6HS500forFIMO.bed   -name -fo BxC_B6_maleHS.500bp.fa

  rm -rf fimo
	fimo --max-stored-scores 1000000 --thresh 1e-3 --o fimo3 ${params.accessorydir}/PRDM9motif/PRBS_CST.MEMEv4.pwm BxC_CST_maleHS.500bp.fa
	perl ${params.codedir}/getHotspotsWithSingleMotif.pl --fimo  ./fimo3/fimo.tsv --w 250 --out B6xCST_CSTHS.oneMotif.500bp.bed
	bedtools slop -l -0.5 -r -0.5 -pct -i B6xCST_CSTHS.oneMotif.500bp.bed -g ${mm10IDX} >B6xCST_CSTHS.oneMotif.1bp.bed
	bedtools slop -l 1500 -r 1500 -pct -i B6xCST_CSTHS.oneMotif.1bp.bed   -g ${mm10IDX} >B6xCST_CSTHS.oneMotif.3Kb.bed

  rm -rf fimo
	fimo --max-stored-scores 1000000 --thresh 1e-3 --o fimo4 ${params.accessorydir}/PRDM9motif/PRBS_B6.MEMEv4.pwm BxC_B6_maleHS.500bp.fa
	perl ${params.codedir}/getHotspotsWithSingleMotif.pl --fimo  ./fimo4/fimo.tsv --w 250 --out B6xCST_B6HS.oneMotif.500bp.bed
	bedtools slop -l -0.5 -r -0.5 -pct -i B6xCST_B6HS.oneMotif.500bp.bed -g ${mm10IDX} >B6xCST_B6HS.oneMotif.1bp.bed
	bedtools slop -l 1500 -r 1500 -pct -i B6xCST_B6HS.oneMotif.1bp.bed   -g ${mm10IDX} >B6xCST_B6HS.oneMotif.3Kb.bed

	sort -k1,1 -k2n,2n -k3n,3n B6xCST_CSTHS.oneMotif.500bp.bed B6xCST_B6HS.oneMotif.500bp.bed >B6xCST.oneMotif.500bp.bed
	sort -k1,1 -k2n,2n -k3n,3n B6xCST_CSTHS.oneMotif.3Kb.bed B6xCST_B6HS.oneMotif.3Kb.bed     >B6xCST.oneMotif.3Kb.bed

  """
  }

process getGENCODE {
  cpus 2
  memory '4 GB'
  time { 1.hour * task.attempt }

  errorStrategy { 'retry' }
  maxRetries 1

  // module 'macs/2.1.2'
  // module 'bedtools/2.27.1'
  // module 'R/3.5.2'
  // module 'meme/5.0.1'
  // module 'ucsc/388'

  input:
  path(mm10FA)
  path(mm10IDX)
  path(mm10w1ks100)

  output:
  path('gencodeTSS.1Kb.bed', emit: tssBEDa)
  path('gencodeTES.1Kb.bed', emit: gencodeTESBED)
  path('gencodeGene.bed', emit: gencodeGeneBED)
	path('refseqTSS.bed', emit: refseqTSSb)

  script:
  """
  ##GENCODE
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz
  ##TSS
  perl ${params.codedir}/gencodeGTFtoTSS.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >gencodeTSS.1bp.noMerge.bed

  bedtools slop -l 500 -r 500   -i gencodeTSS.1bp.noMerge.bed -g ${mm10IDX} >gencodeTSS.1Kb.noMerge.bed
  mergeBed -i gencodeTSS.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,distinct |grep -vP '([\\+\\-],[\\+\\-])' |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX}            >gencodeTSS.1Kb.bed
  mergeBed -i gencodeTSS.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,distinct |grep -vP '([\\+\\-],[\\+\\-])' |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} |cut -f1-3 >gencodeTSS.1Kb3Col.bed

  cat gencodeTSS.1KbDets.bed |perl -lane \'print join("\\t",@F[0..2],0,@F[4..5])\'                           >gencodeTSS.1Kb.forMerge.bed

  ##TES
  perl ${params.codedir}/gencodeGTFtoTES.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} |grep -P "\\s+" >gencodeTES.1bp.noMerge.bed

  bedtools slop -l 500 -r 500   -i gencodeTES.1bp.noMerge.bed -g ${mm10IDX} >gencodeTES.1Kb.noMerge.bed
  mergeBed -i gencodeTES.1Kb.noMerge.bed -c 4,5,6 -o distinct,distinct,first |\
                                          ${params.codedir}/sortBEDByFAI.pl - \
                                          ${mm10IDX} >gencodeTES.1Kb.bed

  ##GENES
  perl ${params.codedir}/gencodeGTFtoCDS.pl gencode.vM20.annotation.gtf.gz |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >gencodeGene.noMerge.bed
  mergeBed -s -i gencodeGene.noMerge.bed -c 4,5,6 \
           -o distinct,distinct,first |${params.codedir}/sortBEDByFAI.pl - \
          ${mm10IDX} |grep -P "\\s+" >gencodeGene.bed

	##REFSEQ - move to annotation section later
	wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ncbiRefSeqCurated.txt.gz
	perl ${params.codedir}/parseRefSeq.pl ncbiRefSeqCurated.txt.gz TSS  |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >refseqTSS.bed
  """
  }

process getMouseMeiosisChromatinMods {
  cpus 2
  memory '4 GB'
  time { 1.hour * task.attempt }

  errorStrategy { 'retry' }
  maxRetries 1

  // module 'macs/2.1.2'
  // module 'bedtools/2.27.1'
  // module 'R/3.5.2'
  // module 'meme/5.0.1'
  // module 'ucsc/388'

  input:
  path(mm10FA)
  path(mm10IDX)
  path(mm10w1ks100)

  output:
	path('Zygotene_H3K4me3.peaks.bed', emit: h3k4m3GL)
	path('Zygotene_H3K4me3.peaks.bedgraph', emit: h3k4m3GLBG)
	path('H3K36m3_B6.bedgraph', emit: h3k36m3B6)
	path('H3K36m3_B6Spo11KO.bedgraph', emit: h3k36m3B6Spo11)

  script:
  """
	## ZYGO H3K4me3 Peaks
	wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734414/suppl/GSM3734414%5FZY%2ER1%2EH3K4me3%2Epeaks%2Ebed%2Egz
	gunzip -c GSM3734414_ZY.R1.H3K4me3.peaks.bed.gz |cut -f1-3 >Zygotene_H3K4me3.peaks.bed

	wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3734nnn/GSM3734414/suppl/GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bigwig
	bigWigToBedGraph GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bigwig GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bedgraph

	mapBed -a Zygotene_H3K4me3.peaks.bed -b GSM3734414_ZY.R1.H3K4me3.monoCorrected.ws25bp.bedgraph -c 4 -o sum |sortBEDByFAI.pl - ${mm10IDX} >Zygotene_H3K4me3.peaks.bedgraph

	## H3K36m3 Peaks and PRDM9 ChIP-Seq data
	wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93955/suppl/GSE93955_CHIP_H3K36me3_B6_coverage.bw
	wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93955/suppl/GSE93955_CHIP_H3K36me3_B6_spo11_coverage.bw
	wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93955/suppl/GSE93955_CHIP_PRDM9_B6_peaks.bed.gz
	bigWigToBedGraph GSE93955_CHIP_H3K36me3_B6_coverage.bw       GSE93955_CHIP_H3K36me3_B6_coverage.mm9.bedgraph
	bigWigToBedGraph GSE93955_CHIP_H3K36me3_B6_spo11_coverage.bw GSE93955_CHIP_H3K36me3_B6_spo11_coverage.mm9.bedgraph

	## GET LIFTOVER CHAIN FILE
	wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz   -O mm9ToMm10.over.chain.gz

	liftOver GSE93955_CHIP_H3K36me3_B6_coverage.mm9.bedgraph        mm9ToMm10.over.chain.gz  H3K36m3_B6.bgtmp na
	liftOver GSE93955_CHIP_H3K36me3_B6_spo11_coverage.mm9.bedgraph  mm9ToMm10.over.chain.gz  H3K36m3_B6Spo11KO.bgtmp na

  cat H3K36m3_B6.bgtmp        |sortBEDByFAI.pl - ${mm10IDX} >H3K36m3_B6.bedgraph
	cat H3K36m3_B6Spo11KO.bgtmp |sortBEDByFAI.pl - ${mm10IDX} >H3K36m3_B6Spo11KO.bedgraph

  """
  }

//
// process getCpGIslands {
//
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:200 '
//   echo true
//   cpus 1
//   memory "8G"
//
//   module 'bedtools/2.27.1'
// 	module 'ucsc/388'
//
//   time { 0.5.hour }
//   errorStrategy { 'retry' }
//   maxRetries 1
//
//   publishDir params.outdirAnnot, mode: 'copy', overwrite: true
//
//   input:
//
//   output:
//   path('CGI.mm10.bedgraph', emit: mouseCGIBG
// 	path('CGI.mm10.bed', emit: (mouseCGIBED, mouseCGIBED_a, mouseCGIBED_spot)
//
//   script:
//   """
//   wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz   -O mm9ToMm10.over.chain.gz
//
//   wget http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-mm9.txt
//
//   cut -f1-3,6 model-based-cpg-islands-mm9.txt  |grep -v start |sort -k1,1 -k2n,2n >CGI.mm9.bg
//
//   liftOver CGI.mm9.bg  mm9ToMm10.over.chain.gz  CGI.mm10.bedgraph na
//
//   cut -f1-3 CGI.mm10.bedgraph >CGI.mm10.bed
//
//   """
//   }
//
// process getSpo11OligoData {
//
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:200 '
//   echo true
//   cpus 1
//   memory "8G"
//
//   module 'bedtools/2.27.1'
// 	module 'ucsc/388'
//
//   time { 2.hour }
//   errorStrategy { 'retry' }
//   maxRetries 2
//
//   publishDir params.outdirAnnot, mode: 'copy', overwrite: true
//
//   input:
//   file(mm10FA)
//   file(mm10IDX)
//   file(mm10w1ks100)
//
//   output:
// 	path('B6_spo11Oligo.bedgraph', emit: spo11BG
// 	path('B6_spo11Oligo.bigwig', emit: spo11BW
//
//  	script:
//  	"""
//   wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2247nnn/GSM2247727/suppl/GSM2247727_B6_unique_as_unique_RPM.txt.gz
//  	gunzip -c GSM2247727_B6_unique_as_unique_RPM.txt.gz |perl -lane \'print join("\\t",\$F[0],\$F[1]-1,\$F[1],\$F[2])\' |grep -P \'^chr[0-9]+\' |${params.codedir}/sortBEDByFAI.pl - ${mm10IDX} >B6_spo11Oligo.bedgraph
//
//   sort -k1,1 -k2n,2n B6_spo11Oligo.bedgraph >spo11.s.bedgraph
//
// 	#bedtools makewindows -g ${mm10IDX} -w 1000 -s 100 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 1000)\' >mm10_w1k_s100.bed
//
//   mapBed -a ${mm10w1ks100} -b spo11.s.bedgraph -c 4 -o sum |perl -lane 'print join("\\t",\$F[0],\$F[1]+450,\$F[2]-451,(\$F[3] eq "."?"0":\$F[3]))' >spo11.w1ks100.bedgraph
// 	bedGraphToBigWig spo11.w1ks100.bedgraph ${mm10IDX} B6_spo11Oligo.bigwig
//  	"""
//   }
