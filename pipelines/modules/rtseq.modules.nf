process getCoverage {
  echo true
  cpus 2

  memory '32 GB'

  time { 3.hour * task.attempt}

  errorStrategy { 'retry' }
  maxRetries 1

  tag { chrom }

  //publishDir params.outdir,  mode: 'copy', overwrite: false

  input:
  path(bam)
  path(bai)
  path(gctab)
  val(chrom)

  output:
  path("*.w*0k.bedgraph", emit: bg)
  path("*.GCdata.tab",    emit: gcData)

  script:

  """
  tbam="${chrom}.tmp.bam"
  s1bam="${chrom}.s1.bam"
  s2bam="${chrom}.s2.bam"
  s3bam="${chrom}.s3.bam"
  s4bam="${chrom}.s4.bam"
  cbam="${chrom}.bam"

  samtools view -hb -q 30 ${bam} ${chrom} >\$tbam
  samtools index \$tbam

  java -jar -Xmx8g \$PICARDJAR SortSam I=\$tbam O=\$s1bam VALIDATION_STRINGENCY=LENIENT SO=queryname TMP_DIR="\$TMPDIR"

  nPE=`samtools view -h \$tbam |head -n 100000 |samtools view -c -f 1 -S /dev/stdin`

  if [ "\$nPE" -eq "0" ]; then
    ln -s \$s1bam \$s3bam
  else
    java -jar -Xmx8g \$PICARDJAR FixMateInformation I=\$s1bam O=\$s2bam VALIDATION_STRINGENCY=LENIENT SO=queryname AS=true TMP_DIR="\$TMPDIR"
    samtools view -hb -f 2 \$s2bam >\$s3bam
  fi

  java -jar -Xmx8g \$PICARDJAR MarkDuplicates I=\$s3bam O=\$s4bam VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true M=metrics.tab TMP_DIR="\$TMPDIR"
  java -jar -Xmx8g \$PICARDJAR SortSam I=\$s4bam O=\$cbam VALIDATION_STRINGENCY=LENIENT SO=coordinate TMP_DIR="\$TMPDIR"
  samtools index \$cbam

  #rm -f *.s?.bam

  readLen=`perl ${params.accessoryDir}/scripts/getClosestReadLength.pl \$tbam 50,150`
  win101="${params.rtData}/genomewin/${chrom}.win101.bed"
  ps101="${params.rtData}/genomewin/${chrom}.psr"\$readLen".tab"
  gc101="${params.rtData}/genomewin/${chrom}.win101.GC.tab"

  grep -w ${chrom} ${params.genome_fai} >idx.fai

  intersectBed -a \$win101 -b \$cbam -c -sorted -g ${params.genome_fai} |perl -pi -e 's/\\./0/g' >${chrom}.cover.tab

  paste ${chrom}.cover.tab \$ps101 \$gc101 >${chrom}.tab

  ## DO GC NORMALIZATION
  expectedSim=`echo \$readLen + 100 |bc`
  mapabilityLowerLim=`echo "(\$readLen + 100)*0.66" |bc`
  mapabilityUpperLim=`echo "(\$readLen + 100)*1.2" |bc`

  perl ${params.accessoryDir}/scripts/normalizeBy2NGC.pl ${gctab} ${chrom}.tab \$expectedSim

  shuf GCData.tab |head -n 1000000 |sort -k1rn,1rn >${chrom}.GCdata.tab

  perl -lane 'print join("\\t",@F[0..2],\$F[4]) if (\$F[4] > '\$mapabilityLowerLim' && \$F[4] <  '\$mapabilityUpperLim' )' ${chrom}.tab >${chrom}.mapability.tab

  bedtools makewindows -g idx.fai -w 500000 -s 50000  | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500000)\' >win500k50k.init.bed
  bedtools makewindows -g idx.fai -w 500000 -s 100000 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500000)\' >win500k100k.init.bed

  mapBed -a win500k50k.init.bed  -b ${chrom}.mapability.tab -c 4 -o count |perl -lane 'print join("\\t",@F[0..2]) if (\$F[3] > 330000)' >win500k50k.bed
  mapBed -a win500k100k.init.bed -b ${chrom}.mapability.tab -c 4 -o count |perl -lane 'print join("\\t",@F[0..2]) if (\$F[3] > 330000)' >win500k100k.bed

  mapBed -a win500k50k.bed  -b ${chrom}.GCcorrected.bedgraph -c 4,4 -o sum,count |perl -M"Math::Round" -lane '\$F[3] = 0 if (\$F[3] eq "." );  \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+250000-25000,\$F[1]+250000+24999,round(\$val*100)/100) if (\$val)' >${chrom}.win500k50k.tmp
  mapBed -a win500k100k.bed -b ${chrom}.GCcorrected.bedgraph -c 4,4 -o sum,count |perl -M"Math::Round" -lane '\$F[3] = 0 if (\$F[3] eq "." );  \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+250000-50000,\$F[1]+250000+49999,round(\$val*100)/100) if (\$val)' >${chrom}.win500k100k.tmp

  intersectBed -c -sorted -a ${chrom}.win500k50k.tmp  -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed  >${chrom}.w500ks50k.bedgraph
  intersectBed -c -sorted -a ${chrom}.win500k100k.tmp -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed  >${chrom}.w500ks100k.bedgraph
  """
  }

process mergeCoverageBedgraphs {
  echo true
  cpus 2
  memory '16 GB'
  time '2h'

  errorStrategy { 'retry' }
  maxRetries 1

  publishDir "${params.outdir}/bedgraph",  mode: 'copy', overwrite: true, pattern: '*bedgraph'
  publishDir "${params.outdir}/fig",       mode: 'copy', overwrite: true, pattern: '*pdf'
  publishDir "${params.outdir}/fig",       mode: 'copy', overwrite: true, pattern: '*png'

  input:
  path(gc)
  path(bg)

  output:
  path("*norm*.bedgraph",    emit: normBG)
  path("*coverage.bedgraph", emit: coverBG)
  path("*png",               emit: png)
  path("*pdf",               emit: pdf)

  script:
  """
  shuf -n 2000000 ${gc} >GCdata.tab

  R --vanilla <${params.accessoryDir}/scripts/R/plotGCcorrectionStats.R
  mv gcCorrection.png ${outName}.gcCorrection.png
  mv gcCorrection.pdf ${outName}.gcCorrection.pdf

  for step in 50 100; do
    sort -k1,1 -k2n,2n *".w500ks"\$step"k.bedgraph"  |perl -lane \'\$F[3] = 0 if (\$F[4]); print join("\\t",@F[0..3])\' >"${outName}.w500ks"\$step"k.coverage.bedgraph"
    perl ${params.accessoryDir}/scripts/normalizeBedgraph.pl "${outName}.w500ks"\$step"k.coverage.bedgraph" "${outName}.w500ks"\$step"k"

    mv "${outName}.w500ks"\$step"k.normByGenomeMean.penultimateBG"    "${outName}.w500ks"\$step"k.normByGenomeMean.bedgraph"
    mv "${outName}.w500ks"\$step"k.normByGenomeMedian.penultimateBG"  "${outName}.w500ks"\$step"k.normByGenomeMedian.bedgraph"
    mv "${outName}.w500ks"\$step"k.normByChromMean.penultimateBG"     "${outName}.w500ks"\$step"k.normByChromMean.bedgraph"
    mv "${outName}.w500ks"\$step"k.normByChromMedian.penultimateBG"   "${outName}.w500ks"\$step"k.normByChromMedian.bedgraph"

    perl ${params.accessoryDir}/scripts/convertToModellingBG.pl "${outName}.w500ks"\$step"k.normByGenomeMean.bedgraph"   ${params.genome_fai} >"${outName}.w500ks"\$step"k.normByGenomeMean.forModel.bedgraph"
    perl ${params.accessoryDir}/scripts/convertToModellingBG.pl "${outName}.w500ks"\$step"k.normByGenomeMedian.bedgraph" ${params.genome_fai} >"${outName}.w500ks"\$step"k.normByGenomeMedian.forModel.bedgraph"
    perl ${params.accessoryDir}/scripts/convertToModellingBG.pl "${outName}.w500ks"\$step"k.normByChromMean.bedgraph"    ${params.genome_fai} >"${outName}.w500ks"\$step"k.normByChromMean.forModel.bedgraph"
    perl ${params.accessoryDir}/scripts/convertToModellingBG.pl "${outName}.w500ks"\$step"k.normByChromMedian.bedgraph"  ${params.genome_fai} >"${outName}.w500ks"\$step"k.normByChromMedian.forModel.bedgraph"
  done
  """
  }

// these are for generating GC files
process generateMpileup {

  echo true
  cpus 2

  memory '32 GB'

  time { 9.hour * task.attempt}

  errorStrategy { 'retry' }
  maxRetries 2

  tag { chrom }

  //publishDir params.outdir,  mode: 'copy', overwrite: false

  input:
  path(bam)
  path(bai)
  val(chrom)

  output:
  //path("*.ALL.tab", emit: allGC)
  path("*.DS.tab",  emit: dsGC)

  script:
  def nregion;
  if (params.test & !params.fullChromTest){
    nregion=":1-15000000"
  }else{
    nregion=""
  }

  """
  tbam="${chrom}.tmp.bam"
  cbam="${chrom}.bam"

  samtools view -hb ${bam} ${chrom} >\$tbam
  samtools index \$tbam

  java -jar -Xmx8g \$PICARDJAR MarkDuplicates I=\$tbam O=\$cbam VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true AS=true M=metrics.tab TMP_DIR="\$TMPDIR"
  samtools index \$cbam

  readLen=`perl \$NXF_PIPEDIR/accessoryFiles/rtSeq/scripts/getClosestReadLength.pl \$tbam 50,150`
  pseudoReadDir=${params.pseudoReadBase}"/"\$readLen"bpReads_1bpStep/"

  samtools mpileup -r ${chrom}${nregion} -f ${params.genome_mask_fa} -q 31 \
                                   \$pseudoReadDir/${chrom}.bam \$cbam \
                                   |perl -lane \'print join("\\t",\$F[0],\$F[1]-1,\$F[1],\$F[3],\$F[6],(\$F[2] =~ /[GATC]/?1:0),(\$F[2] =~ /[GC]/?1:0))\' \
                                   >${chrom}.mpu

  grep -w ${chrom} ${params.genome_fai} >idx.fai
  bedtools makewindows -g idx.fai -w 101 -s 1 | intersectBed -v -sorted -a - -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed \
                                              |perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 101)\' >win101.bed

  nwin=`cat win101.bed |wc -l`
  nDS=`printf %1.0f \$(echo "\$nwin*0.02" |bc)`

  mapBed -a win101.bed -b ${chrom}.mpu -c 4,5,6,7 -o sum,sum,sum,sum |perl -pi -e 's/\\./NA/g' >${chrom}.ALL.tab

  ##Criteria :              NOT NA         && >0 cover &&   Fully mapable  &&  no repeat DNA
  perl -lane 'print \$_ if (\$F[4] !~ /NA/ &&  \$F[4]  && \$F[3] == ('\$readLen' * 101) && \$F[5] eq "101")' ${chrom}.ALL.tab >ok.tab

  shuf ok.tab |head -n \$nDS |sort -k1,1 -k2n,2n >${chrom}.DS.tab
  """
  }

process getGCstats {
  echo true
  cpus 2
  memory '32 GB'

  time '1h'

  errorStrategy { 'retry' }
  maxRetries 1

  publishDir "${params.outdir}/tables",  mode: 'copy', overwrite: true, pattern: '*DS.tab'
  publishDir "${params.outdir}/tables",  mode: 'copy', overwrite: true, pattern: '*table.tab'
  publishDir "${params.outdir}/figs",    mode: 'copy', overwrite: true, pattern: '*Before_v_After*'

  input:
  file(gctab)

  output:
  path("*.wholeGenome.DS.tab"       ,emit: allGCDS)
  path("*GCcorrectiontable.tab"     ,emit: correctionTab)
  path("*chromCorrectiontable.tab"  ,emit: chromMeans)
  path("*Before_v*.png"             ,emit: correctionPNG)
  path("*Before_v*.pdf"             ,emit: correctionPDF)

  script:
  """
  sort -k1,1 -k2n,2n *DS.tab |grep -P "^chr[0123456789XY]+" |cut -f1-3,5,7 >wholeGenome.DS.tab
  R --no-save --no-site-file --no-init-file --no-restore --silent --slave <${params.accessoryFiles}/scripts/R/getGCcorrectionFactors.R ||true

  cp wholeGenome.DS.tab ${outName}.wholeGenome.DS.tab

  mv chromCorrectiontable.tab ${outName}.chromCorrectiontable.tab
  mv GCcorrectiontable.tab    ${outName}.GCcorrectiontable.tab

  for f in `ls *_rtSeq_GCCorrection_Before_v_After.png`; do mv \$f ${outName}.rtSeq_GCCorrection_Before_v_After.png; done
  for f in `ls *_rtSeq_GCCorrection_Before_v_After.pdf`; do mv \$f ${outName}.rtSeq_GCCorrection_Before_v_After.pdf; done

  """
  }

//
// if (params.yeast){
// process correctByGCandSmoothYeast {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:300'
//
//   echo true
//   cpus 2
//
//   memory '16 GB'
//   module 'picard/2.9.2'
//   module 'bedtools/2.27.1'
//   module 'R/3.6.3'
//
//   time '2h'
//
//   errorStrategy { 'retry' }
//   maxRetries 3
//
//   tag { gctab }
//
//   input:
//   file(gctab)    from chrGCall
//   file(cFactors) from correctionTab
//
//   output:
//   file "*w*k*.bedgraph"       into allCSbedgraphs
//
//   shell:
//   def chrom    = gctab.name.replaceFirst(".ALL.tab","")
//
//   """
//   maxSim=`cut -f4 ${gctab} |shuf |head -n 20000 |sort -k1rn,1rn |head -n1`
//
//   ## DO GC NORMALIZATION
//   perl ${params.accessoryDir}/scripts/normalizeByGC.pl ${cFactors} ${gctab} \$maxSim
//
//   ## SMOOTHING STEP
//   grep -w ${chrom} ${params.genome_fai} >idx.fai
//
//   bedtools makewindows -g idx.fai -w 5000 -s 500    | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 5000)\' >win5k500bp.bed
//   bedtools makewindows -g idx.fai -w 5000 -s 1000   | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 5000)\' >win5k1k.bed
//   bedtools makewindows -g idx.fai -w 20000 -s 10000 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 20000)\' >win20k10k.bed
//
//   for bedg in ${chrom}.*bedgraph; do
//     bgNm=\${bedg/ALL./}
//
//     out5=\${bgNm/bedgraph/w5ks500bp.bedgraph}
//     out10=\${bgNm/bedgraph/w5ks1k.bedgraph}
//     out20=\${bgNm/bedgraph/w20ks10k.bedgraph}
//
//     mapBed -a win5k500bp.bed  -b \$bedg -c 4,4 -o sum,count |perl -lane '\$F[3] = 0 if (\$F[3] eq "." || \$F[4] < 500 );  \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+2500-250,\$F[1]+2500+249,\$val)'     >\$out5
//     mapBed -a win5k1k.bed     -b \$bedg -c 4,4 -o sum,count |perl -lane '\$F[3] = 0 if (\$F[3] eq "." || \$F[4] < 1000);  \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+2500-500,\$F[1]+2500+499,\$val)'     >\$out10
//     mapBed -a win20k10k.bed   -b \$bedg -c 4,4 -o sum,count |perl -lane '\$F[3] = 0 if (\$F[3] eq "." || \$F[4] < 10000); \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+10000-5000,\$F[1]+10000+4999,\$val)' >\$out20
//
//   done
//
//   """
//   }
//
// process mergeCSBGsYeast {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:300 --partition=norm'
//
//   echo true
//   cpus 2
//
//   memory '16 GB'
//   module 'bedtools/2.27.1'
//
//   time '2h'
//
//   errorStrategy { 'retry' }
//   maxRetries 3
//
//   //publishDir params.outBG, mode: 'copy', overwrite: true
//
//   input:
//   file(bg) from allCSbedgraphs.collect()
//
//   output:
//   file "*.penultimateBG"       into allPenultimateBGs
//
//   script:
//   //sort -k1,1 -k2n,2n *.GCcorrected.w500ks50k.bedgraph              | intersectBed -v -sorted -a - -b ${params.rtData}/blacklist/mm10_hotspot_blackList.bed |sort -k1,1 -k2n,2n >${outName}.w500ks50k.bedgraph
//   //sort -k1,1 -k2n,2n *.GCcorrected.w500ks100k.bedgraph             | intersectBed -v -sorted -a - -b ${params.rtData}/blacklist/mm10_hotspot_blackList.bed |sort -k1,1 -k2n,2n >${outName}.w500ks100k.bedgraph
//   """
//   sort -k1,1 -k2n,2n *.GCcorrected.w5ks500bp.bedgraph | intersectBed -c -sorted -a - -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed |sort -k1,1 -k2n,2n |perl -lane \'\$F[3] = 0 if (\$F[4]); print join("\\t",@F[0..3])\' >${outName}.w5ks500bp.penultimateBG
//   sort -k1,1 -k2n,2n *.GCcorrected.w5ks1k.bedgraph    | intersectBed -c -sorted -a - -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed |sort -k1,1 -k2n,2n |perl -lane \'\$F[3] = 0 if (\$F[4]); print join("\\t",@F[0..3])\' >${outName}.w5ks1k.penultimateBG
//   sort -k1,1 -k2n,2n *.GCcorrected.w20ks10k.bedgraph  | intersectBed -c -sorted -a - -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed |sort -k1,1 -k2n,2n |perl -lane \'\$F[3] = 0 if (\$F[4]); print join("\\t",@F[0..3])\' >${outName}.w2ms10k.penultimateBG
//
//   perl ${params.accessoryDir}/scripts/normalizeBedgraph.pl ${outName}.w5ks500bp.penultimateBG ${outName}.w5ks500bp
//   perl ${params.accessoryDir}/scripts/normalizeBedgraph.pl ${outName}.w5ks1k.penultimateBG    ${outName}.w5ks1k
//   perl ${params.accessoryDir}/scripts/normalizeBedgraph.pl ${outName}.w20kms10k.penultimateBG ${outName}.w20ks10k
//
//   """
//   }
//
// process finalizeBGsYeast {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:300 --partition=norm'
//
//   echo true
//   cpus 2
//
//   memory '16 GB'
//   module 'bedtools/2.27.1'
//
//   time '2h'
//
//   errorStrategy { 'retry' }
//   maxRetries 3
//
//   publishDir params.outBG, mode: 'copy', overwrite: true
//
//   input:
//   file(bg) from allPenultimateBGs.collect()
//
//   output:
//   file "*.bedgraph"       into allFinalBGs
//
//   script:
//   """
//   ## SMOOTHING STEP
//   cat ${params.genome_fai} >idx.fai
//
//   for bg in *.penultimateBG; do
//     outBG=\${bg/penultimateBG/bedgraph}
//
//     bedtools makewindows -g idx.fai -w 5000 -s 500    | perl -lane \'print join("\\t",\$F[0],\$F[1]+2250,\$F[2]-2251) if ((\$F[2]-\$F[1]) == 5000)\'  |sort -k1,1 -k2n,2n >win5k500bp.bed
//     bedtools makewindows -g idx.fai -w 5000 -s 1000   | perl -lane \'print join("\\t",\$F[0],\$F[1]+2000,\$F[2]-2001) if ((\$F[2]-\$F[1]) == 5000)\'  |sort -k1,1 -k2n,2n >win5k1k.bed
//     bedtools makewindows -g idx.fai -w 20000 -s 10000 | perl -lane \'print join("\\t",\$F[0],\$F[1]+5000,\$F[2]-5001) if ((\$F[2]-\$F[1]) == 20000)\' |sort -k1,1 -k2n,2n >win20k10k.bed
//
//     if [[ \$bg =~ "w5ks500bp" ]]; then
//       mapBed -a win5k500bp.bed -b \$bg -c 4 -o sum |perl -pi -e \'s/\\t\\./\\t0/g\' >\$outBG
//     fi
//
//     if [[ \$bg =~ "w5ks1k" ]]; then
//       mapBed -a win5k1k.bed -b \$bg -c 4 -o sum |perl -pi -e \'s/\\t\\./\\t0/g\' >\$outBG
//     fi
//
//     if [[ \$bg =~ "w20ks10k" ]]; then
//       mapBed -a win20k10k.bed -b \$bg -c 4 -o sum |perl -pi -e \'s/\\t\\./\\t0/g\' >\$outBG
//     fi
//   done
//   """
//   }
// }else{
// process correctByGCandSmooth {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:300'
//
//   echo true
//   cpus 2
//
//   memory '16 GB'
//   module 'picard/2.9.2'
//   module 'bedtools/2.27.1'
//   module 'R/3.6.3'
//
//   time '2h'
//
//   errorStrategy { 'retry' }
//   maxRetries 3
//
//   tag { gctab }
//
//   input:
//   file(gctab)    from chrGCall
//   file(cFactors) from correctionTab
//
//   output:
//   file "*00ks*0k*.bedgraph"       into allCSbedgraphs
//
//   shell:
//   def chrom    = gctab.name.replaceFirst(".ALL.tab","")
//
//   """
//   maxSim=`cut -f4 ${gctab} |shuf |head -n 20000 |sort -k1rn,1rn |head -n1`
//
//   ## DO GC NORMALIZATION
//   perl ${params.accessoryDir}/scripts/normalizeByGC.pl ${cFactors} ${gctab} \$maxSim
//
//   ## SMOOTHING STEP
//   grep -w ${chrom} ${params.genome_fai} >idx.fai
//
//   bedtools makewindows -g idx.fai -w 500000 -s 50000    | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500000)\' >win500k50k.bed
//   bedtools makewindows -g idx.fai -w 500000 -s 100000   | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500000)\' >win500k100k.bed
//   bedtools makewindows -g idx.fai -w 2000000 -s 1000000 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 2000000)\' >win2m1m.bed
//
//   for bedg in ${chrom}.*bedgraph; do
//     bgNm=\${bedg/ALL./}
//
//     out50=\${bgNm/bedgraph/w500ks50k.bedgraph}
//     out100=\${bgNm/bedgraph/w500ks100k.bedgraph}
//     out2000=\${bgNm/bedgraph/w2000ks1000k.bedgraph}
//
//     mapBed -a win500k50k.bed  -b \$bedg -c 4,4 -o sum,count |perl -lane '\$F[3] = 0 if (\$F[3] eq "." || \$F[4] < 5000 );  \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+250000-25000,\$F[1]+250000+24999,\$val)'     >\$out50
//     mapBed -a win500k100k.bed -b \$bedg -c 4,4 -o sum,count |perl -lane '\$F[3] = 0 if (\$F[3] eq "." || \$F[4] < 10000);  \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+250000-50000,\$F[1]+250000+49999,\$val)'     >\$out100
//     mapBed -a win2m1m.bed     -b \$bedg -c 4,4 -o sum,count |perl -lane '\$F[3] = 0 if (\$F[3] eq "." || \$F[4] < 100000); \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+1000000-500000,\$F[1]+1000000+499999,\$val)' >\$out2000
//
//   done
//
//   """
//   }
//
// process mergeCSBGs {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:300 --partition=norm'
//
//   echo true
//   cpus 2
//
//   memory '16 GB'
//   module 'bedtools/2.27.1'
//
//   time '2h'
//
//   errorStrategy { 'retry' }
//   maxRetries 3
//
//   //publishDir params.outBG, mode: 'copy', overwrite: true
//
//   input:
//   file(bg) from allCSbedgraphs.collect()
//
//   output:
//   file "*.penultimateBG"       into allPenultimateBGs
//
//   script:
//   //sort -k1,1 -k2n,2n *.GCcorrected.w500ks50k.bedgraph              | intersectBed -v -sorted -a - -b ${params.rtData}/blacklist/mm10_hotspot_blackList.bed |sort -k1,1 -k2n,2n >${outName}.w500ks50k.bedgraph
//   //sort -k1,1 -k2n,2n *.GCcorrected.w500ks100k.bedgraph             | intersectBed -v -sorted -a - -b ${params.rtData}/blacklist/mm10_hotspot_blackList.bed |sort -k1,1 -k2n,2n >${outName}.w500ks100k.bedgraph
//   """
//   sort -k1,1 -k2n,2n *.GCcorrected.w500ks50k.bedgraph    | intersectBed -c -sorted -a - -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed |sort -k1,1 -k2n,2n |perl -lane \'\$F[3] = 0 if (\$F[4]); print join("\\t",@F[0..3])\' >${outName}.w500ks50k.penultimateBG
//   sort -k1,1 -k2n,2n *.GCcorrected.w500ks100k.bedgraph   | intersectBed -c -sorted -a - -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed |sort -k1,1 -k2n,2n |perl -lane \'\$F[3] = 0 if (\$F[4]); print join("\\t",@F[0..3])\' >${outName}.w500ks100k.penultimateBG
//   sort -k1,1 -k2n,2n *.GCcorrected.w2000ks1000k.bedgraph | intersectBed -c -sorted -a - -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed |sort -k1,1 -k2n,2n |perl -lane \'\$F[3] = 0 if (\$F[4]); print join("\\t",@F[0..3])\' >${outName}.w2ms1m.penultimateBG
//
//   perl ${params.accessoryDir}/scripts/normalizeBedgraph.pl ${outName}.w500ks50k.penultimateBG  ${outName}.w500ks50k
//   perl ${params.accessoryDir}/scripts/normalizeBedgraph.pl ${outName}.w500ks100k.penultimateBG ${outName}.w500ks100k
//   perl ${params.accessoryDir}/scripts/normalizeBedgraph.pl ${outName}.w2ms1m.penultimateBG     ${outName}.w2ms1m
//
//   """
//   }
//
// process finalizeBGs {
//   scratch '/lscratch/$SLURM_JOBID'
//   clusterOptions ' --gres=lscratch:300 --partition=norm'
//
//   echo true
//   cpus 2
//
//   memory '16 GB'
//   module 'bedtools/2.27.1'
//
//   time '2h'
//
//   errorStrategy { 'retry' }
//   maxRetries 3
//
//   publishDir params.outBG, mode: 'copy', overwrite: true
//
//   input:
//   file(bg) from allPenultimateBGs.collect()
//
//   output:
//   file "*.bedgraph"       into allFinalBGs
//
//   script:
//   """
//   ## SMOOTHING STEP
//   cat ${params.genome_fai} >idx.fai
//
//   for bg in *.penultimateBG; do
//     outBG=\${bg/penultimateBG/bedgraph}
//
//     bedtools makewindows -g idx.fai -w 500000 -s 50000    | perl -lane \'print join("\\t",\$F[0],\$F[1]+225000,\$F[2]-225001) if ((\$F[2]-\$F[1]) == 500000)\'  |sort -k1,1 -k2n,2n >win500k50k.bed
//     bedtools makewindows -g idx.fai -w 500000 -s 100000   | perl -lane \'print join("\\t",\$F[0],\$F[1]+200000,\$F[2]-200001) if ((\$F[2]-\$F[1]) == 500000)\'  |sort -k1,1 -k2n,2n >win500k100k.bed
//     bedtools makewindows -g idx.fai -w 2000000 -s 1000000 | perl -lane \'print join("\\t",\$F[0],\$F[1]+500000,\$F[2]-500001) if ((\$F[2]-\$F[1]) == 2000000)\' |sort -k1,1 -k2n,2n >win2m1m.bed
//
//     if [[ \$bg =~ "w500ks50k" ]]; then
//       mapBed -a win500k50k.bed -b \$bg -c 4 -o sum |perl -pi -e \'s/\\t\\./\\t0/g\' >\$outBG
//     fi
//
//     if [[ \$bg =~ "w500ks100k" ]]; then
//       mapBed -a win500k100k.bed -b \$bg -c 4 -o sum |perl -pi -e \'s/\\t\\./\\t0/g\' >\$outBG
//     fi
//
//     if [[ \$bg =~ "w2ms1m" ]]; then
//       mapBed -a win2m1m.bed -b \$bg -c 4 -o sum |perl -pi -e \'s/\\t\\./\\t0/g\' >\$outBG
//     fi
//   done
//   """
//   }
