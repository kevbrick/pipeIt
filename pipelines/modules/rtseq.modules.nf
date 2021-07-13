process getCoverage {

  tag { chrom }

  input:
  tuple(path(bam), path(bai))
  path(gctab)
  val(chrom)

  output:
  path("*.w*0k.bedgraph", emit: bg)
  path("*.GCdata.tab",    emit: gcData)

  script:
  
  def win101="${chrom}.win101ALL.bed"
  """
  tbam="${chrom}.tmp.bam"
  s1bam="${chrom}.s1.bam"
  s2bam="${chrom}.s2.bam"
  s3bam="${chrom}.s3.bam"
  s4bam="${chrom}.s4.bam"
  cbam="${chrom}.bam"

  samtools view -hb -q 30 ${bam} ${chrom} >\$tbam
  samtools index \$tbam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam I=\$tbam O=\$s1bam VALIDATION_STRINGENCY=LENIENT SO=queryname TMP_DIR="\$TMPDIR"

  nPE=`samtools view -h \$tbam |head -n 100000 |samtools view -c -f 1 -S /dev/stdin`

  if [ "\$nPE" -eq "0" ]; then
    ln -s \$s1bam \$s3bam
  else
    java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR FixMateInformation I=\$s1bam O=\$s2bam VALIDATION_STRINGENCY=LENIENT SO=queryname AS=true TMP_DIR="\$TMPDIR"
    samtools view -hb -f 2 \$s2bam >\$s3bam
  fi

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MarkDuplicates I=\$s3bam O=\$s4bam VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true M=metrics.tab TMP_DIR="\$TMPDIR"
  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam I=\$s4bam O=\$cbam VALIDATION_STRINGENCY=LENIENT SO=coordinate TMP_DIR="\$TMPDIR"
  samtools index \$cbam

  readLen=`perl ${params.rtData}/scripts/getClosestReadLength.pl \$tbam 50,150`
  #win101="${params.rtData}/genomewin/${chrom}.win101.bed"
  #ps101="${params.rtData}/genomewin/${chrom}.psr"\$readLen".tab"
  #gc101="${params.rtData}/genomewin/${chrom}.win101.GC.tab"

  ### NEW ###
  grep -w ${chrom} ${params.genome_fai} >chrom.fai
  
  bedtools makewindows -g chrom.fai -w 101 -s 1 |
    perl -lane 'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 101)' >${win101}
    
  bedtools nuc -fi ${params.genome_fasta} -bed ${win101} -C | grep ^chr | \
                   cut -f1-3,5 | perl -M"Math::Round" -pi -e 's/(0\\.\\d\\d+)/round(\$1*100)/e' |cut -f4 >GC.tab
  
  cp "${params.rtData}/mapability/${params.genome}/${chrom}.mapability_"\$readLen"bp.tab.gz" ./mapability.txt.gz
  gunzip mapability.txt.gz
      
  intersectBed -a ${win101} -b \$cbam -c -sorted -g ${params.genome_fai} |perl -pi -e 's/\\./0/g' >${chrom}.cover.tab

  paste ${chrom}.cover.tab mapability.txt GC.tab |perl -pi -e 's/\\s\\./0(\\s|\$)/g' >${chrom}.tmp
  
  intersectBed -v -sorted -a ${chrom}.tmp -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed >${chrom}.tab
  ### END NEW ###
  
  ## DO GC NORMALIZATION
  expectedSim=`echo \$readLen + 100 |bc`
  mapabilityLowerLim=`echo "(\$readLen + 100)*0.66" |bc`
  mapabilityUpperLim=`echo "(\$readLen + 100)*1.2" |bc`

  perl ${params.rtData}/scripts/normalizeBy2NGC.pl ${gctab} ${chrom}.tab \$expectedSim

  shuf GCData.tab |head -n 1000000 |sort -k1rn,1rn >${chrom}.GCdata.tab

  perl -lane 'print join("\\t",@F[0..2],\$F[4]) if (\$F[4] > '\$mapabilityLowerLim' && \$F[4] <  '\$mapabilityUpperLim' )' ${chrom}.tab >${chrom}.mapability.tab
  nCov=`echo "500000 * ${params.covPC} /100" |bc`
  
  bedtools makewindows -g chrom.fai -w 500000 -s 50000  | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500000)\' >win500k50k.init.bed
  bedtools makewindows -g chrom.fai -w 500000 -s 100000 | perl -lane \'print join("\\t",@F) unless ((\$F[2]-\$F[1]) != 500000)\' >win500k100k.init.bed

  mapBed -a win500k50k.init.bed  -b ${chrom}.mapability.tab -c 4 -o count |perl -lane 'print join("\\t",@F[0..2]) if (\$F[3] > int(500000*'${params.covPC}'/100))' >win500k50k.bed
  mapBed -a win500k100k.init.bed -b ${chrom}.mapability.tab -c 4 -o count |perl -lane 'print join("\\t",@F[0..2]) if (\$F[3] > int(500000*'${params.covPC}'/100))' >win500k100k.bed

  mapBed -a win500k50k.bed  -b ${chrom}.GCcorrected.bedgraph -c 4,4 -o sum,count |perl -M"Math::Round" -lane '\$F[3] = 0 if (\$F[3] eq "." );  \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+250000-25000,\$F[1]+250000+24999,round(\$val*100)/100) if (\$val)' >${chrom}.win500k50k.tmp
  mapBed -a win500k100k.bed -b ${chrom}.GCcorrected.bedgraph -c 4,4 -o sum,count |perl -M"Math::Round" -lane '\$F[3] = 0 if (\$F[3] eq "." );  \$val=\$F[4]?\$F[3]/\$F[4]:0; print join("\\t",\$F[0],\$F[1]+250000-50000,\$F[1]+250000+49999,round(\$val*100)/100) if (\$val)' >${chrom}.win500k100k.tmp

  intersectBed -c -sorted -a ${chrom}.win500k50k.tmp  -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed  >${chrom}.w500ks50k.bedgraph
  intersectBed -c -sorted -a ${chrom}.win500k100k.tmp -b ${params.rtData}/blacklist/${params.genome}.blacklist.bed  >${chrom}.w500ks100k.bedgraph
  """
  }

process mergeCoverageBedgraphs {

  publishDir "${params.outdir}/bedgraph/RTmodelling/raw"                ,  mode: 'copy', overwrite: true, pattern: '*age.forModel.bedgraph'
  publishDir "${params.outdir}/bedgraph/RTmodelling/medianNorm/byChr"   ,  mode: 'copy', overwrite: true, pattern: '*ChromMedian.forModel.bedgraph'
  publishDir "${params.outdir}/bedgraph/RTmodelling/meanNorm/byChr"     ,  mode: 'copy', overwrite: true, pattern: '*ChromMean.forModel.bedgraph'
  publishDir "${params.outdir}/bedgraph/RTmodelling/medianNorm/byGenome",  mode: 'copy', overwrite: true, pattern: '*GenomeMedian.forModel.bedgraph'
  publishDir "${params.outdir}/bedgraph/RTmodelling/meanNorm/byGenome"  ,  mode: 'copy', overwrite: true, pattern: '*GenomeMean.forModel.bedgraph'
  publishDir "${params.outdir}/bedgraph/raw"                            ,  mode: 'copy', overwrite: true, pattern: '*age.bedgraph'
  publishDir "${params.outdir}/bedgraph/medianNorm/byChr"               ,  mode: 'copy', overwrite: true, pattern: '*ChromMedian.bedgraph'
  publishDir "${params.outdir}/bedgraph/meanNorm/byChr"                 ,  mode: 'copy', overwrite: true, pattern: '*ChromMean.bedgraph'  
  publishDir "${params.outdir}/bedgraph/medianNorm/byGenome"            ,  mode: 'copy', overwrite: true, pattern: '*GenomeMedian.bedgraph'
  publishDir "${params.outdir}/bedgraph/meanNorm/byGenome"              ,  mode: 'copy', overwrite: true, pattern: '*GenomeMean.bedgraph'
  publishDir "${params.outdir}/bedgraph"                                ,  mode: 'copy', overwrite: true, pattern: '*w500ks50k.normByGenomeMedian.forModel.bedgraph'
  publishDir "${params.outdir}/fig"                                     ,  mode: 'copy', overwrite: true, pattern: '*pdf'
  publishDir "${params.outdir}/fig"                                     ,  mode: 'copy', overwrite: true, pattern: '*png'

  input:
  val(outName)
  path(gc)
  path(bg)

  output:
  path("*norm*.bedgraph",    emit: normBG)
  path("*coverage.bedgraph", emit: coverBG)
  path("*png",               emit: png)
  path("*pdf",               emit: pdf)

  script:
  """
  shuf -n 2000000 ${bg} |awk NF >GCdata.tab
  
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
process getGCcorrFile{
  
  publishDir "${params.outdir}/gcCorrectionFile",  mode: 'copy', overwrite: false

  output:
  path("gcCorrectionTable.tab",  emit: tab)

  script:
  """
  if [ -f "${params.gcCorrection}" ]; then
    cp ${params.gcCorrection} gcCorrectionTable.tab
  else
    cp /usr/local/rtseq/gcCorrectionFiles/${params.gcCorrection.toUpperCase()}.tab gcCorrectionTable.tab
  fi
  """  
  }

process generateMpileup {

  tag { chrom }

  //publishDir params.outdir,  mode: 'copy', overwrite: false

  input:
  tuple(path(bam),path(bai))
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

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MarkDuplicates I=\$tbam O=\$cbam \
                     VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true AS=true M=metrics.tab TMP_DIR="\$TMPDIR"
  samtools index \$cbam

  readLen=`perl ${params.accessoryDir}/rtSeq/scripts/getClosestReadLength.pl \$tbam 50,150`
  pseudoReadDir=${params.pseudoReadBase}"/"\$readLen"bpReads_1bpStep/"
  
  genomeCoverageBed -ibam \$cbam -d -g >coverage.gcb.tab
  
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
