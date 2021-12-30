process ssdsAlign {

  time { fqR1.size()<300000000? 3.hour : 3.hour * fqR1.size()/300000000 * task.attempt }

  tag {fqR1}

  input:
  tuple file(fqR1),file(fqR2)

  output:
  path('*_split_*bam')

  script:
  // generate random filenames
  def nm = new Random().with {(1..12).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}
  def tmp = 'tmp_' + nm
  def tmpR1fq = 'tmp_' + nm + '.R1.fq'
  def tmpR2fq = 'tmp_' + nm + '.R2.fq'
  def tmpR1sai = 'tmp_' + nm + '.R1.sai'
  def tmpR2sai = 'tmp_' + nm + '.R2.sai'
  def bamOut = params.name + '_split_' + nm + '.bam'
  def bamAll = 'all_' + nm + '.bam'

  println("R1 len: ${params.r1Len}")

  if (params.original){
    """
	  fastx_trimmer -Q 33 -f 1 -l ${params.r1Len} -i ${fqR1} -o ${tmpR1fq}
	  fastx_trimmer -Q 33 -f 1 -l ${params.r2Len} -i ${fqR2} -o ${tmpR2fq}

	  bwa_0.7.12 aln \
	    -t ${task.cpus} \
	    ${params.genome_bwaidx} \
	    ${tmpR1fq} >${tmpR1sai}

	  bwa_ra_0.7.12 aln \
	    -t ${task.cpus} \
	    ${params.genome_bwaidx} \
	    ${tmpR2fq} >${tmpR2sai}

	  bwa_0.7.12 sampe \
	    ${params.genome_bwaidx} \
	    ${tmpR1sai} ${tmpR2sai} \
	    ${tmpR1fq} ${tmpR2fq} >${tmp}".unsorted.sam"

	  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SamFormatConverter \
                 I=${tmp}.unsorted.sam \
                 O=${tmp}.unsorted.tmpbam \
                 TMP_DIR="\$TMPDIR" \
                 VALIDATION_STRINGENCY=LENIENT

	  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam \
                 I=${tmp}.unsorted.tmpbam \
                 O=${bamOut} \
                 SO=coordinate \
                 TMP_DIR="\$TMPDIR" \
                 VALIDATION_STRINGENCY=LENIENT

    samtools index ${bamOut}
    """
  }else{
    """
    fqLen=`cat ${fqR1} |sed -n '2~4p' |perl -lane 'print length(\$_)' |head -n 10000 |sort -k1n,1n |tail -n1`

    if [ \$fqLen -ge 99 ]; then
      echo "FASTQ Max Length = \$fqLen : Using minimap2 for initial alignment ..."
      minimap2 -ax sr -I 12g \
        -t ${task.cpus} \
        ${params.genome_mm2idx} \
        ${fqR1} ${fqR2} >${tmp}.unsorted.sam
    else
      echo "FASTQ Max Length = \$fqLen : Using bwa mem for initial alignment ..."
      bwa_0.7.12 mem \
        -t ${task.cpus} \
        ${params.genome_bwaidx} \
        ${fqR1} ${fqR2} >${tmp}.unsorted.sam
    fi

    java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SamFormatConverter \
              I=${tmp}.unsorted.sam \
              O=${tmp}.unsorted.tmpbam \
              TMP_DIR="\$TMPDIR" \
              VALIDATION_STRINGENCY=LENIENT

    java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam \
              I=${tmp}.unsorted.tmpbam \
              O=${tmp}.Qsorted.sam \
              SO=queryname \
              TMP_DIR="\$TMPDIR" \
              VALIDATION_STRINGENCY=LENIENT

    flagReadPairsWithMultipleSupplementaryMappings.pl ${tmp}.Qsorted.sam ${tmp}.ok.sam

    java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam \
              I=${tmp}.ok.sam \
              O=${bamAll} \
              SO=coordinate \
              TMP_DIR="\$TMPDIR" \
              VALIDATION_STRINGENCY=LENIENT

    samtools index ${bamAll}

    samtools view -hb ${bamAll} >${bamOut}
    samtools index ${bamOut}
    """
    }
  }

process mergeBAMssds {

  publishDir "${params.outdir}/bam",     mode: 'copy', overwrite: true, pattern: '*unparsed.bam*'
  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: '*MDmetrics.txt*'

  tag {params.name}

  input:
  path(bams)

  output:
  tuple path('*unparsed.bam'),path('*unparsed.bam.bai'), emit: bam
  tuple path('*ignments.bam'),path('*ignments.bam.bai'), emit: suppbam, optional: true
  path('*MDmetrics.txt', emit: mdmetrics)

  script:
  // get INPUT files as string
  def input_args = bams.collect{ "I=$it" }.join(" ")
  def name = "${params.name}.${params.genome}"
  """
  java -jar -Xmx32g \$PICARDJAR MergeSamFiles \
                 ${input_args} \
                 O=allREADS.bam \
                 AS=true \
                 SO=coordinate \
                 TMP_DIR="\$TMPDIR" \
                 VALIDATION_STRINGENCY=LENIENT

  java -jar -Xmx32g \$PICARDJAR MarkDuplicatesWithMateCigar \
               I=allREADS.bam \
               O=${name}.allReads.bam \
               PG=Picard2.9.2_MarkDuplicates \
               M=${name}.MDmetrics.txt \
               MINIMUM_DISTANCE=400 \
               TMP_DIR="\$TMPDIR" \
		   CREATE_INDEX=false \
		   ASSUME_SORT_ORDER=coordinate \
       VALIDATION_STRINGENCY=LENIENT

  samtools index ${name}.allReads.bam

  samtools view -F 2048 -hb ${name}.allReads.bam > ${name}.SSDSunparsed.bam
  samtools index ${name}.SSDSunparsed.bam

  #samtools view -f 2 -hb ${name}.SSDSall.bam >${name}.SSDSunparsed_andPaired.bam
  #samtools index ${name}.SSDSunparsed_andPaired.bam

  samtools view -f 2048 -hb ${name}.allReads.bam > ${name}.SSDSunparsed.suppAlignments.bam
  samtools index ${name}.SSDSunparsed.suppAlignments.bam

  """
  }

process parseITRs {

  tag {bam}

  cpus 4
  memory '12 GB'
  time { 6.hour * task.attempt}

  errorStrategy 'retry'
  maxRetries 1

  input:
  path(bam)

  output:
  path('*.bam', emit: pBam)
  path('*.bai', emit: pBai)
  path('*.bed', emit: pBed)

  script:
  def parseScript = params.original ? 'ITR_id_v3.pl' : 'ITR_id_v3_long.pl';
  def name = bam.name.replaceFirst(".bam",".pairs.bam")
  """
  samtools view -f 2 -hb ${bam} >${name}
  samtools index ${name}

  ${parseScript} ${name} ${params.genome} 2>/dev/null

  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${name}.ssDNA_type1.bed  -o ${name}.ssDNA_type1.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${name}.ssDNA_type2.bed  -o ${name}.ssDNA_type2.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${name}.dsDNA.bed        -o ${name}.dsDNA.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${name}.dsDNA_strict.bed -o ${name}.dsDNA_strict.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${name}.unclassified.bed -o ${name}.unclassified.bed

  samtools view -H ${name}   >header.txt
  echo -e '@PG\\tID:SSDSpipeline PN:SSDSpipeline VN:2.0_nextflowDSL2' >>header.txt

  cat header.txt ${name}.ssDNA_type1.sam  >${name}.ssDNA_type1.RH.sam
  cat header.txt ${name}.ssDNA_type2.sam  >${name}.ssDNA_type2.RH.sam
  cat header.txt ${name}.dsDNA.sam        >${name}.dsDNA.RH.sam
  cat header.txt ${name}.dsDNA_strict.sam >${name}.dsDNA_strict.RH.sam
  cat header.txt ${name}.unclassified.sam >${name}.unclassified.RH.sam

  samtools view -Shb ${name}.ssDNA_type1.RH.sam  >${name}.ssDNA_type1.US.bam
  samtools view -Shb ${name}.ssDNA_type2.RH.sam  >${name}.ssDNA_type2.US.bam
  samtools view -Shb ${name}.dsDNA.RH.sam        >${name}.dsDNA.US.bam
  samtools view -Shb ${name}.dsDNA_strict.RH.sam >${name}.dsDNA_strict.US.bam
  samtools view -Shb ${name}.unclassified.RH.sam >${name}.unclassified.US.bam

  java -jar -Xmx8g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" I=${name}.ssDNA_type1.US.bam  O=${name}.ssDNA_type1.bam  SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar -Xmx8g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" I=${name}.ssDNA_type2.US.bam  O=${name}.ssDNA_type2.bam  SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar -Xmx8g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" I=${name}.dsDNA.US.bam        O=${name}.dsDNA.bam        SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar -Xmx8g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" I=${name}.dsDNA_strict.US.bam O=${name}.dsDNA_strict.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar -Xmx8g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" I=${name}.unclassified.US.bam O=${name}.unclassified.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT

  samtools index ${name}.ssDNA_type1.bam
  samtools index ${name}.ssDNA_type2.bam
  samtools index ${name}.dsDNA.bam
  samtools index ${name}.dsDNA_strict.bam
  samtools index ${name}.unclassified.bam
  """
  }

process gatherITROutputs {

  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: '*txt'
  publishDir "${params.outdir}/bam",     mode: 'copy', overwrite: true, pattern: '*bam*'
  publishDir "${params.outdir}/bed",     mode: 'copy', overwrite: true, pattern: '*bed*'

  tag {bam}

  time { bam.size() < 100 ? 2.hour : 2.hour * (1 + bam.size()/100 * task.attempt) }

  input:
  path(bam)
  path(bed)
  val(type)

  output:
  tuple(path('*bam'), path('*bai'), emit: bambai)
  tuple(path('*bam'), path('*bai'), path('*bed'), path('*metrics.txt'), emit: itrFinal)
  path('*.txt*', emit: report)
  path('*.bam*', emit: bam)
  path('*.bed',  emit: bed)

  script:
  def name="${params.name}.${params.genome}.${type}"
  """
  java -jar -Xmx8g \$PICARDJAR MergeSamFiles TMP_DIR="\$TMPDIR" O=${name}.US.BAM `ls *${type}.bam | sed 's/.*\$/I=& /'`
  java -jar -Xmx8g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" I=${name}.US.BAM   O=${name}.S.BAM   SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar -Xmx8g \$PICARDJAR MarkDuplicatesWithMateCigar TMP_DIR="\$TMPDIR" I=${name}.S.BAM  O=${name}.bam  PG=Picard2.9.2_MarkDuplicates M=${name}.MDmetrics.txt  CREATE_INDEX=false VALIDATION_STRINGENCY=LENIENT

  samtools index ${name}.bam

  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 -m *${type}.bed  >${name}.bed
  rm -f ${name}.US.BAM
  rm -f ${name}.S.BAM
    """
  }

process makeFRBWssds {

  publishDir "${params.outdir}/bigwig", mode: 'copy', overwrite: true, pattern: '*bigwig*'

  cpus 2
  memory '12 GB'

  tag {bam}
  time { bam.size()< 2000000000 ? 1.hour * task.attempt : 1.hour * bam.size()/2000000000 * task.attempt }

  errorStrategy { 'retry' }
  maxRetries 1

  input:
  tuple path(bam), path(idx)

  output:
  path('*.bigwig', emit: bw)

  script:
  iName = bam.name.replaceAll(/.bam/,".out")

  """
  if [[ `samtools view -F 4 ${bam} |head -n 5001 |wc -l` -gt 5000 ]]; then
    idx="${params.genome_fai}";

    bedtools makewindows  -g \$idx -w 1000 -s 100 |sort -k1,1 -k2n,2n |perl -lane 'print \$_ if ((\$F[2]-\$F[1]) == 1000)' >win.bed
    java -jar -Xmx8g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" I=${bam} O=qSort.bam SO=queryname VALIDATION_STRINGENCY=LENIENT
    bedtools bamtobed -mate1 -i qSort.bam -bedpe >frags.bedpe

    perl -lane 'print join("\\t", \$F[0], \$F[1], \$F[5], @F[6..8]) if (\$F[8] eq "+")' frags.bedpe |sort -k1,1 -k2n,2n >frags.POS.bed
    perl -lane 'print join("\\t", \$F[0], \$F[4], \$F[2], @F[6..8]) if (\$F[8] eq "-")' frags.bedpe |sort -k1,1 -k2n,2n >frags.NEG.bed

    mapBed -a win.bed -b frags.POS.bed -c 1 -o count |perl -lane 'use Math::Round; \$p=round((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$p-49,\$p+50,\$F[3])' >pos.bg
    mapBed -a win.bed -b frags.NEG.bed -c 1 -o count |perl -lane 'use Math::Round; \$p=round((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$p-49,\$p+50,\$F[3])' >neg.bg

    paste pos.bg neg.bg |perl -lane '\$F[3]+=0.5; \$F[7]+=0.5; \$v=(\$F[3]/\$F[7]); print join("\\t",@F[0..2],(log(\$v)/log(2)))' >fr.bg
    paste pos.bg neg.bg |perl -lane '\$v=(\$F[3]+\$F[7]); print join("\\t",@F[0..2],\$v)' >tot.bg

    bedGraphToBigWig pos.bg \$idx ${iName}.F.bigwig
    bedGraphToBigWig neg.bg \$idx ${iName}.R.bigwig
    bedGraphToBigWig fr.bg  \$idx ${iName}.FR.bigwig
    bedGraphToBigWig tot.bg \$idx ${iName}.Tot.bigwig
  else
    touch "EMPTY_TOOFEWREADS_${iName}.F.bigwig
    touch "EMPTY_TOOFEWREADS_${iName}.R.bigwig
    touch "EMPTY_TOOFEWREADS_${iName}.FR.bigwig
    touch "EMPTY_TOOFEWREADS_${iName}.Tot.bigwig
  fi
  """
  }

process makeSSreport {

  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

  tag {bam}

  input:
  tuple path(bam), path(bai)
  path(ssfiles)

  output:
  path('*SSDSreport*', emit: report)

  script:
  """
  echo "makeSSMultiQCReport_nextFlow.pl ${bam} "`ls *bed`" --g ${params.genome}"
  makeSSMultiQCReport_nextFlow.pl ${bam} `ls *bed` --g ${params.genome}

  """
  }

process multiQCssds {

  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

  tag {reports}

  input:
  path(reports)

  output:
  path('*ultiQC*', emit: mqcReport)

  script:
  """
  multiqc -f -m ssds -n ${params.name}.multiQC .
  """
  }
