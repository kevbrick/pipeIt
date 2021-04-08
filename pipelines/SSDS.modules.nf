process ssdsAlign {
  cpus 4
  memory 6

  time { fqR1.size()<100000000? 1.hour : 1.hour * fqR1.size()/100000000 * task.attempt }

  input:
  tuple file(fqR1),file(fqR2)

  output:
  path('bar*bam')

  script:
  // generate random filenames
  def nm = new Random().with {(1..30).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}
  def tmp = 'tmp_' + nm
  def tmpR1fq  = 'tmp_' + nm + 'A.R1.fq'
  def tmpR2fq  = 'tmp_' + nm + 'A.R2.fq'
  def tmpR1fqB = 'tmp_' + nm + 'B.R1.fq'
  def tmpR2fqB = 'tmp_' + nm + 'B.R2.fq'
  def tmpR1sai = 'tmp_' + nm + '.R1.sai'
  def tmpR2sai = 'tmp_' + nm + '.R2.sai'
  def bamOut = 'bar_' + nm + '.bam'
  def bamAll = 'all_' + nm + '.bam'

  println("R1 len: ${params.r1Len}")

  if (params.long){

    """
    minimap2 -ax sr \
      -t ${task.cpus} \
      ${params.genome_mm2idx} \
      ${fqR1} ${fqR2} >${tmp}.unsorted.sam

    java -jar \$PICARDJAR SamFormatConverter \
              I=${tmp}.unsorted.sam \
              O=${tmp}.unsorted.tmpbam \
              VALIDATION_STRINGENCY=LENIENT

    java -jar \$PICARDJAR SortSam \
              I=${tmp}.unsorted.tmpbam \
              O=${tmp}.Qsorted.sam \
              SO=queryname \
              VALIDATION_STRINGENCY=LENIENT

    flagReadPairsWithMultipleSupplementaryMappings.pl ${tmp}.Qsorted.sam ${tmp}.ok.sam

    java -jar \$PICARDJAR SortSam \
              I=${tmp}.ok.sam \
              O=${bamAll} \
              SO=coordinate \
              VALIDATION_STRINGENCY=LENIENT

    samtools index ${bamAll}

    samtools view -f 2 -hb ${bamAll} >${bamOut}
    samtools index ${bamOut}
    """
  }else{
    """
	  fastx_trimmer -Q 33 -f 1 -l ${params.r1Len} -i ${fqR1} -o ${tmpR1fq}
	  fastx_trimmer -Q 33 -f 1 -l ${params.r2Len} -i ${fqR2} -o ${tmpR2fq}

    ## Remove paired reads where either R1 / R2 is shorter than 32bp
    ## Helps to minimaize the likelihood of a very rare, but very problematic
    ## bwa bug. This bug occurred (rarely) when short reads had a long ITR ...
    ## but killed the pipeline when it did occur.

    filter_paired_end_FQ_by_size --fqr1 ${tmpR1fq} --fqr2 ${tmpR2fq} --size 31 --outr1 ${tmpR1fqB} --outr2 ${tmpR2fqB}

    rm ${tmpR1fq} ${tmpR2fq}

    bwa_0.7.12 aln \
	    -t ${task.cpus} \
	    ${params.genome_bwaidx} \
	    ${tmpR1fqB} >${tmpR1sai}

	  bwa_ra_0.7.12 aln \
	    -t ${task.cpus} \
	    ${params.genome_bwaidx} \
	    ${tmpR2fqB} >${tmpR2sai}

	  bwa_0.7.12 sampe \
	    ${params.genome_bwaidx} \
	    ${tmpR1sai} ${tmpR2sai} \
	    ${tmpR1fqB} ${tmpR2fqB} >${tmp}".unsorted.sam"

    rm ${tmpR1fqB} ${tmpR2fqB}

	  java -jar \$PICARDJAR SamFormatConverter \
                 I=${tmp}.unsorted.sam \
                 O=${tmp}.unsorted.tmpbam \
                 VALIDATION_STRINGENCY=LENIENT

	  java -jar \$PICARDJAR SortSam \
                 I=${tmp}.unsorted.tmpbam \
                 O=${bamOut} \
                 SO=coordinate \
                 VALIDATION_STRINGENCY=LENIENT

    samtools index ${bamOut}
    """

  }
  }

process mergeBAMssds {

  publishDir "${params.outdir}/bam",     mode: 'copy', overwrite: true, pattern: '*unparsed*bam*'
  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: '*MDmetrics.txt*'

  cpus 4
  memory '4g'

  time { 4.hour * task.attempt }

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
  java -jar \$PICARDJAR MergeSamFiles \
                 ${input_args} \
                 O=allREADS.bam \
                 AS=true \
                 SO=coordinate \
                 VALIDATION_STRINGENCY=LENIENT

  samtools view -f 2 -F 2048 -hb allREADS.bam    >allREADS.ok.bam

  java -jar \$PICARDJAR MarkDuplicatesWithMateCigar \
                 I=allREADS.ok.bam \
                 O=${name}.SSDSunparsed.bam \
                 PG=Picard2.9.2_MarkDuplicates \
                 M=${name}.MDmetrics.txt \
                 MINIMUM_DISTANCE=400 \
			   CREATE_INDEX=false \
			   ASSUME_SORT_ORDER=coordinate \
         VALIDATION_STRINGENCY=LENIENT

  samtools index ${name}.SSDSunparsed.bam

  samtools view -f 2048 -hb allREADS.bam >allREADS.supp.bam

  java -jar \$PICARDJAR MarkDuplicatesWithMateCigar \
                 I=allREADS.supp.bam \
                 O=${name}.SSDSunparsed.suppAlignments.bam \
                 PG=Picard2.9.2_MarkDuplicates \
                 M=${name}.suppAlignments.MDmetrics.txt \
                 MINIMUM_DISTANCE=400 \
			   CREATE_INDEX=false \
			   ASSUME_SORT_ORDER=coordinate \
         VALIDATION_STRINGENCY=LENIENT

  samtools index ${name}.SSDSunparsed.suppAlignments.bam
  """
  }

process parseITRs {

  cpus 4
  memory '6g'
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
  def parseScript = params.long ? 'ITR_id_v3_long.pl' : 'ITR_id_v3.pl'
  """
  ${parseScript} ${bam} ${params.genome}

  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.ssDNA_type1.bed  -o ${bam}.ssDNA_type1.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.ssDNA_type2.bed  -o ${bam}.ssDNA_type2.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.dsDNA.bed        -o ${bam}.dsDNA.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.dsDNA_strict.bed -o ${bam}.dsDNA_strict.bed
  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 ${bam}.unclassified.bed -o ${bam}.unclassified.bed

  samtools view -H ${bam}   >header.txt
  echo '@SSDSpipeline_2.0' >>header.txt

  cat header.txt ${bam}.ssDNA_type1.sam  >${bam}.ssDNA_type1.RH.sam
  cat header.txt ${bam}.ssDNA_type2.sam  >${bam}.ssDNA_type2.RH.sam
  cat header.txt ${bam}.dsDNA.sam        >${bam}.dsDNA.RH.sam
  cat header.txt ${bam}.dsDNA_strict.sam >${bam}.dsDNA_strict.RH.sam
  cat header.txt ${bam}.unclassified.sam >${bam}.unclassified.RH.sam

  samtools view -Shb ${bam}.ssDNA_type1.RH.sam  >${bam}.ssDNA_type1.US.bam
  samtools view -Shb ${bam}.ssDNA_type2.RH.sam  >${bam}.ssDNA_type2.US.bam
  samtools view -Shb ${bam}.dsDNA.RH.sam        >${bam}.dsDNA.US.bam
  samtools view -Shb ${bam}.dsDNA_strict.RH.sam >${bam}.dsDNA_strict.US.bam
  samtools view -Shb ${bam}.unclassified.RH.sam >${bam}.unclassified.US.bam

  java -jar \$PICARDJAR SortSam I=${bam}.ssDNA_type1.US.bam  O=${bam}.ssDNA_type1.bam  SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar \$PICARDJAR SortSam I=${bam}.ssDNA_type2.US.bam  O=${bam}.ssDNA_type2.bam  SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar \$PICARDJAR SortSam I=${bam}.dsDNA.US.bam        O=${bam}.dsDNA.bam        SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar \$PICARDJAR SortSam I=${bam}.dsDNA_strict.US.bam O=${bam}.dsDNA_strict.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar \$PICARDJAR SortSam I=${bam}.unclassified.US.bam O=${bam}.unclassified.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}.ssDNA_type1.bam
  samtools index ${bam}.ssDNA_type2.bam
  samtools index ${bam}.dsDNA.bam
  samtools index ${bam}.dsDNA_strict.bam
  samtools index ${bam}.unclassified.bam
  """
  }

process gatherITROutputs {

  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: '*txt'
  publishDir "${params.outdir}/bam",     mode: 'copy', overwrite: true, pattern: '*bam*'
  publishDir "${params.outdir}/bed",     mode: 'copy', overwrite: true, pattern: '*bed*'

  cpus 4
  memory '6g'

  time { bam.size() < 2000000000 ? 1.hour : 1.hour * bam.size()/2000000000 * task.attempt }

  input:
  path(bam)
  path(bed)
  val(type)

  output:
  tuple(path('*bam'), path('*bai'), path('*bed'), path('*metrics.txt'), emit: itrFinal)
  path('*.txt*', emit: report)
  path('*.bam*', emit: bam)
  path('*.bed',  emit: bed)

  script:
  def name="${params.name}.${params.genome}.${type}"
  """
  java -jar \$PICARDJAR MergeSamFiles O=${name}.US.BAM `ls *${type}.bam | sed 's/.*\$/I=& /'`
  java -jar \$PICARDJAR SortSam I=${name}.US.BAM   O=${name}.S.BAM   SO=coordinate VALIDATION_STRINGENCY=LENIENT
  java -jar \$PICARDJAR MarkDuplicatesWithMateCigar I=${name}.S.BAM  O=${name}.bam  PG=Picard2.9.2_MarkDuplicates M=${name}.MDmetrics.txt  CREATE_INDEX=false VALIDATION_STRINGENCY=LENIENT

  samtools index ${name}.bam

  sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5,5 -k6,6 -m *${type}.bed  >${name}.bed
  rm -f ${name}.US.BAM
  rm -f ${name}.S.BAM
    """
  }

process makeFRBWssds {

  publishDir "${params.outdir}/bigwig", mode: 'copy', overwrite: true, pattern: '*bigwig*'

  cpus 2
  memory 6

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
    java -jar \$PICARDJAR SortSam I=${bam} O=qSort.bam SO=queryname VALIDATION_STRINGENCY=LENIENT
    bedtools bamtobed -mate1 -i qSort.bam -bedpe >frags.bedpe

    perl -lane 'print join("\\t", \$F[0], \$F[1], \$F[5], @F[6..8]) if (\$F[8] eq "+")' frags.bedpe |sort -k1,1 -k2n,2n >frags.POS.bed
    perl -lane 'print join("\\t", \$F[0], \$F[4], \$F[2], @F[6..8]) if (\$F[8] eq "-")' frags.bedpe |sort -k1,1 -k2n,2n >frags.NEG.bed

    mapBed -a win.bed -b frags.POS.bed -c 1 -o count |perl -lane 'use Math::Round; \$p=round((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$p-49,\$p+50,\$F[3])' >pos.bg
    mapBed -a win.bed -b frags.NEG.bed -c 1 -o count |perl -lane 'use Math::Round; \$p=round((\$F[1]+\$F[2])/2); print join("\\t",\$F[0],\$p-49,\$p+50,\$F[3])' >neg.bg

    paste pos.bg neg.bg |perl -lane '\$F[3]+=0.5; \$F[7]+=0.5; \$v=(\$F[3]/\$F[7]); print join("\\t",@F[0..2],(log(\$v)/log(2)))' >fr.bg
    paste pos.bg neg.bg |perl -lane '\$v=(\$F[3]+\$F[7]); print join("\\t",@F[0..2],\$v)' >tot.bg

    bedGraphToBigWig pos.bg \$idx ${iName}.F.bigwig
    bedGraphToBigWig neg.bg \$idx ${iName}.F.bigwig
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

  cpus 1
  memory '4g'

  time '1h'

  input:
  tuple path(bam), path(bai)
  tuple path(ssbam), path(ssbai), path(ssbed)

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

  cpus 1
  memory '4g'

  time '1h'

  input:
  path(reports)

  output:
  path('*ultiQC*', emit: mqcReport)

  script:
  """
  multiqc -f -m ssds -n ${params.name}.multiQC .
  """
  }
