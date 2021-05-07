process minimap2SR {
  label 'aligner'

  time { fq.size() < 1000000000 ? 0.5.hour : 1.hour * (1+fq.size()/1000000000 * task.attempt) }

  input:
  file(fq)

  output:
  tuple path('*.bam'), path('*.bai'), emit: bam

  script:
  def bam   = fq.name.replaceFirst("(.R1.fastq|.R1.fastq.gz|.fastq.gz|.fastq)",".bam")
  """
  minimap2 -ax sr \
    -t ${task.cpus} \
    ${params.genome_mm2idx} \
    ${fq} >tmp.sam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam \
            I=tmp.sam \
            O=${bam} \
            SO=coordinate \
						TMP_DIR="\$TMPDIR" \
            VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}
  """
  }

process minimap2PE {
  label 'aligner'

  time { fq1.size() < 3.GB ? 2.hour : 2.hour + 1.hour * (fq1.size()/3000000000) * task.attempt }

  input:
  tuple path(fq1), path(fq2)

  output:
  tuple path('*.bam'), path('*.bai'), emit: bam

  script:
  def bam   = fq1.name.replaceFirst("(.R1.fastq|.R1.fastq.gz|.fastq.gz|.fastq)",".bam")

  """
  minimap2 -ax sr \
    -t ${task.cpus} \
    ${params.genome_mm2idx} \
    ${fq1} ${fq2} >tmp.sam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" \
            I=tmp.sam \
            O=${bam} \
            SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}
  """
  }

process bwaAlnSR {
  label 'aligner'

  time { fq.size() < 1.GB ? 2.hour : 2.hour + 3.hour * (fq.size()/3000000000) * task.attempt }

  tag { fq }

  input:
  path(fq)

  output:
  tuple path('*.bam'), path('*.bai'), emit: bam

  script:
  def bam   = fq.name.replaceFirst("(.R1.fastq|.R1.fastq.gz|.fastq.gz|.fastq)",".bam")
  """
  bwa aln \
    -t ${task.cpus} \
    ${params.genome_bwaidx} \
    ${fq} >r1.sai

  bwa samse \
    ${params.genome_bwaidx} \
    r1.sai ${fq} >tmp.sam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam  TMP_DIR="\$TMPDIR" \
            I=tmp.sam \
            O=${bam} \
            SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}
  """
  }

process bwaAlnPE {
  label 'aligner'

  time { fq1.size() < 1.GB ? 2.hour : 2.hour + 3.hour * (fq1.size()/3000000000) * task.attempt }

  tag { fq1 }

  input:
  tuple path(fq1), path(fq2)

  output:
  tuple path('*.bam'), path('*.bai'), emit: bam

  script:
  def bam   = fq1.name.replaceFirst("(.R1.fastq|.R1.fastq.gz|.fastq.gz|.fastq)",".bam")
  """
  bwa aln \
    -t ${task.cpus} \
    ${params.genome_bwaidx} \
    ${fq1} >r1.sai

  bwa aln \
    -t ${task.cpus} \
    ${params.genome_bwaidx} \
    ${fq2} >r2.sai

  bwa sampe \
    ${params.genome_bwaidx} \
    r1.sai r2.sai \
    ${fq1} ${fq2} >tmp.sam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" \
            I=tmp.sam \
            O=${bam} \
            SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}
  """
  }

process bwaMemSR {
  label 'aligner'

  time { fq.size() < 3.GB ? 2.hour : 2.hour + 3.hour * (fq.size()/3000000000) * task.attempt }

  tag { fq }

  input:
  path(fq)

  output:
  tuple path('*.bam'), path('*.bai'), emit: bam

  script:
  def bam   = fq.name.replaceFirst("(.R1.fastq|.R1.fastq.gz|.fastq.gz|.fastq)",".bam")
  """
  bwa mem \
    -t ${task.cpus} \
    ${params.genome_bwaidx} ${fq} >tmp.sam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" \
            I=tmp.sam \
            O=${bam} \
            SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}
  """
  }

process bwaMemPE {
  label 'aligner'

  time { fq1.size() < 3.GB ? 2.hour : 2.hour + 3.hour * (fq1.size()/3000000000) * task.attempt }

  tag { fq1 }

  input:
  tuple path(fq1), path(fq2)

  output:
  tuple path('*.bam'), path('*.bai'), emit: bam

  script:
  def bam   = fq1.name.replaceFirst("(.R1.fastq|.R1.fastq.gz|.fastq.gz|.fastq)",".bam")

  """
  bwa mem \
    -t ${task.cpus} \
    ${params.genome_bwaidx} \
    ${fq1} ${fq2} >tmp.sam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" \
            I=tmp.sam \
            O=${bam} \
            SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}
  """
  }

process bowtie2SR {
  label 'aligner'

  time { fq.size() < 3.GB ? 2.hour : 2.hour + 1.hour * (fq.size()/3000000000) * task.attempt }

  tag { fq }

  input:
  file(fq)

  output:
  tuple path('*.bam'), path('*.bai'), emit: bam

  script:
  def bam   = fq.name.replaceFirst("(.R1.fastq|.R1.fastq.gz|.fastq.gz|.fastq)",".bam")
  """
  bowtie2 -x ${params.genome_bt2idx} -1 ${fq} -b init.tmpbam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" \
            I=init.tmpbam \
            O=${bam} \
            SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}
  """
  }

process bowtie2PE {
  label 'aligner'

  time { fq1.size() < 3.GB ? 2.hour : 2.hour + 1.hour * (fq1.size()/3000000000) * task.attempt }

  tag { fq1 }

  input:
  tuple path(fq1), path(fq2)

  output:
  tuple path('*.bam'), path('*.bai'), emit: bam

  script:
  def bam   = fq1.name.replaceFirst("(.R1.fastq|.R1.fastq.gz|.fastq.gz|.fastq)",".bam")
  """
  bowtie2 -x ${params.genome_bt2idx} -1 ${fq1} -2 ${fq2} -b init.tmpbam

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" \
            I=init.tmpbam \
            O=${bam} \
            SO=coordinate \
            VALIDATION_STRINGENCY=LENIENT

  samtools index ${bam}
  """
  }

process mergeBAM {
  label 'mergeBAM'

  publishDir "${params.outdir}/bam",     mode: 'copy', overwrite: true, pattern: '*bam*'
  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: '*txt'

  time { bams[0].size() < 200000000 ? 6.hour : bams[0].size()/200000000 * task.attempt * 6.hour }

  tag { bams.size() }
  input:
  path(bams)

  output:
  tuple(path('*.bam'),path('*.bai'), emit: bam)
  path('*MDmetrics.txt', emit: mdreport)

  script:
  // get INPUT files as string
  def input_args = bams.findAll{ it =~ ".bam\$" }.collect{ "I=$it" }.join(" ")
  def name = "${params.name}.${params.genome}"
  """
  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MergeSamFiles TMP_DIR="\$TMPDIR" \
                 ${input_args} \
                 O=merged.tmpbam \
                 AS=false \
                 SO=coordinate \
                 VALIDATION_STRINGENCY=LENIENT

  if [[ `samtools view -h merged.bam |head -n 100000 |samtools view -f 2 ` ]]; then
	  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MarkDuplicatesWithMateCigar TMP_DIR="\$TMPDIR" \
                   I=merged.tmpbam \
                   O=${name}.bam \
                   PG=Picard2.9.2_MarkDuplicatesWithMateCigar \
                   M=${name}.MDmetrics.txt \
                   MINIMUM_DISTANCE=400 \
			     CREATE_INDEX=false \
			     ASSUME_SORT_ORDER=coordinate \
           VALIDATION_STRINGENCY=LENIENT
  else
    java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MarkDuplicates TMP_DIR="\$TMPDIR" \
                   I=merged.tmpbam \
                   O=${name}.bam \
                   PG=Picard2.9.2_MarkDuplicates \
                   M=${name}.MDmetrics.txt \
			     CREATE_INDEX=false \
			     ASSUME_SORT_ORDER=coordinate \
           VALIDATION_STRINGENCY=LENIENT
  fi

  samtools index ${name}.bam
  """
  }

process mergeBAMv2 {

  label 'mergeBAM'

  time { 5.hour * task.attempt }

  tag { "${name} : ${bams.size()}" }

  input:
  tuple(val(name), path(bams))

  output:
  tuple(val(name), path('*.bam'), path('*.bai'), emit: bam)

  script:
  // get INPUT files as string
  //def input_args = bams.findAll{ it =~ ".bam\$" }.collect{ "I=$it" }.join(" ")
  """
  ## STUPID PAIRTOOLS CAN ADD DUPLICATE PG FIELDS TO BAM
  ## THIS FIXES THE ISSUE (WHICH KILLS PICARD)
  for bam in *.bam; do
    newbam=\${bam/bam/okheader.tmpbam}
    samtools view -h \$bam |perl -lane 'if (\$_ =~/ID:(\\S+)/){\$prog=\$1; if(\$PG{\$prog}++){\$_ =~ s/ID:(\\S+)/ID:\$1\\.\$PG{\$prog}/}}; if (\$_ =~/PN:(\\S+)/){\$prog=\$1; if(\$PN{\$prog}++){\$_ =~ s/PN:(\\S+)/PN:\$1\\.\$PN{\$prog}/}};print \$_;' |samtools view -Shb - >\$newbam
    #samtools view -h \$bam |perl -lane 'if (\$_ =~/ID:(\\S+)/){\$prog=\$1; if(\$PG{\$prog}++){\$_ =~ s/ID:(\\S+)/ID:\$1\\.\$PG{\$prog}/}}; if (\$_ =~/PN:(\\S+)/){\$prog=\$1; if(\$PN{\$prog}++){\$_ =~ s/PN:(\\S+)/PN:\$1\\.\$PN{\$prog}/}}; print \$_;' |samtools view -Shb - >S3ok.bam
  done

  input_args=`ls *okheader.tmpbam |perl -pi -e 's/\\n/ /g' |perl -pi -e 's/(\\S+)/I=\$1/g'`

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MergeSamFiles TMP_DIR="\$TMPDIR" \
                 \$input_args \
                 O=merged.tmpbam \
                 AS=false \
                 VALIDATION_STRINGENCY=LENIENT

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam TMP_DIR="\$TMPDIR" \
                I=merged.tmpbam \
                O=${name}.bam \
                SO=coordinate \
                VALIDATION_STRINGENCY=LENIENT

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR BuildBamIndex \
                I=${name}.bam VALIDATION_STRINGENCY=LENIENT
  """
  }

process markBAMduplicates {

  label 'mergeBAM'

  publishDir "${params.outdir}/bam",     mode: 'copy', overwrite: true, pattern: '*bam*'
  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: '*txt'

  time { bam.size() < 2000000000 ? 6.hour : bam.size()/2000000000 * task.attempt * 6.hour }

  tag { bam }

  input:
  tuple(val(name), path(bam), path(bai))

  output:
  tuple(val(name), path('*.bam'), path('*.bai'), emit: bam)
  path('*MDmetrics.txt', emit: mdreport)

  script:
  """
  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MergeSamFiles TMP_DIR="\$TMPDIR" \
                 I=${bam} \
                 O=merged.tmpbam \
                 AS=true \
                 SO=coordinate \
                 VALIDATION_STRINGENCY=LENIENT

  if [[ `samtools view -h merged.bam |head -n 100000 |samtools view -f 2 ` ]]; then
	  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MarkDuplicatesWithMateCigar TMP_DIR="\$TMPDIR" \
                   I=merged.tmpbam \
                   O=${name}.MD.bam \
                   PG=Picard_MarkDuplicatesWithMateCigar \
                   M=${name}.MDmetrics.txt \
                   MINIMUM_DISTANCE=400 \
			     CREATE_INDEX=false \
			     ASSUME_SORT_ORDER=coordinate \
           VALIDATION_STRINGENCY=LENIENT
  else
    java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR MarkDuplicates TMP_DIR="\$TMPDIR" \
                   I=merged.tmpbam \
                   O=${name}.MD.bam \
                   PG=Picard_MarkDuplicates \
                   M=${name}.MDmetrics.txt \
			     CREATE_INDEX=false \
			     ASSUME_SORT_ORDER=coordinate \
           VALIDATION_STRINGENCY=LENIENT
  fi

  java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR BuildBamIndex \
                I=${name}.MD.bam VALIDATION_STRINGENCY=LENIENT
  """
  }

process getPicardMetrics {

  publishDir "${params.outdir}/reports",  mode: 'copy', overwrite: true

  time { bam.size()< 200000000 ? 2.hour * task.attempt: 2.hour + 2.hour * bam.size()/200000000 * task.attempt }

  tag {bam}

  input:
  tuple path(bam), path(bai)
  val(srpe)

  output:
  path('*Metrics*', emit: report)

  script:
  def name=bam[0].name.replaceFirst(".bam","")
  def picardMem=task.memory.toGiga() - 2

  if (srpe == 'PE'){
	  """
    picard -Xmx${picardMem}g CreateSequenceDictionary \
	      R=${params.genome_fasta} \
	      O=genome.dict \
	      TMP_DIR=\$TMPDIR \
	      VALIDATION_STRINGENCY=LENIENT

    picard -Xmx${picardMem}g CollectAlignmentSummaryMetrics \
	    VALIDATION_STRINGENCY=LENIENT \
	    REFERENCE_SEQUENCE=${params.genome_fasta} \
	    I=${bam} \
	    O=${name}.AlignmentSummary.picardMetrics.tab \
	    TMP_DIR=\$TMPDIR

    picard -Xmx${picardMem}g CollectInsertSizeMetrics \
	      VALIDATION_STRINGENCY=LENIENT \
	      I=${bam} \
	      O=${name}.InsertSize.picardMetrics.tab \
	      H=${name}.InsertSize.picardMetrics.pdf \
	      M=0.5 \
	      TMP_DIR=\$TMPDIR

    picard -Xmx${picardMem}g MeanQualityByCycle \
	      VALIDATION_STRINGENCY=LENIENT \
	      I=${bam} \
	      O=${name}.QByCycle.picardMetrics.tab \
	      TMP_DIR=\$TMPDIR \
	      CHART=${name}.QByCycle.picardMetrics.pdf

    picard -Xmx${picardMem}g QualityScoreDistribution \
	      VALIDATION_STRINGENCY=LENIENT \
	      I=${bam} \
	      O=${name}.QScoreDist.picardMetrics.tab \
	      TMP_DIR=\$TMPDIR \
	      CHART=${name}.QScoreDist.picardMetrics.pdf
	  """
  }else{
	  """
    picard -Xmx${picardMem}g CreateSequenceDictionary \
	      R=${params.genome_fasta} \
	      O=genome.dict \
	      TMP_DIR=\$TMPDIR \
	      VALIDATION_STRINGENCY=LENIENT

    picard -Xmx${picardMem}g MeanQualityByCycle \
	      VALIDATION_STRINGENCY=LENIENT \
	      I=${bam} \
	      O=${name}.QByCycle.picardMetrics.tab \
	      TMP_DIR=\$TMPDIR \
	      CHART=${name}.QByCycle.picardMetrics.pdf

    picard -Xmx${picardMem}g QualityScoreDistribution \
	      VALIDATION_STRINGENCY=LENIENT \
	      I=${bam} \
	      O=${name}.QScoreDist.picardMetrics.tab \
	      TMP_DIR=\$TMPDIR \
	      CHART=${name}.QScoreDist.picardMetrics.pdf
	  """
  }
  }

process bamToBW {

  publishDir "${params.outdir}/bigwig",  mode: 'copy', overwrite: true, pattern: '*bigwig'

  time { bam.size()< 1000000000 ? 0.5.hour * task.attempt: 1.hour * bam.size()/1000000000 * task.attempt }

	tag { bam }

  input:
  tuple path(bam), path(bai)

  output:
  path('*bigwig', emit: bw)

  script:
  def name=bam.name.replaceFirst(".bam","")
  """
  if [[ `samtools view -F 4 ${bam} |head -n 5001 |wc -l` -gt 5000 ]]; then
    genomeCoverageBed -ibam ${bam} -bg >${name}.tmpbg
    ## SORT BY INDEX TO ASSURE WE ONLY HAVE CHROMOSOMES MENTIONED IN THE INDEX THEN DO A LINUX SORT
    ## THIS IS A BIT ASSWAYS; FIX LATER
    sortBEDByFAI.pl ${name}.tmpbg ${params.genome_fai} |sort -k1,1 -k2n,2n -k3n,3n >${name}.bedgraph
    bedGraphToBigWig ${name}.bedgraph ${params.genome_fai} ${name}.bedtools.bigwig
  else
    touch EMPTY_TOOFEWREADS_${name}.bedtools.bigwig
  fi
  """
  }

process trimFASTQsr {

  label 'trimFQ'

	time { fq.size() < 3.GB ? 1.hour : 1.hour + 1.hour * (fq.size()/3000000000) * task.attempt }

  tag {fq}

  input:
  path(fq)

  output:
  path "${params.name}.trimmed.R1.fastq", emit: fq
  //path "*_report.txt"                   , emit: report

  script:
  """
  if [[ "${fq}" =~ .gz\$ ]]; then gz=${fq}; fq=\${fq/.gz/}; zcat ${fq} >\$fq; fi
  fqr1=`ls *.R1.fastq`

  #trimFQ1=`echo "\$fqr1" |perl -pi -e 's/.fastq/_val_1.fq/'`
  trimFQ1=\${fqr1/.fastq/_trimmed.fq}

  trim_galore -q 10 --dont_gzip --stringency 6 --length 25 \$fqr1

  mv \$trimFQ1 ${params.name}.trimmed.R1.fastq

  mv \$fqr1"_trimming_report.txt" ${params.name}.R1_trimgalore_trimming_report.txt

  """
  }

process trimFASTQpe {

  label 'trimFQ'

  time { fq[0].size() < 3.GB ? 1.hour : 1.hour + 1.hour * (fq[0].size()/3000000000) * task.attempt }

  tag { fq[0] }

  input:
  path(fq)

  output:
  tuple path("${params.name}.trimmed.R1.fastq"), path("${params.name}.trimmed.R2.fastq"), emit: fq
  //path "*_report.txt"                                                                   , emit: report

  script:
  def fqR1=fq[0]
  def fqR2=fq[1]
  """
  if [[ "${fqR1}" =~ .gz\$ ]]; then gz=${fqR1}; fq=\${gz/.gz/}; zcat ${fqR1} >\$fq; fi
  if [[ "${fqR2}" =~ .gz\$ ]]; then gz=${fqR2}; fq=\${gz/.gz/}; zcat ${fqR2} >\$fq; fi

  fqr1=`ls *.R1.fastq`
  fqr2=`ls *.R2.fastq`
  trimFQ1=`echo "\$fqr1" |perl -pi -e 's/.fastq/_val_1.fq/'`
  trimFQ2=`echo "\$fqr2" |perl -pi -e 's/.fastq/_val_2.fq/'`

  trim_galore -q 10 --paired --dont_gzip --stringency 6 --length 25 \$fqr1 \$fqr2

  mv \$trimFQ1 ${params.name}.trimmed.R1.fastq
  mv \$trimFQ2 ${params.name}.trimmed.R2.fastq

  mv \$fqr1"_trimming_report.txt" ${params.name}.R1_trimgalore_trimming_report.txt
  mv \$fqr2"_trimming_report.txt" ${params.name}.R2_trimgalore_trimming_report.txt

  """
  }

process makeDeeptoolsBW {

  publishDir "${params.outdir}/bigwig",  mode: 'copy', overwrite: true, pattern: '*bigwig'
  publishDir "${params.outdir}/plots",   mode: 'copy', overwrite: true, pattern: '*png'
  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: '*tab'

  time { bam.size()< 500000000 ? 1.hour : 1.hour + 1.hour * bam.size()/500000000 * task.attempt }

  tag {bam}

  input:
  tuple(path(bam), path(bai))

  output:
  path('*bigwig', emit: bw)
  path('*png',    emit: png)
  path('*tab',    emit: tab)

  script:
  def name            = bam.name.replaceAll(/.bam/,'')
  """
  if [[ `samtools view -F 4 ${bam} |head -n 5001 |wc -l` -gt 5000 ]]; then

   ## Added KB 2003-03-19:
   ## Split BAM by chromosome
   for cs in `samtools idxstats ${bam} |cut -f1 |grep -vP '\\*'`; do
     csbam=\$cs".bam"
     csBW=\$cs".ALL.BW"
     csBWND=\$cs".ND.BW"
     cswig=\$cs".ALL.wig"
     cswigND=\$cs".ND.wig"

     if [ `samtools view -F 4 ${bam} \$cs |head -n 10 |wc -l` -eq 10 ]; then
       samtools view -F 4 -hb ${bam} \$cs >\$csbam
       samtools index \$csbam

       bamCoverage --bam \$csbam --binSize 150 --normalizeUsing RPKM \
         -p max -v -o \$csBW
       bigWigToWig \$csBW \$cswig
       grep -w \$cs \$cswig >>all.wig

       bamCoverage --bam \$csbam --binSize 150 --normalizeUsing RPKM \
           --ignoreDuplicates -p max -v -o \$csBWND
       bigWigToWig \$csBWND \$cswigND
       grep -w \$cs \$cswigND >>ND.wig
     fi
   done

   wigToBigWig all.wig ${params.genome_fai} ${name}.deeptools.150bp.RPKM.bigwig
   wigToBigWig ND.wig  ${params.genome_fai} ${name}.deeptools.150bp.RPKM.noDups.bigwig

   plotCoverage --bamfiles ${bam} --numberOfProcessors ${task.cpus} \
     -o ${name}.deeptools.coveragePlot.png

   plotFingerprint --bamfiles ${bam} --labels ${bam} --numberOfProcessors ${task.cpus} \
     --minMappingQuality 30 --skipZeros \
     --plotFile ${name}.deeptools.fingerprints.png \
     --outRawCounts ${name}.deeptools.fingerprints.tab

  else
    touch EMPTY_TOOFEWREADS_${name}.bigwig
    touch EMPTY_TOOFEWREADS_${name}.deeptools.coveragePlot.png
    touch EMPTY_TOOFEWREADS_${name}.deeptools.fingerprints.png
    touch EMPTY_TOOFEWREADS_${name}.deeptools.fingerprints.tab
  fi
  """
  }

process samStats {

  publishDir "${params.outdir}/reports",  mode: 'copy', overwrite: true

  time { bam.size()< 10000000000 ? 1.hour : 1.hour * bam.size()/10000000000 * task.attempt }

  tag {bam}

  input:
  tuple path(bam),  path(idx)

  output:
  path('*stats.tab', emit: report)

  script:
  iStat = bam.name.replaceAll(/.bam/,".idxstats.tab")
  sStat = bam.name.replaceAll(/.bam/,".samstats.tab")
  """
  samtools idxstats ${bam} >${iStat}
  samtools stats ${bam} >${sStat}
  """
  }


  process makeFRBW {
    publishDir "${params.outdir}/bigwig",  mode: 'copy', overwrite: true, pattern: '*bigwig'

    cpus 2
    memory 6

    time { bam.size()< 5000000000 ? 1.hour : 1.hour * bam.size()/5000000000 * task.attempt }

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
      java -jar -Xmx${task.memory.toGiga()}g \$PICARDJAR SortSam I=${bam} O=qSort.bam SO=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR="\$TMPDIR"
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

// OK ... let's start
workflow alignFromFQ {
  take:
  fastq

  main:
  if (params.isPE) {
    switch (params.aligner) {
      case 'bwa':
        trimFASTQpe(fastq) | splitFastq (by: params.splitSz, pe: true, file:true) | bwaMemPE |collect |mergeBAM
        break
      case 'bwaaln':
        trimFASTQpe(fastq) | splitFastq (by: params.splitSz, pe: true, file:true) | bwaAlnPE |collect |mergeBAM
        break
      case 'mm2':
        trimFASTQpe(fastq) | splitFastq (by: params.splitSz, pe: true, file:true) | minimap2PE |collect |mergeBAM
        break
      case 'bowtie' :
        trimFASTQpe(fastq) | splitFastq (by: params.splitSz, pe: true, file:true) | bowtie2PE |collect |mergeBAM
        break
    }
  }else{
    switch (params.aligner) {
      case 'bwa':
        trimFASTQsr(fastq) | splitFastq (by: params.splitSz, file:true) | bwaMemSR |collect |mergeBAM
        break
      case 'bwaaln':
        trimFASTQsr(fastq) | splitFastq (by: params.splitSz, file:true) | bwaAlnSR |collect |mergeBAM
        break
      case 'mm2':
        trimFASTQsr(fastq) | splitFastq (by: params.splitSz, file:true) | minimap2SR |collect |mergeBAM
        break
      case 'bowtie' :
        trimFASTQsr(fastq) | splitFastq (by: params.splitSz, file:true) | bowtie2SR |collect |mergeBAM
        break
    }
  }

  getPicardMetrics(mergeBAM.out.bam,params.isPE?"PE":"SR")
  bamToBW(mergeBAM.out.bam)
  makeDeeptoolsBW(mergeBAM.out.bam)
  samStats(mergeBAM.out.bam)

  emit:
  bam    = mergeBAM.out.bam
  repMD  = mergeBAM.out.mdreport
  repAln = getPicardMetrics.out.report
  repDT  = makeDeeptoolsBW.out.tab
  repST  = samStats.out.report
  bwBT   = bamToBW.out.bw
  bwDT   = makeDeeptoolsBW.out.bw

  }
