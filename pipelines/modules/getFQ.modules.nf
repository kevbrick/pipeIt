nextflow.preview.dsl=2

process SRAtoFQ {
  cpus 1
  memory '4 GB'

  time { 4.hour * task.attempt }
  errorStrategy { 'retry' }
  maxRetries 3

  input:
  val sraNum
  output:
  path('*.R[12].fastq', emit: fqTuple)

  script:
  """
  sra=\$(echo "${sraNum}" | tr "," "\\n")

  for s in \$sra
  do
    fastq-dump \$s --split-files >/dev/null 2>/dev/null

    if [ -f \$s"_2.fastq" ] ; then
      mv \$s"_1.fastq" \$s".R1.fastq"
      mv \$s"_2.fastq" \$s".R2.fastq"
    else
      mv \$s"_1.fastq" \$s".R1.fastq"
    fi
  done
  """
  }

process BAMtoFQ {
  cpus 4
  memory '12 GB'

  time { bam.size() < 3.GB ? 1.hour : 1.hour + 1.hour * (bam.size()/3000000000) * task.attempt }

  errorStrategy { 'retry' }
  maxRetries 2

  tag { bam.size() }

  input:
  path(bam)

  output:
  path('*.R[12].fastq', emit: fqTuple)

  script:
  """
  n=`samtools view -h ${bam} |head -n 100000 |samtools view -f 1 -S /dev/stdin |wc -l`

  if [ \$n -eq 0 ]; then
    isSRPE="SR"

    java -jar -Xmx8g \$PICARDJAR SortSam \
                   I=${bam}  \
                   O=querynameSort.bam \
                   SO=queryname \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null

    java -jar -Xmx8g \$PICARDJAR SamToFastq I=querynameSort.bam \
                   F=nxf.R1.fastq \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null
  else
    isSRPE="PE"
    java -jar -Xmx8g \$PICARDJAR FixMateInformation \
                   I=${bam} \
                   O=fixMate.bam \
                   SORT_ORDER=queryname \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null

    java -jar -Xmx8g \$PICARDJAR SortSam \
                   I=fixMate.bam  \
                   O=querynameSort.bam \
                   SO=queryname \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null

    java -jar -Xmx8g \$PICARDJAR SamToFastq I=querynameSort.bam \
                   F=nxf.R1.fastq F2=nxf.R2.fastq \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null
    fi
  """
  }

process FQtoFQpe {
  cpus 2
  memory '4 GB'

  time { inFQ1.size() < 3.GB ? 1.hour : 1.hour + 1.hour * (inFQ1.size()/3000000000) * task.attempt }

  errorStrategy { 'retry' }
  maxRetries 1

  tag { inFQ1 }

  input:
  path(inFQ1)
  path(inFQ2)

  output:
  path('*.R[12].fastq', emit: fqs)

  script:
  if (inFQ1 =~ /.gz$/)
    """
    gunzip --stdout ${inFQ1} >nxf.R1.fastq 2>/dev/null
    gunzip --stdout ${inFQ2} >nxf.R2.fastq 2>/dev/null
    isSRPE="PE"
    """
  else
    """
    ln -s ${inFQ1} nxf.R1.fastq 2>/dev/null
    ln -s ${inFQ2} nxf.R2.fastq 2>/dev/null
    isSRPE="PE"
    """
  }

process FQtoFQsr {
  cpus 2
  memory '4 GB'

  time { inFQ1.size() < 3.GB ? 1.hour : 1.hour + 1.hour * (inFQ1.size()/3000000000) * task.attempt }

  errorStrategy { 'retry' }
  maxRetries 1

  tag { inFQ1 }

  input:
  path(inFQ1)

  output:
  path('*.R[12].fastq', emit: fqs)

  script:
  if (inFQ1 =~ /.gz$/)
    """
    gunzip --stdout ${inFQ1} >nxf.R1.fastq 2>/dev/null
    isSRPE="SR"
    """
  else
    """
    ln -s ${inFQ1} nxf.R1.fastq
    isSRPE="PE"
    """
  }

process OBJtoFQ{
  cpus 1
  memory '4 GB'

  time { 3.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  tag { objName }

  input:
  val(objName)

  output:
  path('*.R[12].fastq', emit: fqTuple)

  script:
  """
  getFromBiowulfObj --v RDCO --p ${objName} --get --f
  gunzip *fastq.gz
  """
  }

process mergeFQ {

  cpus 4
  memory '6 GB'

  time { 6.hour }

  errorStrategy { 'retry' }
  maxRetries 1

  tag {params.name}

  input:
  path(fqs)

  output:
  path("*.R[12].fastq*", emit: fqs)

  script:
  // Modified to work better with huge files
  // Now, we make soft links if at all possible instead of copying stuff
  // This also speeds things up quite a bit
  def name="${params.name}.merged"
  """
  R1files=\$(ls *R1*fastq 2> /dev/null | wc -l)
  R2files=\$(ls *R2*fastq 2> /dev/null | wc -l)

  if [ "\$R2files" != "0" ]
  then
    if [ "\$R2files" != "1" ]
    then
      cat *R2*fastq >mergeR2.fastq
      rm *.R2.fastq
    else
      fqR2=`ls *R2*fastq`
      ln -s \$fqR2 mergeR2.fastq
    fi

    if [ "${params.sortFQ}" == "true" ]; then
      fastq-sort --id mergeR2.fastq >${name}.R2.fastq
      rm mergeR2.fastq
    else
      mv mergeR2.fastq ${name}.R2.fastq
    fi
  fi

  if [ "\$R1files" != "1" ]
  then
    cat *.R1.fastq >mergeR1.fastq
    rm  *.R1.fastq
  else
    fqR1=`ls *R1*fastq`
    ln -s \$fqR1 mergeR1.fastq
  fi

  if [ "${params.sortFQ}" == "true" ]; then
    fastq-sort --id mergeR1.fastq >${name}.R1.fastq
    rm mergeR1.fastq
  else
    mv mergeR1.fastq ${name}.R1.fastq
  fi

  FQsize=`du -k "${name}.R1.fastq" | cut -f1`

  if [ "${params.gzipoutput}" == "true" ]
  then
      gzip ${name}.R1.fastq -c >${name}.R1.fastq.gz
      gzip ${name}.R2.fastq -c >${name}.R2.fastq.gz|| true

      rm -f ${name}.R1.fastq ${name}.R2.fastq
  fi
  """
  }

process commitToObj {

  cpus 1
  memory '6 GB'

  time { 1.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  input:
  path(fqs)

  script:
  """
  for f in `ls *fastq.gz`; do
    #obj_rm  -v RDCO \$f 2>/dev/null || true
    obj_put --force -v RDCO \$f
  done
  """
  }

process fastqC {

  publishDir "${params.outdir}/reports",  mode: 'copy', overwrite: true

  cpus 1
  memory '4 GB'

  time { 1.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  tag {fq}

  input:
  path(fq)

  output:
  path('*fastqc*', emit: fqc)

  script:
  """
  for gz in *fastq.gz; do
    fq=\${gz/.gz/}
    gunzip -c \$gz |head -n 100000000 |tail -n 1000000 >\$fq
  done

  for f in *fastq; do
    head -n 10000000 \$f >subset.fastq

    fastqc -t 4 subset.fastq

    mv subset_fastqc.html \$f"c_report.html"
    mv subset_fastqc.zip  \$f"c_report.zip"
  done
  """
  }

process fastqScreen {

  publishDir "${params.outdir}/reports",  mode: 'copy', overwrite: true

  cpus 4
  memory '12 GB'

  time { 1.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  tag {fq}

  input:
  path(fq)

  output:
  path('*_screen*', emit: report)

  script:
  """
  ##Make FastqScreen config file
  echo "# This is a config file for FastQ Screen\n\n" >conf.txt
  echo "THREADS ${task.cpus}\n\n" >>conf.txt

  gList="${params.genomes2screen}"
  g2use=\${gList//,/\$'\n'}

  for g in \$g2use; do
    echo "DATABASE \$g \$NXF_GENOMES/\$g/BWAIndex/version0.7.10/genome.fa" >>conf.txt
  done

  for f in *.R1.fastq*; do
    fastq_screen \
      --threads ${task.cpus} --force \
      --aligner bwa \$f \
      --conf conf.txt
  done
  """
  }

// OK ... let's start
workflow getFQs {

  switch (params.inputType) {
    case 'sra':
      SRAtoFQ(params.sra) |mergeFQ |flatten | collect | (fastqC & fastqScreen)
      break
    case 'bam':
      BAMtoFQ(file(params.bam)) |mergeFQ | (fastqC & fastqScreen)
      break
    case 'obj':
      OBJtoFQ(params.obj) |mergeFQ | (fastqC & fastqScreen)
      break
    case 'fqsr' :
      FQtoFQsr(file(params.fq1)) |mergeFQ | (fastqC & fastqScreen)
      break
    case 'fqpe' :
      FQtoFQpe(file(params.fq1),file(params.fq2)) |mergeFQ |flatten | collect | (fastqC & fastqScreen)
      break
  }

  emit:
  fq    = mergeFQ.out.fqs
  fqc   = fastqC.out.fqc
  fqscr = fastqScreen.out.report
  }
