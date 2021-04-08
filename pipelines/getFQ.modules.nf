nextflow.preview.dsl=2

process SRAtoFQ {
  cpus 1
  memory '4g'

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
  memory '4g'

  time { 1.hour * bam.size()/4000000000 * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 2

  tag { bam.size() }

  input:
  file(bam)

  output:
  path('*.R[12].fastq', emit: fqTuple)

  script:
  """
  n=`samtools view -h ${bam} |head -n 100000 |samtools view -f 1 -S /dev/stdin |wc -l`

  if [ \$n -eq 0 ]; then
    isSRPE="SR"

    java -jar \$PICARDJAR SortSam \
                   I=${bam}  \
                   O=querynameSort.bam \
                   SO=queryname \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null

    java -jar \$PICARDJAR SamToFastq I=querynameSort.bam \
                   F=nxf.R1.fastq \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null
  else
    isSRPE="PE"
    java -jar \$PICARDJAR FixMateInformation \
                   I=${bam} \
                   O=fixMate.bam \
                   SORT_ORDER=queryname \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null

    java -jar \$PICARDJAR SortSam \
                   I=fixMate.bam  \
                   O=querynameSort.bam \
                   SO=queryname \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null

    java -jar \$PICARDJAR SamToFastq I=querynameSort.bam \
                   F=nxf.R1.fastq F2=nxf.R2.fastq \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null
    fi
  """
  }

process FQtoFQpe {
  cpus 2
  memory '4g'

  time { 1.hour * inFQ1.size()/5000000000 * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  tag { inFQ1 }

  input:
  file inFQ1
  file inFQ2

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
  memory '4g'

  time { 1.hour * inFQ1.size()/5000000000 * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  tag { inFQ1 }

  input:
  file inFQ1

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
  memory '4g'

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
  memory '6g'

  time { 2.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  input:
  file(fqs)

  output:
  path("*.R[12].fastq*", emit: fqs)

  script:
  // Modified to work better with huge files
  // Now, we make soft links if at all possible instead of copying stuff
  // This also speeds things up quite a bit
  """
  R1files=\$(ls *R1*fastq 2> /dev/null | wc -l)
  R2files=\$(ls *R2*fastq 2> /dev/null | wc -l)

  if [ "\$R2files" != "0" ]
  then
    if [ "\$R2files" != "1" ]
    then
      cat *R2*fastq >mergeR2.fastq
      rm -f *.R2.fastq
    else
      fqR2=`ls *R2*fastq`
      ln -s \$fqR2 mergeR2.fastq
    fi
    fastq-sort --id mergeR2.fastq >${params.name}.R2.fastq
    rm -f mergeR2.fastq
  fi

  if [ "\$R1files" != "1" ]
  then
    cat *.R1.fastq >mergeR1.fastq
    rm  -f *.R1.fastq
  else
    fqR1=`ls *R1*fastq`
    ln -s \$fqR1 mergeR1.fastq
  fi

  fastq-sort --id mergeR1.fastq >${params.name}.R1.fastq
  rm -f mergeR1.fastq

  FQsize=`du -k "${params.name}.R1.fastq" | cut -f1`

  if [ "${params.gzipoutput}" == "true" ]
  then
      gzip ${params.name}.R1.fastq
      gzip ${params.name}.R2.fastq || true
  fi
  """
  }

process commitToObj {

  cpus 1
  memory '6g'

  time { 1.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  input:
  file(fqs)

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
  memory { 4.gb }

  time { 1.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

  input:
  path(fq)

  output:
  path('*fastqc*', emit: fqc)

  script:
  """
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
  memory '4g'

  time { 1.hour * task.attempt}
  errorStrategy { 'retry' }
  maxRetries 1

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

  for f in *.R1.fastq; do
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
