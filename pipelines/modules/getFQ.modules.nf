nextflow.preview.dsl=2

process SRAtoFQ {

  label 'getFQ'

  time { 4.hour * task.attempt }

  input:
  val sraNum

  output:
  tuple(val("gz"), path('*.R1.fastq.gz'), emit: R1)
  tuple(val("gz"), path('*.R2.fastq.gz'), emit: R2, optional: true)

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

  gzip *R?.fastq
  """
  }

process BAMtoFQ {

  label 'getFQ'

  time { 1.hour * bam.size()/4000000000 * task.attempt}

  tag { bam.size() }

  input:
  file(bam)

  output:
  tuple(val("gz"), path('*.R1.fastq.gz'), emit: R1)
  tuple(val("gz"), path('*.R2.fastq.gz'), emit: R2, optional: true)

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
                   F=nxf.R1.fastq.gz \
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
                   F=nxf.R1.fastq.gz F2=nxf.R2.fastq.gz \
                   TMP_DIR=\$TMPDIR \
                   VALIDATION_STRINGENCY=LENIENT >/dev/null 2>/dev/null
    fi
  """
  }

process OBJtoFQ{

  label 'getFQ'

  time { 3.hour * task.attempt}

  tag { objName }

  input:
  val(objName)

  output:
  tuple(val("gz"), path('*.R1.fastq.gz'), emit: R1)
  tuple(val("gz"), path('*.R2.fastq.gz'), emit: R2, optional: true)

  script:
  """
  echo "OBJ"
  getFromBiowulfObj --v RDCO --p ${objName} --get --f
  """
  }

process mergeFQ {

  tag {"${params.name}"}

  input:
  tuple(val(type),file(R1),file(R2))

  output:
  path("*.R[12].fastq*", emit: fqs)

  script:
  def multiR1 = R1[1]
  def multiR2 = R2[1]
  def name = "${params.name}.merged"
  // Modified to work better with huge files
  // Now, we make soft links if at all possible instead of copying stuff
  // This also speeds things up quite a bit
  """
  if [ "${type}" == "gz" ]; then
    if [ "${multiR1}" == "null" ]; then
      ln -s ${R1} ${name}.R1.fastq.gz
    else
      cat ${R1} >${name}.R1.fastq.gz
    fi

    #if [ ! -z "${R2}" ]; then
    if [ "${R2}" != "input.2" ]; then
      if [ "${multiR2}" == "null" ]; then
        ln -s ${R2} ${name}.R2.fastq.gz
      else
        cat ${R2} >${name}.R2.fastq.gz
      fi
    fi
  else
    cat ${R1} |gzip -c >${name}.R1.fastq

    #if [ ! -z "${R2}" ]; then
    if [ "${R2}" != "input.2" ]; then
      cat ${R2} |gzip -c >${name}.R2.fastq
    fi
  fi

  if [ "${params.sortFQ}" == "true" ]; then
    if [ "${type}" == "gz" ]; then
      gunzip -c ${name}.R1.fastq.gz >tmp.fq
      fastq-sort --id tmp.fq |gzip -c >${name}.R1.fastq.gz
      rm tmp.fq

      #if [ ! -z "${R2}" ]; then
      if [ "${R2}" != "input.2" ]; then
        gunzip -c ${name}.R2.fastq.gz >tmp.fq
        fastq-sort --id tmp.fq |gzip -c >${name}.R2.fastq.gz
        rm tmp.fq
      fi
    else
      fastq-sort --id ${name}.R1.fastq R1.fastq
      mv R1.fastq ${name}.R1.fastq
      gzip ${name}.R1.fastq

      #if [ ! -z "${R2}" ]; then
      if [ "${R2}" != "input.2" ]; then
        fastq-sort --id ${name}.R2.fastq R2.fastq
        mv R2.fastq ${name}.R2.fastq
        gzip ${name}.R2.fastq
      fi
    fi
  fi

  """
  }

process commitToObj {

  cpus 1
  memory '6g'

  time { 5.hour * task.attempt}
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
  tag {"${params.name}"}

  publishDir "${params.outdir}/reports",  mode: 'copy', overwrite: true

  input:
  path(fq)

  output:
  path('*fastqc_report.*', emit: rep)

  when:
  !("${workflow.scriptName}" =~ /commit/)

  script:
  """
  for f in *fastq*; do
    if [[ "\$f" =~ ".gz" ]]; then
      zcat \$f |head -n 10000000 >subset.fastq
      name=\${f/.(fastq|fq).gz/}
    else
      head -n 10000000 \$f >subset.fastq
      name=\${f/.(fastq|fq).gz/}
    fi

    fastqc -t ${task.cpus} subset.fastq

    mv subset_fastqc.html \$name".fastqc_report.html"
    mv subset_fastqc.zip  \$name".fastqc_report.zip"
  done
  """
  }

process fastqScreen {
  tag {"${params.name}"}

  publishDir "${params.outdir}/reports",  mode: 'copy', overwrite: true

  when:
  !("${workflow.scriptName}" =~ /commit/)

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
    if [[ "\$f" =~ "q.gz" ]]; then
      fqs=`echo \$f |perl -pi -e 's/(fq|fastq).gz/ds.fastq/' 2>/dev/null`
      zcat \$f |head -n 4000000 |tail -n 400000 >\$fqs
    else
      fqs=`echo \$f |perl -pi -e 's/(fq|fastq)\$/ds.fastq/' 2>/dev/null`
      cat \$f |head -n 4000000 |tail -n 400000 >\$fqs
    fi
    
    fastq_screen \
      --threads ${task.cpus} --force \
      --aligner bwa \$fqs \
      --conf conf.txt
  done
  """
  }

// OK ... let's start
workflow getFQs {

  switch (params.inputType) {
    case 'sra':
      init  = SRAtoFQ(params.sra)
      fqs   = init.R1.join(init.R2, by: 0, remainder: true)
      merge = mergeFQ(fqs) | (fastqC & fastqScreen)
      break
    case 'bam':
      init  = BAMtoFQ(params.sra)
      fqs   = init.R1.join(init.R2, by: 0, remainder: true)
      merge = mergeFQ(fqs) | (fastqC & fastqScreen)
      break
    case 'obj':
      init  = OBJtoFQ(params.obj)
      fqs   = init.R1.join(init.R2, by: 0, remainder: true)
      merge = mergeFQ(fqs) | (fastqC & fastqScreen)
      break
    case 'fqsr' :
      fqs = Channel.from(params.fq1 =~ /gz/ ? 'gz' : 'fastq',file(params.fq1)).collect()
      fqs.view()
      merge = mergeFQ(fqs) | (fastqC & fastqScreen)
      break
    case 'fqpe' :
      fqs = Channel.from(params.fq1 =~ /gz/ ? 'gz' : 'fastq',file(params.fq1),file(params.fq2)).collect()
      merge = mergeFQ(fqs) | (fastqC & fastqScreen)
      break
  }

  emit:
  fq    = mergeFQ.out.fqs
  fqc   = fastqC.out.rep
  fqscr = fastqScreen.out.report
  }
