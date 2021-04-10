// KB: March 11 2021
// These modules are based on the 4D nucleome protocols
// Detailed information is at:
// https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline#recheck

process bwa4D {
  tag { fq1 }

  publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true, pattern: '*bam'

  input:
  tuple(val(id), val(sample), val(rep), val(bio), path(fq1), path(fq2))

  output:
  tuple(val(sample), val(rep), val(bio), path('*BWAinit.bam'))

  script:
  def nm = fq1.name.replaceAll("R1.(\\d+).fastq.gz","\$1")
  //println txt.replaceAll("(kevin)_(\\d+)_","\$2--\$1")
  """
  ln -s ${fq1} fastq1.gz
  ln -s ${fq2} fastq2.fz

  bwa mem \
    -t ${task.cpus} \
    -SP5M \
    ${params.bwaidx} \
    ${fq1} ${fq2} | samtools view -Shb - >${nm}.BWAinit.bam
  """
  }

process bowtie2_end2end {
  tag { fq1 }

  input:
  tuple(val(id), val(sample), val(rep), val(bio), path(fq1), path(fq2))
  val(trim)
  val(removeUnaligned)

  output:
  tuple(val(sample), val(rep), val(bio), val(fq1.name), val('R1'), path('*R1*endtoend.bam'), emit: bamR1)
  tuple(val(sample), val(rep), val(bio), val(fq1.name), val('R2'), path('*R2*endtoend.bam'), emit: bamR2)
  tuple(val(sample), val(rep), val(bio), val(fq1.name), val('R1'), path('*R1*.unmap.fastq'), emit: unmappedR1)
  tuple(val(sample), val(rep), val(bio), val(fq1.name), val('R2'), path('*R2*.unmap.fastq'), emit: unmappedR2)

  script:
  def nm = fq1.name.replaceAll("R1.(\\d+).fastq","\$1")
  """
  if [ "${trim}" -gt 0 ]; then
    zcat ${fq1} |perl -lane 'chomp; \$l++; if (\$l == 1 || \$l == 3){\$out = \$_}else{\$out = substr(\$_,1,${trim})}; print \$out; \$l=(\$l==4)?0:\$l' >R1.fastq
    zcat ${fq2} |perl -lane 'chomp; \$l++; if (\$l == 1 || \$l == 3){\$out = \$_}else{\$out = substr(\$_,1,${trim})}; print \$out; \$l=(\$l==4)?0:\$l' >R2.fastq
    name="${nm}_BOWTIE_trim${trim}"
  else
    ln -s ${fq1} R1.fastq.gz
    ln -s ${fq2} R2.fastq.gz
    name="${nm}_BOWTIE${trim}"
  fi

  if [ ${removeUnaligned} -eq 0 ]; then
    bowtie2 --rg-id BMG --rg SM:bowtiealn \\
        	--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder \\
          -p ${task.cpus} \\
          -x ${params.bt2idx} \\
          --un \${name}".bowtie.R1.unmap.fastq" \\
        	-U R1.fastq  | samtools view -bS - > \${name}".R1.bowtie_endtoend.bam"

    bowtie2 --rg-id BMG --rg SM:bowtiealn \\
        	--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder \\
          -p ${task.cpus} \\
          -x ${params.bt2idx} \\
          --un \${name}".bowtie.R2.unmap.fastq" \\
        	-U R2.fastq  | samtools view -bS - > \${name}".R2.bowtie_endtoend.bam"
  else
    bowtie2 --rg-id BMG --rg SM:bowtiealn \\
          --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder \\
          -p ${task.cpus} \\
          -x ${params.bt2idx} \\
          --un \${name}".bowtie.R1.unmap.fastq" \\
          -U R1.fastq  | samtools view -F 4 -bS - > \${name}".R1.bowtie_endtoend.bam"

    bowtie2 --rg-id BMG --rg SM:bowtiealn \\
          --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder \\
          -p ${task.cpus} \\
          -x ${params.bt2idx} \\
          --un \${name}".bowtie.R2.unmap.fastq" \\
          -U R2.fastq  | samtools view -F 4 -bS - > \${name}".R2.bowtie_endtoend.bam"
  fi
  rm R1.fastq R2.fastq

  """
  }

process trim_hic_reads {

  input:
  tuple(val(sample), val(rep), val(bio), val(name), val(read), path(fq))

  output:
  tuple(val(sample), val(rep), val(bio), val(name), val(read), path("*trimmed.fastq"))

  script:
  def outfq=fq.name.replaceFirst('unmap.fastq.*','trimmed.fastq')
  """
  /HiC-Pro_3.0.0/scripts/cutsite_trimming --fastq ${fq} \\
                   --cutsite  ${params.ligation_site} \\
                   --out ${outfq}
  """
  }

process bowtie2_on_trimmed_reads {

  input:
  tuple(val(sample), val(rep), val(bio), val(name), val(read), path(fq))

  output:
  tuple(val(sample), val(rep), val(bio), val(name), val(read), path("*bam"))

  script:
  def bam=fq.name.replaceFirst('trimmed.fastq.*','trimmed.bam')
  """
  bowtie2 --rg-id BMG --rg SM:${bam} \
          --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder \
          -p ${task.cpus} \
          -x ${params.bt2idx} \
          -U ${fq} | samtools view -bS - > ${bam}
  """
  }

process bowtie2_mergeR1R2{

  input:
  tuple(val(sample), val(rep), val(bio), val(name), val(read), path(bams))

  output:
  tuple(val(sample), val(rep), val(bio), val(name), path("*merged.bam"), emit: bam)

  script:
  def nm = name.replaceAll("R1.(\\d+).fastq.*","\$1")
  def outbam="${nm}.${read}.merged.bam"
  """
  samtools merge -@ ${task.cpus} \
               -f tmp.bam \
             ${bams[0]} ${bams[1]}

  samtools sort -@ ${task.cpus} -m 800M \
  	              -n -T \$TMPDIR \
                  -o ${outbam} \
                  tmp.bam
  """
  }

process bowtie2_make_paired_bam{

  tag { fqName.replaceFirst(".fastq.*","") }

  //publishDir "${params.outdir}/bowtie", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), val(rep), val(bio), val(fqName), path(read1), path(read2))

  output:
  tuple(val(sample), val(rep), val(bio), path('*bt2.bam'))

  script:
  def opts = " -t --single --multi -q 10 "
  def name = fqName.replaceFirst(".fastq.*","")
  """
  /HiC-Pro_3.0.0/scripts/mergeSAM.py -f ${read1} -r ${read2} -o ${name}.bt2.bam ${opts}

  #/HiC-Pro_3.0.0/scripts/mergeSAM.py -f \${name}".R1.bowtie_endtoend.bam" -r \${name}".R2.bowtie_endtoend.bam" -o \$name".bwt2.bam" -q 0 -t bowtie.stat
  """
  }

process fixMDtags {
  tag { bam.name.replaceFirst(".bt2.bam","") }

  //publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), val(rep), val(bio), path(bam))

  output:
  tuple(val(sample), val(rep), val(bio), path('*qSort.bam'), emit:bam)

  script:
  def bamname = bam.name.replaceFirst(".bam","")
  """
  ## NOTE LENIENT IS REQUIRED BECAUSE THE FLAGS INTRODUCED BY HiCPro BOWTIE ALIGNMENT
  ## ARE FLAGGED AS ERRORS BY PICARD

  java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar /usr/picard/picard.jar SortSam \
                -SO coordinate \
                -I ${bam} \
                --VALIDATION_STRINGENCY LENIENT \
                -O ${bamname}.sorted.bam

  java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar /usr/picard/picard.jar SetNmMdAndUqTags \
                -R ${params.fasta} \
                -I ${bamname}.sorted.bam \
                -O ${bamname}.MDfixed.bam \
                --VALIDATION_STRINGENCY LENIENT

  java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar /usr/picard/picard.jar SortSam \
                -SO queryname \
                -I ${bamname}.MDfixed.bam \
                -O ${bamname}.qSort.bam \
                --VALIDATION_STRINGENCY LENIENT
  """
  }

process markAlleleOfOrigin {
  tag { bam.name.replaceFirst(".bt2.qSort.bam","") }

  publishDir "${params.outdir}/bam",     mode: 'copy', overwrite: true, pattern: "*phased.bam", enabled: params.saveBAM
  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: "*markallelicstatus.txt"

  input:
  tuple(val(sample), val(rep), val(bio), path(bam))

  output:
  tuple(val(sample), val(rep), val(bio), path('*bam'), emit: bam)
  path('*markallelicstatus.txt', emit: rep)
  tuple(val("${sample}_${bio}"), path('*allelic_stats.txt'), emit: report)

  script:
  def outbam = bam.name.replaceFirst(".bam",".phased.bam")
  def tmprep = bam.name.replaceFirst(".bam",".phased.allelstat")
  def outrep = bam.name.replaceFirst(".bam",".phased.markallelicstatus.txt")
  def outmqc = bam.name.replaceFirst(".bam",".allelic_stats.txt")
  """
  echo "${sample}${rep}${bio}"

  ##Bugs in reporting - use custom script instead
  ##/HiC-Pro_3.0.0/scripts/markAllelicStatus.py -i ${bam} -s ${params.vcf} -r |samtools view -Shb - > ${outbam}

  python /usr/local/hicpro/scripts/markAllelicStatus.py -i ${bam} -s ${params.vcf} -r -o ${outbam}

  grep    allelic_status ${tmprep}                      >${outmqc}
  grep -v allelic_status ${tmprep} |grep -v FOR-MULTIQC >${outrep}
  """
  }

process pairsam {

  tag { bam.name.replaceFirst(".bt2.qSort.phased.bam","") }

  input:
  tuple(val(sample), val(rep), val(bio), path(bam))

  output:
  tuple(val("${sample}_${bio}"), val(sample), val(bio), path('*.sam.pairs.gz'), emit: pairs)
  tuple(val("${sample}_${bio}"), val(sample), val(bio), path('*.stat.txt'), emit: rep)
  script:
  def initpairs   = bam.name.replaceFirst(".bam",".initpairs.gz")
  def initpairsOK = bam.name.replaceFirst(".bam",".initpairsOK")
  def pairs       = bam.name.replaceFirst(".bam",".sam.pairs.gz")
  def pairstat    = bam.name.replaceFirst(".bam",".sam.pairsam.stat.txt")
  def samcols     = (params.phased)?"mapq,XA,MD":"mapq"
  """
  #!/bin/bash

  # Classify Hi-C molecules as unmapped/single-sided/multimapped/chimeric/etc
  # and output one line per read, containing the following, separated by \\v:
  #  * true-flipped pairs
  #  * read id
  #  * type of a Hi-C molecule
  #  * corresponding sam entries

  samtools view -h ${bam} | pairtools parse -c ${params.fai} \
    --output-stats ${pairstat} \
    --add-columns ${samcols} \
    -o ${initpairs}

  ## The output of pairtools parse does not account for missing fields. This messes
  ## things up later. Here, we do a small cleanup, replacing missing fields with NA

  ## First, get the header - use perl to prevent the need to go through whole file with grep
  zcat ${initpairs} |perl -lane 'exit if (\$_ !~ /^#/);print \$_' >${initpairsOK}

  ## Get the number of expected columns
  ncols=`grep ^#columns ${initpairsOK} |perl -lane 'print \$#F'`

  ## Pad out gaps with NA
  zcat ${initpairs}|grep -v ^# |
      perl  -pe 's/\\t(?=\\t)/\\tNA/g' | ##This line does adjacent gaps
      perl -lane 'print \$_.("\\tNA" x ('\$ncols'-\$#F-1))' >>${initpairsOK} ## This line does gaps at the end of the line

  cat ${initpairsOK} | pairtools sort --nproc ${task.cpus} \
    --compress-program lz4c \
    --tmpdir ${TMPDIR} \
    --output ${pairs}
  """
  }

process mergesampairs {

  tag { grp }

  publishDir "${params.outdir}/dedup",   mode: 'copy', overwrite: true, pattern: '*dedup.pair*'

  input:
  tuple(val(grp), val(sample), val(bio), path(pairsGZ))

  output:
  //tuple(val(sample), val(grp), path('*.dedup.all.pairs.gz'), path('*.dedup.all.pairs.gz.px2'))
  tuple(val(sample), val(grp), path('*.all.pairs.gz'), path('*.all.pairs.gz.px2'))

  script:
  """
  #!/bin/bash
  ### MERGE SEQUENCING REPLICATES AND
  ### MARK DUPLICATES
  nmax=`ls *pairs.gz |wc -l`

  pairtools merge --max-nmerge \$nmax --nproc ${task.cpus} --memory  ${task.memory} --output ${grp}.merged.sam.pairs.gz ${pairsGZ}
  pairtools dedup --mark-dups --output-dups - --output-unmapped - --output ${grp}.all.pairs.gz ${grp}.merged.sam.pairs.gz

  pairix ${grp}.all.pairs.gz
  """
  }

process filter_pairs {

  tag { grp }

  publishDir "${params.outdir}/filtered",   mode: 'copy', overwrite: true, pattern: '*filtered.pair*'

  input:
  tuple(val(sample), val(grp), path(pairs), path(idx))

  output:
  tuple(val(sample), path('*.filtered.pairs.gz'), path('*.filtered.pairs.gz.px2'))

  script:
  """
  #!/bin/bash
  ### FILTER PAIRSSAM FILE
  ## Generate lossless bam
  pairtools split --output-sam ${grp}.lossless.bam ${pairs}

  # Select UU, UR, RU reads
  pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
      --output-rest ${grp}.unmapped.sam.pairs.gz \
      --output temp.gz \
      ${pairs}

  pairtools split --output-pairs temp1.gz temp.gz
  pairtools select 'True' --chrom-subset ${params.fai} -o ${grp}.filtered.pairs.gz temp1.gz
  pairix ${grp}.filtered.pairs.gz

  """
  }

process pairtools_stats {

  tag { grp }

  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

  input:
  tuple(val(id), val(grp), path(pairsGZ))

  output:
  path("*pairtools_report*", emit: report, optional: true)

  script:
  def pairs_report="${grp}.pairtools_report.txt";
  """
  #!/bin/bash
  ## AVOID GENERATING PAIRS REPORTS FOR LESS THAN 1000 READS
  ## CAN BE PROBLEMATIC
  nreads=`zcat ${pairsGZ} |wc -l`
  if [[ "\$nreads" -gt 1000 ]]; then
    pairtools stats -o ${pairs_report} ${pairsGZ}
  fi
  """
  }

process splitByPhase {

  tag { grp }

  publishDir "${params.outdir}/dedup",  mode: 'copy', overwrite: true, pattern: '*pairs*'
  publishDir "${params.outdir}/report", mode: 'copy', overwrite: true, pattern: '*pairstats.txt'

  input:
  tuple(val(sample), val(grp), path(pairs), path(pairsIDX))

  output:
  tuple(val(sample), val(grp), path('*.pairs.gz'), path('*.pairs.gz.px2'), emit: pairs)
  tuple(val(grp), path('*.pairstats.txt'), emit: report)

  script:
  def pairsGT1="${grp}.hom1.pairs.gz"
  def pairsGT2="${grp}.hom2.pairs.gz"
  def pairsHET="${grp}.het.pairs.gz"
  def pairsHOM="${grp}.hom.pairs.gz"
  def pairs00="${grp}.00.gz"
  def pairs01="${grp}.01.gz"
  def pairs02="${grp}.02.gz"
  def pairsSTAT="${grp}.pairstats.txt"
  """
  #!/bin/bash
  pairtools select 'XA1+XA2 == "11"'            ${pairs} --output ${pairsGT1}
  pairtools select 'XA1+XA2 == "22"'            ${pairs} --output ${pairsGT2}
  pairtools select 'XA1+XA2 in ("12","21")'     ${pairs} --output ${pairsHET}
  pairtools select 'XA1+XA2 not in ("12","21")' ${pairs} --output ${pairsHOM}
  pairtools select 'XA1+XA2 in ("00")'          ${pairs} --output ${pairs00}
  pairtools select 'XA1+XA2 in ("01","10")'     ${pairs} --output ${pairs01}
  pairtools select 'XA1+XA2 in ("02","20")'     ${pairs} --output ${pairs02}

  for pairsGZ in ${pairsGT1} ${pairsGT2} ${pairsHET} ${pairsHOM}; do
    n=`zcat \$pairsGZ |grep -v ^# |wc -l`
    if [[ "\$n" -gt 0 ]]; then
      pairix \$pairsGZ
    else
      echo "**** WARNING **** : No valid pairs for \$pairsGZ"
      rm \$pairsGZ
    fi
  done

  ## Pair types report
  n00=`zcat ${pairs00} |grep -v ^# |wc -l`; echo "allelic_pairs/0-0\\t\$n00" >>${pairsSTAT}
  n01=`zcat ${pairs01} |grep -v ^# |wc -l`; echo "allelic_pairs/0-1\\t\$n01" >>${pairsSTAT}
  n02=`zcat ${pairs02} |grep -v ^# |wc -l`; echo "allelic_pairs/0-2\\t\$n02" >>${pairsSTAT}
  n11=`zcat ${pairsGT1} |grep -v ^# |wc -l`; echo "allelic_pairs/1-1\\t\$n11" >>${pairsSTAT}
  n12=`zcat ${pairsHET} |grep -v ^# |wc -l`; echo "allelic_pairs/1-2\\t\$n12" >>${pairsSTAT}
  n22=`zcat ${pairsGT2} |grep -v ^# |wc -l`; echo "allelic_pairs/2-2\\t\$n22" >>${pairsSTAT}
  """
  }

process pairsQC {

  tag { grp }

  publishDir "${params.outdir}/reports",     mode: 'copy', overwrite: true

  input:
  tuple(val(sample), val(grp), path(pairs), path(idx))

  output:
  path('*.zip')

  script:
  """
  scriptdir="/usr/local/bin/pairsqc/"

  #if [ "${params.genome}" == "mm10" ]; then
  #  max_distance=8.2
  #else
  #  max_distance=8.4
  #fi

  sort -k1,1 -V -s ${params.fai} >chrom_sizes.txt

  max_distance=`sort -k2n,2n chrom_sizes.txt |tail -n1 |cut -f2 |perl -lane 'print sprintf("%2.1f",((log(\$_)/log(10))-0.05))'`

  python3 \$scriptdir/pairsqc.py -p ${pairs} -c chrom_sizes.txt -tP -s ${grp} -O ${grp} -M \$max_distance

  ## NOTE : Should be coded properly ... assumes a 4-cutter (4)
  Rscript \$scriptdir/plot.r 4 ${grp}_report

  zip -r ${grp}_pairsQC_report.zip ${grp}_report

  """
  }

process get_RE_file {

  tag { params.re }

  publishDir "${params.outdir}/re", mode: 'copy', overwrite: true

  output:
  path('*.txt')

  script:
  """
  #!/bin/bash
  python /usr/local/bin/juicer/misc/generate_site_positions.py ${params.re} genome ${params.fasta}

  """
  }

process merge_biological_replicates {

  tag { sample }

  publishDir "${params.outdir}/pairs", mode: 'copy', overwrite: true, enabled: params.dnase

  input:
  tuple(val(sample), path(pairs), path(idx))

  output:
  tuple(val(sample), path('*.pairs.gz'),path('*.pairs.gz.px2'))

  script:
  """
  #!/bin/bash

  n=`ls *.pairs.gz | wc -l`

  if [ \$n -eq 1 ]
  then
      cp ${pairs} ${sample}.pairs.gz
      pairix -f ${sample}.pairs.gz
  else
      # unzipping to named pipes
      arg=''
      k=1
      for f in *.pairs.gz
      do
        mkfifo pp.\$k
        arg="\$arg pp.\$k"
        gunzip -c \$f | grep -v '^#' > pp.\$k &
        let "k++"
      done

      # header
      gunzip -c ${pairs[0]} | grep "^#" | grep -v '^#command:'  > ${sample}.pairs

      for f in *.pairs.gz
      do
        gunzip -c \$f | grep '^#command:' >> ${sample}.pairs
      done

      # merging
      sort -m -k2,2 -k4,4 -k3,3g -k5,5g \$arg >> ${sample}.pairs

      # compressing
      bgzip -f ${sample}.pairs

      # indexing
      pairix -f ${sample}.pairs.gz

  fi
  """
  }

process addfrag2pairs{
  tag { sample }

  publishDir "${params.outdir}/pairs", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), path(pairs), path(idx), path(re))

  output:
  tuple(val(sample), path('*.ff.pairs.gz'),path('*.ff.pairs.gz.px2'))

  script:
  """
  #!/bin/bash
  if [ "${params.dnase}" == "true" ]; then
    cp ${pairs} ${sample}.ff.pairs.gz
    cp ${idx}   ${sample}.ff.pairs.gz.px2
  else
    run-addfrag2pairs.sh -i ${pairs} -r ${re} -o ${sample}
  fi
  """
  }

process run_cooler{
  tag { sample }

  publishDir "${params.outdir}/cooler", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), path(pairs), path(idx))

  output:
  tuple(val(sample),path('*.cool'))

  script:
  """
  #!/bin/bash
  sort -k1,1 -V -s ${params.fai} >chrom_sizes.txt

  cooler cload pairix -p ${task.cpus} -s 10 chrom_sizes.txt:1000 ${pairs} ${sample}.cool

  """
  }

process balance_cool_matrix{
  tag { sample }

  publishDir "${params.outdir}/cooler", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), path(cool))

  output:
  tuple(val(sample),path('*.balanced.cool'))

  script:
  def balanced_cool=cool.name.replaceFirst(".cool",".balanced.cool")
  """
  cp ${cool} ${bal}

  if [ -f "${params.blacklist}" ]; then
    grep -vP '(rand|Un)' ${params.blacklist} |sort -k1,1 -k2n,2n -V >bl.bed
    cooler balance -p ${task.cpus} --blacklist bl.bed --max-iters 500 --force ${bal}
  else
    cooler balance -p ${task.cpus} --max-iters 500 --force ${bal}
  fi

  """
  }

process call_compartments_from_cool{
  tag { cool }

  publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), path(cool))

  output:
  tuple(val(sample),path('*.balanced.cool'))

  script:
  def balanced_cool=cool.name.replaceFirst(".cool",".compartments")
  """
  cooltools call-compartments -v --bigwig -o ${compartments} ${cool}
  """
  }

process run_juicebox_pre{
  tag { sample }

  publishDir "${params.outdir}/juicer", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), path(pairs), path(idx))

  output:
  tuple(val(sample),path('*hic'))

  script:
  """
  #!/bin/bash

  sort -k1,1 -V -s ${params.fai} |cut -f1,2 >chrom_sizes.txt

  java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar /usr/local/bin/juicer_tools.jar pre \
            -n ${pairs} ${sample}.hic \
             -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 \
             chrom_sizes.txt -q 0

  java -Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g -jar /usr/local/bin/juicer_tools.jar addNorm \
            -w 1000 -d -F ${sample}.hic

  """
  }

process cool2multirescool{
  tag { sample }

  publishDir "${params.outdir}/cooler", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), path(cool))

  output:
  tuple(val(sample),path('*cool*'))

  script:
  """
  #!/bin/bash

  run-cool2multirescool.sh \
     -i ${cool} \
     -p ${task.cpus} \
     -o ${sample} \
     -c 10000000 \
     -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000

  """
  }

process hicnormvector_to_mcool{
  tag { sample }

  publishDir "${params.outdir}/cooler", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), path(cool_file), path(hic_file))

  output:
  tuple(val(sample),path('*mcool*'))

  script:
  def outcool=cool_file.name.replaceFirst("cool","mcool")
  """
  #!/bin/bash
  cp ${cool_file} ${outcool}

  scriptdir=/usr/local/bin/
  hic2cool extract-norms -e ${hic_file} ${outcool}
  """
  }

process concatenate_phasing_reports {
  tag { id }

  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

  input:
  tuple(val(id), path(reports))

  output:
  path('*phasing_report.txt')

  script:
  """
  #!/bin/bash
  for v in `cut -f1 *stats.txt |sort |uniq`; do n=`grep \$v *txt |cut -f2 |
       awk -F',' '{sum+=\$1;} END{print sum;}'`; echo -e "\$v\\t\$n"; done |sort -k1,1 >${id}.phasing_report.txt
  """
  }

process multiqc{
  tag { reports[0] }

  publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

  input:
  path(reports)

  output:
  path('*multiqc*')

  script:
  """
  #!/bin/bash
  nm=`ls *.pairtools_report.txt`
  name=\${nm/.all.pairtools_report.txt/.multiqc}

  for z in *zip; do
    cp \$z report.zip
    unzip report.zip
  done

  multiqc -n \$name .
  """
  }
