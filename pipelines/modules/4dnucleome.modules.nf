// KB: March 11 2021
// These modules are based on the 4D nucleome protocols
// Detailed information is at:
// https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline#recheck

process bwa4D {
  tag { fq1 }

  publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true, pattern: '*bam', enabled: params.saveBAM

  input:
  tuple(val(id), val(sample), val(rep), val(bio), path(fq1), path(fq2))

  output:
  tuple(val(sample), val(rep), val(bio), path('*BWAinit.bam'))

  script:
  def nm = fq1.name.replaceAll("R1.(\\d+).fastq.gz","\$1")
  //println txt.replaceAll("(kevin)_(\\d+)_","\$2--\$1")
  """
  #ln -s ${fq1} fastq1.gz
  #ln -s ${fq2} fastq2.gz

  bwa mem \
    -t ${task.cpus} \
    -SP5M \
    ${params.bwaidx} \
    ${fq1} ${fq2} | samtools view -Shb -o ${nm}.BWAinit.bam -
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

  publishDir "${params.outdir}/bam", mode: 'copy', overwrite: true, pattern: '*bam', enabled: params.saveBAM

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
  //publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true, pattern: "*markallelicstatus.txt"

  input:
  tuple(val(sample), val(rep), val(bio), path(bam))

  output:
  tuple(val(sample), val(rep), val(bio), path('*bam'), emit: bam)
  tuple(val("${sample}_${bio}"), path('*markallelicstatus.txt'), emit: rep)
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
  tuple(val({ (bio == 'NA') ? "${sample}" : "${sample}_${bio}"} ), path('*.sam.pairs.gz'), emit: pairs)
  tuple(val({ (bio == 'NA') ? "${sample}" : "${sample}_${bio}"} ), path('*.stat.txt'), emit: rep)

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

process mergepairs_dedup {

  tag { grp }

  publishDir "${params.outdir}/pairs/dedup",   mode: 'copy', overwrite: true, pattern: '*pairs*', enabled: params.saveAllPairs

  input:
  tuple(val(grp), path(pairsGZ))

  output:
  //tuple(val(sample), val(grp), path('*.dedup.all.pairs.gz'), path('*.dedup.all.pairs.gz.px2'))
  tuple(val(grp), path('*.all.pairs.gz'), path('*.all.pairs.gz.px2'))

  script:
  """
  #!/bin/bash
  ### MERGE SEQUENCING REPLICATES AND
  ### MARK DUPLICATES
  nmax=`ls *pairs.gz |wc -l`

  pairtools merge --max-nmerge \$nmax --nproc ${task.cpus} --memory  ${task.memory} --output ${grp}.merged.sam.pairs.gz ${pairsGZ}
  pairtools dedup --mark-dups --output-dups - --output-unmapped - --output ${grp}.dedup.all.pairs.gz ${grp}.merged.sam.pairs.gz

  pairix ${grp}.dedup.all.pairs.gz
  """
  }

process pairs_to_bam {

  tag { grp }

  input:
  tuple(val(id), val(grp), path(pairsGZ))

  output:
  tuple(val(id), val(grp), path('*.bam'))

  script:
  def bam=pairsGZ.name.replaceFirst("pairs.gz","bam")
  """
  #!/bin/bash
  pairtools split ${pairsGZ} --output-pairs pairs.gz --output-sam pairs.sam
  samtools view -Shb pairs.sam >${bam}
  """
  }

process clean_dnase_bam {

  tag { grp }

  publishDir "${params.outdir}/dnase",   mode: 'copy', overwrite: true, pattern: '*dnase.bam*'

  input:
  tuple(val(sample), val(grp), path(pairsBAM))

  output:
  tuple(val(sample), val(grp), path('*.dnase.bam'), path('*.dnase.bam.bai'))

  shell:
  outBAM = pairsBAM.name.replaceFirst("bam","dnase.bam")
  tmpBAM = pairsBAM.name.replaceFirst("bam","tmp.bam")
  '''
  #!/usr/bin/env python
  import pysam
  import argparse
  import re

  inBAM  ="!{pairsBAM}"
  tmpBAM ="!{tmpBAM}"
  outBAM ="!{outBAM}"

  if re.search(".bam",inBAM):
      samfile = pysam.AlignmentFile(inBAM, "rb")
  else:
      if re.search(".sam",inBAM):
          samfile = pysam.AlignmentFile(inBAM, "r")

  bamfile = pysam.AlignmentFile(tmpBAM, "wb", template=samfile)

  for read in samfile:
      if read.is_paired:
          if (read.is_read1):
              read.is_paired=False
              read.is_proper_pair=False
              read.is_read1=False
              read.is_read2=False
              read.mate_is_reverse=False
              read.mate_is_unmapped=False
              read.next_reference_id=-1
              read.next_reference_start=-1
              read.template_length=-1

              bamfile.write(read)

  samfile.close()
  bamfile.close()

  pysam.sort("-o",outBAM,tmpBAM)
  pysam.index(outBAM)

  '''
  }

process make_mnase_bigwig {

  tag { grp }
  publishDir "${params.outdir}/dnase",   mode: 'copy', overwrite: true, pattern: '*bigwig'

  input:
  tuple(val(id), val(grp), path(bam), path(bai))

  output:
  path('*bigwig', emit: bigwig)

  script:
  def bigwig=bam.name.replaceFirst("bam","bigwig")
  """
  #!/bin/bash
  bamCoverage -b ${bam} --centerReads -p ${task.cpus} --minMappingQuality 20 --ignoreDuplicates --extendReads 150 --centerReads --binSize 10 -o ${bigwig}
  """
  }

process filter_pairs {

  tag { grp }

  publishDir "${params.outdir}/pairs/filtered",   mode: 'copy', overwrite: true, pattern: '*filtered.pair*'

  input:
  tuple(val(grp), path(pairs), path(idx))

  output:
  tuple(val(grp), path('*.filtered.pairs.gz'), path('*.filtered.pairs.gz.px2'))

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
  tuple(val(grp), path(pairsGZ))

  output:
  path("*pairtools_report*", emit: reports, optional: true)

  script:
  //def pairs_report=pairsGZ.name.replaceFirst("pairs.gz","pairtools_report.txt");
  """
  #!/bin/bash
  ## AVOID GENERATING PAIRS REPORTS FOR LESS THAN 1000 READS
  ## CAN BE PROBLEMATIC
  nreads=`zcat ${pairsGZ} |wc -l`
  if [[ "\$nreads" -gt 1000 ]]; then
    pairtools stats -o ${grp}.pairtools_report.txt ${pairsGZ}
  fi
  """
  }

process splitByPhase {

  tag { grp }

  publishDir "${params.outdir}/pairs/phased",  mode: 'copy', overwrite: true, pattern: '*pairs*'
  publishDir "${params.outdir}/reports"      , mode: 'copy', overwrite: true, pattern: '*pairstats.txt'

  input:
  tuple(val(grp), path(pairs), path(pairsIDX))

  output:
  tuple(val(grp), path('*.phased.*.pairs.gz'), path('*.phased.*.pairs.gz.px2'), emit: pairs)
  tuple(val(grp), path('*.pairstats.txt'), emit: report)

  script:
  def pairsGT1="${grp}.phased.hom1.pairs.gz"
  def pairsGT2="${grp}.phased.hom2.pairs.gz"
  def pairsHET="${grp}.phased.het.pairs.gz"
  def pairsHOM="${grp}.phased.hom.pairs.gz"
  def pairs00="${grp}.00.gz"
  def pairs01="${grp}.01.gz"
  def pairs02="${grp}.02.gz"
  def pairsSTAT="${grp}.pairstats.txt"
  """
  #!/bin/bash
  pairtools select '(chrom1 != "!") and (chrom2 != "!")' ${pairs}         --output aligned.pairs.gz
  pairtools select 'XA1+XA2 == "11"'                     aligned.pairs.gz --output ${pairsGT1}
  pairtools select 'XA1+XA2 == "22"'                     aligned.pairs.gz --output ${pairsGT2}
  pairtools select 'XA1+XA2 in ("12","21")'              aligned.pairs.gz --output ${pairsHET}
  pairtools select 'XA1+XA2 not in ("12","21")'          aligned.pairs.gz --output ${pairsHOM}
  pairtools select 'XA1+XA2 in ("00")'                   aligned.pairs.gz --output ${pairs00}
  pairtools select 'XA1+XA2 in ("01","10")'              aligned.pairs.gz --output ${pairs01}
  pairtools select 'XA1+XA2 in ("02","20")'              aligned.pairs.gz --output ${pairs02}

  for pairsGZ in ${pairsGT1} ${pairsGT2} ${pairsHET} ${pairsHOM}; do
    n=`zcat \$pairsGZ |grep -v ^# |head -n 1000 |wc -l`
    if [[ "\$n" -gt 0 ]]; then
      pairix \$pairsGZ
    else
      echo "**** WARNING **** : No valid pairs for \$pairsGZ"
      rm \$pairsGZ
    fi
  done

  ## Pair types report
  n00=`zcat ${pairs00}  |grep -v ^# |wc -l`; echo -e "allelic_pairs/0-0\\t\$n00" >>${pairsSTAT}
  n01=`zcat ${pairs01}  |grep -v ^# |wc -l`; echo -e "allelic_pairs/0-1\\t\$n01" >>${pairsSTAT}
  n02=`zcat ${pairs02}  |grep -v ^# |wc -l`; echo -e "allelic_pairs/0-2\\t\$n02" >>${pairsSTAT}
  n11=`zcat ${pairsGT1} |grep -v ^# |wc -l`; echo -e "allelic_pairs/1-1\\t\$n11" >>${pairsSTAT}
  n12=`zcat ${pairsHET} |grep -v ^# |wc -l`; echo -e "allelic_pairs/1-2\\t\$n12" >>${pairsSTAT}
  n22=`zcat ${pairsGT2} |grep -v ^# |wc -l`; echo -e "allelic_pairs/2-2\\t\$n22" >>${pairsSTAT}
  """
  }

  // This PROCESS generates the same report but only runs through the file once.
  // If runtime becomes an issue, we can use this, but for now, too much
  // value in using established tools
  //--------------------------------------------------------------------------------------------------------------
  // process perlSplitByPhase {
  //
  //   tag { grp }
  //
  //   publishDir "${params.outdir}/pairs/perlphased",  mode: 'copy', overwrite: true, pattern: '*pairs*'
  //   publishDir "${params.outdir}/reports"      , mode: 'copy', overwrite: true, pattern: '*pairstatsperl.txt'
  //
  //   input:
  //   tuple(val(sample), val(grp), path(pairs), path(pairsIDX))
  //
  //   output:
  //   tuple(val(sample), val(grp), path('*.phased.*.pairs.gz'), path('*.phased.*.pairs.gz.px2'), emit: pairs)
  //   tuple(val(grp), path('*.pairstatsperl.txt'), emit: report)
  //
  //   shell:
  //   '''
  //   #!/usr/bin/env perl
  //   use strict;
  //
  //   my ($name, $pairsGZ) = ("!{grp}","!{pairs}");
  //
  //   open my $ref, '>', $name.'.phased.hom1.pairs';
  //   open my $alt, '>', $name.'.phased.hom2.pairs';
  //   open my $het, '>', $name.'.phased.het.pairs';
  //   open my $hom, '>', $name.'.phased.hom.pairs';
  //   open my $r00, '>', $name.'.phased.00.pairs';
  //   open my $r01, '>', $name.'.phased.01.pairs';
  //   open my $r02, '>', $name.'.phased.02.pairs';
  //
  //   open REPORT, '>', $name.".pairstatsperl.txt";
  //
  //   my ($colXA1, $colXA2, $colChrom1, $colChrom2, @colnames, %countReport);
  //
  //   open my $IN, "-|", "gunzip -c $pairsGZ";
  //
  //   while (<$IN>){
  //       chomp;
  //
  //       if ($_ =~ /^\\#columns:\\s+(.+)$/){
  //           my $cols = $1;
  //           @colnames=split(/\\s+/,$cols);
  //           for my $n(0..$#colnames){
  //               $colXA1    = $n if ($colnames[$n] eq 'XA1');
  //               $colXA2    = $n if ($colnames[$n] eq 'XA2');
  //               $colChrom1 = $n if ($colnames[$n] eq 'chrom1');
  //               $colChrom2 = $n if ($colnames[$n] eq 'chrom2');
  //           }
  //
  //           print $ref join("\\t",'#samheader:','@PG',
  //                           'ID:custom_parsing_script_like_pairtools_select',
  //                           'PN:custom_parsing_script_like_pairtools_select',
  //                           'CL:filter XA1=1 and XA2=1')."\\n";
  //
  //           print $alt join("\\t",'#samheader:','@PG',
  //                           'ID:custom_parsing_script_like_pairtools_select',
  //                           'PN:custom_parsing_script_like_pairtools_select',
  //                           'CL:filter XA1=2 and XA2=2')."\\n";
  //
  //           print $het join("\\t",'#samheader:','@PG',
  //                           'ID:custom_parsing_script_like_pairtools_select',
  //                           'PN:custom_parsing_script_like_pairtools_select',
  //                           'CL:filter ((XA1=1 and XA2=2) or (XA1=2 and XA2=1))')."\\n";
  //
  //           print $hom join("\\t",'#samheader:','@PG',
  //                           'ID:custom_parsing_script_like_pairtools_select',
  //                           'PN:custom_parsing_script_like_pairtools_select',
  //                           'CL:filter !((XA1=1 and XA2=2) or (XA1=2 and XA2=1))')."\\n";
  //
  //           print $r00 join("\\t",'#samheader:','@PG',
  //                           'ID:custom_parsing_script_like_pairtools_select',
  //                           'PN:custom_parsing_script_like_pairtools_select',
  //                           'CL:filter XA1=0 and XA2=0')."\n";
  //
  //           print $r01 join("\\t",'#samheader:','@PG',
  //                           'ID:custom_parsing_script_like_pairtools_select',
  //                           'PN:custom_parsing_script_like_pairtools_select',
  //                           'CL:filter XA1=0 and XA2=1')."\\n";
  //
  //           print $r02 join("\\t",'#samheader:','@PG',
  //                           'ID:custom_parsing_script_like_pairtools_select',
  //                           'PN:custom_parsing_script_like_pairtools_select',
  //                           'CL:filter XA1=0 and XA2=2')."\\n";
  //       }
  //
  //       if ($_ =~ /^\\#/){
  //           print $ref $_."\\n";
  //           print $alt $_."\\n";
  //           print $het $_."\\n";
  //           print $hom $_."\\n";
  //           print $r00 $_."\\n";
  //           print $r01 $_."\\n";
  //           print $r02 $_."\\n";
  //           next;
  //       }
  //
  //       my @row = split(/\\t/,$_);
  //
  //       ## Only use mappable read pairs
  //       next unless ($row[$colChrom1] ne "!" && $row[$colChrom2] ne "!");
  //
  //       my $gt = join(":",$row[$colXA1],$row[$colXA2]);
  //
  //       if ($gt eq "1:1")                   {print $ref $_."\\n"; $countReport{"1-1"}++}
  //       if ($gt eq "2:2")                   {print $alt $_."\\n"; $countReport{"2-2"}++}
  //       if ($gt eq "1:2" || $gt eq "2:1")   {print $het $_."\\n"; $countReport{"1-2"}++}
  //       if (!($gt eq "1:2" || $gt eq "2:1")){print $hom $_."\\n"; }
  //       if ($gt eq "0:0")                   {print $r00 $_."\\n"; $countReport{"0-0"}++}
  //       if ($gt eq "1:0" || $gt eq "0:1")   {print $r01 $_."\\n"; $countReport{"0-1"}++}
  //       if ($gt eq "2:0" || $gt eq "0:2")   {print $r02 $_."\\n"; $countReport{"0-2"}++}
  //   }
  //
  //   for my $n("0-0","0-1","0-2","1-1","1-2","2-2"){
  //       print REPORT join("\\t","allelic_pairs/$n",$countReport{$n})."\\n";
  //   }
  //   close REPORT;
  //   close $ref;
  //   close $alt;
  //   close $het;
  //   close $hom;
  //   close $r00;
  //   close $r01;
  //   close $r02;
  //
  //   for my $f ("hom1","hom2","hom","het"){
  //
  //       my $pairs = $name.'.phased.'.$f.'.pairs';
  //       system("bgzip $pairs");
  //
  //       my $GZ    = $pairs.'.gz';
  //       system("pairix $GZ");
  //
  //   }
  //
  //   '''
  //   }

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

process mergepairs_nodedup {

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

  //publishDir "${params.outdir}/pairs", mode: 'copy', overwrite: true

  input:
  tuple(val(sample), path(pairs), path(idx), path(re))

  output:
  tuple(val(sample), path('*.final.pairs.gz'),path('*.final.pairs.gz.px2'))

  script:
  """
  #run-addfrag2pairs.sh -i ${pairs} -r ${re} -o ${sample}

  if [ "${params.dnase}" == "true" ]; then
    cp ${pairs} ${sample}.final.pairs.gz
    cp ${idx}   ${sample}.final.pairs.gz.px2
  else
    gunzip -c ${pairs} | /usr/local/bin/pairix/util/fragment_4dnpairs.pl -a - ${sample}.final.pairs ${re}
    bgzip -f ${sample}.final.pairs
    pairix -f ${sample}.final.pairs.gz
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

  #For some reason, this is EXTREMELY SLOOOOOOOOOOOOOOOOOOOOOOOOOW
  #cooler cload pairix -p ${task.cpus} -s 10 chrom_sizes.txt:1000 ${pairs} ${sample}.cool
  ## Use non-indexed version instead !
  cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 chrom_sizes.txt:1000 ${pairs} ${sample}.cool
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
  cp ${cool} ${balanced_cool}

  if [ -f "${params.blacklist}" ]; then
    grep -vP '(rand|Un)' ${params.blacklist} |sort -k1,1 -k2n,2n -V >bl.bed
    cooler balance -p ${task.cpus} --blacklist bl.bed --max-iters 500 --force ${balanced_cool}
  else
    cooler balance -p ${task.cpus} --max-iters 500 --force ${balanced_cool}
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

process concatenate_allelicstatus_reports {
  tag { id }

  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

  input:
  tuple(val(name), path(report))

  output:
  path('*allelicstatusreport.txt')

  shell:
  '''
  #!/usr/bin/perl
  my $nord = 0;
  my (%ord);

  open my $IN, "-|", "cat *markallelicstatus.txt";
  open OUT, '>', "!{name}.allelicstatusreport.txt";

  my ($totSNP, $totSNPN, $totREADS, $totN, $totOneN, $totUn);
  my ($nRef  ,   $nAlt,   $nUnN,   $nUnO,   $nOth,   $nConf, $tot, $cnttot);
  my ($cntRef, $cntAlt, $cntUnN, $cntUnO, $cntOth, $cntConf);

  while (<$IN>){
    chomp;

    if ($_ =~ /^##/ && $_ !~ /^##\s.+========/){
      $_ =~ s/^(.+?)\\..+bam$/$1/;
      $header .= $_."\n" unless ($headerDets{$_}++);
    }

    $totSNP   += $1 if ($_ =~ /Total number of snps loaded\\s+([\\d\\.]+)/);
    $totSNPN  ++    if ($_ =~ /Total number of snps loaded\\s+([\\d\\.]+)/);
    $totREADS += $1 if ($_ =~ /Total number of reads\\s+([\\d\\.]+)/);
    $totN     += $1 if ($_ =~ /Total number of \\'N\\'s\\s+([\\d\\.]+)/);
    $totOneN  += $1 if ($_ =~ /Total number of reads with at least one \\'N\\'\\s+([\\d\\.]+)/);
    $totUn    += $1 if ($_ =~ /Total number of unassigned reads\\s+([\\d\\.]+)/);

    if ($_ =~ /Number of reads assigned to ref genome\\s+([\\d\\.]+)\\s+([\\d\\.]+)/){
        $nRef += $1;
        $tot  += int($1/$2*100) if ($1>0);
        $cnttot++ if ($1>0);
        $cntRef++;
    }

    if ($_ =~ /Number of reads assigned to alt genome\\s+([\\d\\.]+)\\s+([\\d\\.]+)/){
        $nAlt += $1;
        $tot  += int($1/$2*100) if ($1>0 && int($2)>0);
        $cnttot++ if ($1>0);
        $cntAlt++;
    }

    if ($_ =~ /Number of unassigned N-containing reads\\s+([\\d\\.]+)\\s+([\\d\\.]+)/){
        $nUnN += $1;
        $tot  += int($1/$2*100) if ($1>0 && int($2)>0);
        $cnttot++ if ($1>0);
        $cntUnN++;
    }

    if ($_ =~ /Number of other type unassigned reads\\s+([\\d\\.]+)\\s+([\\d\\.]+)/){
        $nUnO += $1;
        $tot  += int($1/$2*100) if ($1>0 && int($2)>0);
        $cnttot++ if ($1>0);
        $cntUnO++;
    }

    if ($_ =~ /Number of conflicting reads\\s+([\\d\\.]+)\\s+([\\d\\.]+)/){
        $nConf += $1;
        $tot   += int($1/$2*100) if ($1>0 && int($2)>0);
        $cnttot++ if ($1>0);
        $cntConf++;
    }

    if ($_ =~ /Number of other reads\\s+([\\d\\.]+)\\s+([\\d\\.]+)/){
        $nOth += $1;
        $tot   += int($1/$2*100) if ($1>0 && int($2)>0);
        $cnttot++ if ($1>0);
        $cntOth++;
    }
  }

  print OUT $header;

  print OUT "## =========================\n";
  print OUT join("\t","Total number of snps loaded",$totSNP?int($totSNP/$totSNPN):"0")."\n";
  print OUT "## =========================\n";
  print OUT join("\t","Total number of reads",$totREADS?$totREADS:"0")."\n";
  print OUT join("\t","Number of \'N\'s",$totN?$totN:"0")."\n";
  print OUT join("\t","Number of reads with at least one \'N\'",$totOneN?$totOneN:"0")."\n";
  print OUT join("\t","Number of unassigned reads",$totUn?$totUn:"0")."\n";
  print OUT join("\t","Total number of snps loaded",$totSNP?$totSNP:"0")."\n";
  print OUT join("\t","## = PERCENTAGES ========================")."\n";
  print OUT join("\t","Number of reads assigned to ref genome", ($nRef?$nRef:"0"),($nRef?sprintf("%4.2f",$nRef/$totREADS*100):"0"))."\n";
  print OUT join("\t","Number of reads assigned to alt genome", ($nAlt?$nAlt:"0"),($nAlt?sprintf("%4.2f",$nAlt/$totREADS*100):"0"))."\n";
  print OUT join("\t","Number of unassigned N-containing reads",($nUnN?$nUnN:"0"),($nUnN?sprintf("%4.2f",$nUnN/$totREADS*100):"0"))."\n";
  print OUT join("\t","Number of other type unassigned reads",  ($nUnO?$nUnO:"0"),($nUnO?sprintf("%4.2f",$nUnO/$totREADS*100):"0"))."\n";
  print OUT join("\t","Number of conflicting reads",            ($nCon?$nCon:"0"),($nCon?sprintf("%4.2f",$nCon/$totREADS*100):"0"))."\n";
  print OUT join("\t","Number of other reads",                  ($nOth?$nOth:"0"),($nOth?sprintf("%4.2f",$nOth/$totREADS*100):"0"))."\n";

  '''
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
