process_length_factor = 1

profiles{
  test{
    process_length_factor = 0.1
  }

  singlecell{
    process_length_factor = 0.25
  }
}

process {
  withName:bwa4D{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 16.h * task.attempt * process_length_factor}
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:bowtie2_end2end{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 16.h * task.attempt * process_length_factor}
    container = 'docker://nservant/hicpro:latest'
  }

  withName:trim_hic_reads{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 16.h * task.attempt * process_length_factor}
    container = 'docker://nservant/hicpro:latest'
  }

  withName:bowtie2_on_trimmed_reads{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 16.h * task.attempt * process_length_factor }
    container = 'docker://nservant/hicpro:latest'
  }

  withName:bowtie2_mergeR1R2{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 16.h * task.attempt * process_length_factor }
    container = 'docker://nservant/hicpro:latest'
  }

  withName:bowtie2_make_paired_bam{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 16.h * task.attempt * process_length_factor }
    container = 'docker://nservant/hicpro:latest'
  }

  withName:fixMDtags{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 8.h * task.attempt * process_length_factor }
    container = "docker://broadinstitute/picard:2.25.0"
  }

  withName:markAlleleOfOrigin{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 16.h * task.attempt * process_length_factor }
    container = 'docker://kevbrick/markallelicstatus:1.0'
    //container = 'docker://nservant/hicpro:latest'
  }

  withName:pairsam{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 8.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:mergepairs_dedup{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 8.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:filter_pairs{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 8.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:pairtools_stats{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 8.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:splitByPhase{
    cpus = { 4 }
    memory = { 16.GB }
    time = { 6.hour * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:pairs_to_bam{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 4.hour * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:clean_dnase_bam{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 4.hour * task.attempt * process_length_factor }
    container = "quay.io/biocontainers/pysam:0.15.2--py36h02877da_7"
  }

  withName:make_mnase_bigwig{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 4.hour * task.attempt * process_length_factor }
    container = "quay.io/biocontainers/deeptools:3.5.1--py_0"
  }

  // withName:perlSplitByPhase{
  //   cpus = { 4 }
  //   memory = { 16.GB }
  //   time = { 6.hour * task.attempt * process_length_factor }
  //   container = "docker://duplexa/4dn-hic:v43"
  // }

  withName:mergepairs_nodedup{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 8.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:get_RE_file{
    cpus = { 8 }
    memory = { 8.GB }
    time = { 2.h * task.attempt  }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:pairsQC{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 4.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:addfrag2pairs{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 8.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:run_cooler{
    cpus = { 32 }
    memory = { 64.GB }
    time = { 48.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:balance_cool_matrix{
    cpus = { 16 }
    memory = { 16.GB }
    time = { 10.h * task.attempt * process_length_factor }
    container = "docker://kevbrick/cooltools:1.0"
  }

  withName:call_compartments_from_cool{
    cpus = { 32 }
    memory = { 1500.GB }
    time = { 12.h * task.attempt * process_length_factor }
    container = "docker://kevbrick/cooltools:1.0"
  }

  withName:run_juicebox_pre{
    cpus = { 8 }
    memory = { 64.GB }
    time = { 24.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:cool2multirescool{
    cpus = { 16 }
    memory = { 32.GB }
    time = { 8.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:hicnormvector_to_mcool{
    cpus = { 8 }
    memory = { 32.GB }
    time = { 16.h * task.attempt * process_length_factor }
    container = "docker://duplexa/4dn-hic:v43"
  }

  withName:concatenate_phasing_reports{
    cpus = { 1 }
    memory = { 4.GB }
    time = { 0.25.h * task.attempt }
  }

  withName:multiqc{
    cpus = { 2 }
    memory = { 8.GB }
    time = { 2.h * task.attempt }
    container = "kevbrick/multiqc_pairsqc:1.0"
  }
}
