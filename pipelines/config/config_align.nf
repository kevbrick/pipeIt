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

  errorStrategy = 'retry'
  maxRetries = 1

  withLabel: aligner{
    cpus = { 12 }
    memory = { 32.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withLabel: mergeBAM{
    cpus = { 8 }
    memory = { 32.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withLabel: mergeBAMsv2{
    cpus = { 8 }
    memory = { 32.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:getPicardMetrics{
    cpus = { 4 }
    memory = { 16.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    //container = "$DSL2DIR/singularity/RDCO.sif"
    container = 'quay.io/biocontainers/picard:2.25.3--hdfd78af_0'
  }

  withName:bamToBW{
    cpus = { 4 }
    memory = { 12.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withLabel: trimFQ{
    cpus = { 4 }
    memory = { 8.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    //container = "$DSL2DIR/singularity/RDCO.sif"
    container = "quay.io/biocontainers/trim-galore:0.6.6--hdfd78af_1"
  }

  withName:makeDeeptoolsBW{
    cpus = { 12 }
    memory = { 24.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:samStats{
    cpus = { 2 }
    memory = { 8.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    //container = "$DSL2DIR/singularity/RDCO.sif"
    container = "quay.io/biocontainers/samtools:1.3.1--h1b8c3c0_8"
  }
}
