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

  withName:getCoverage{
    cpus = { 2 }
    memory = { 32.GB }
    time = { 3.hour * task.attempt * process_length_factor }
    //container = "$DSL2DIR/singularity/RDCO.sif"
    container = "docker://kevbrick/rtseq:2.0"
  }

  withName:mergeCoverageBedgraphs{
    cpus = { 2 }
    memory = { 16.GB }
    time = { 2.hour * task.attempt * process_length_factor }
    //container = "$DSL2DIR/singularity/RDCO.sif"
    container = "docker://kevbrick/rtseq:2.0"
  }

  withName:generateMpileup{
    cpus = { 2 }
    memory = { 32.GB }
    time = { 9.hour * task.attempt * process_length_factor }
    //container = "$DSL2DIR/singularity/RDCO.sif"
    container = "docker://kevbrick/rtseq:2.0"
  }

  withName:getGCstats{
    cpus = { 2 }
    memory = { 32.GB }
    time = { 1.hour * task.attempt * process_length_factor }
    //container = "$DSL2DIR/singularity/RDCO.sif"
    container = "docker://kevbrick/rtseq:2.0"
  }

  withName:getGCcorrFile{
    cpus = { 1 }
    memory = { 4.GB }
    time = { 0.25.hour * task.attempt * process_length_factor }
    //container = "$DSL2DIR/singularity/RDCO.sif"
    container = "docker://kevbrick/rtseq:2.0"
  }
}
