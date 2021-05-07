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

  withLabel: getFQ{
    cpus = { 4 }
    memory = { 8.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName: mergeFQ{
    cpus = { 4 }
    memory = { 8.GB }
    time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:commitToObj{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 1.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:fastqC{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 1.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:fastqScreen{
    cpus = { 1 }
    memory = { 8.GB }
    time = { 1.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }
}
