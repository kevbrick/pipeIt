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
  withName:ssdsAlign{
    cpus = { 12 }
    memory = { 32.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:mergeBAMssds{
    cpus = { 4 }
    memory = { 32.GB }
    time = { 16.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:parseITRs{
    cpus = { 4 }
    memory = { 12.GB }
    time = { 6.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:gatherITROutputs{
    cpus = { 4 }
    memory = { 12.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:makeFRBWssds{
    cpus = { 2 }
    memory = { 12.GB }
    //Defined in process
    //time = { 2.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:makeSSreport{
    cpus = { 1 }
    memory = { 4.GB }
    time = { 3.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }

  withName:multiQCssds{
    cpus = { 1 }
    memory = { 4.GB }
    time = { 1.hour * task.attempt * process_length_factor }
    container = "$DSL2DIR/singularity/RDCO.sif"
  }
}
