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
  withName:bamCoverage{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 2.hour * task.attempt * process_length_factor }
    container = "quay.io/biocontainers/deeptools:3.5.1--py_0"
  }

  withName:computeMatrixRefPoint{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 2.hour * task.attempt * process_length_factor }
    container = "quay.io/biocontainers/deeptools:3.5.1--py_0"
  }

  withName:computeMatrixRefPointMultiple{
    cpus = { 12 }
    memory = { 16.GB }
    time = { 2.hour * task.attempt * process_length_factor }
    container = "quay.io/biocontainers/deeptools:3.5.1--py_0"
  }
}
