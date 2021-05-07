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

  withName:multiQC{
    cpus = { 1 }
    memory = { 4.GB }
    time = { 1.hour * task.attempt * process_length_factor }
    container = "quay.io/biocontainers/multiqc:1.10.1--py_0"
  }
}
