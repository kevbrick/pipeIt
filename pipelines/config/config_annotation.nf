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
  withName:getB6DSBhotspots{
    cpus = { 8 }
    memory = { 16.GB }
    time = { 2.hour * task.attempt * process_length_factor }
    container = "docker://kevbrick/getb6xcastf1dsbhotspots:1.0"
  }
}
