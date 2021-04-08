process multiQC {
  cpus 1
  memory '4g'

  time '1h'

  input:
  path(reports)

  output:
  path('*ultiQC*', emit: mqcReport)

  script:
  """
  multiqc -f -n ${params.name}.multiQC .
  """
  }
