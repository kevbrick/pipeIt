process multiQC {

  input:
  path(reports)

  output:
  path('*ultiQC*', emit: mqcReport)

  script:
  """
  multiqc -f -n ${params.name}.multiQC .
  """
  }
