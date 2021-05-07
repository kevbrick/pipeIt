process bamCoverage {
    tag { id }

    publishDir "${params.outdir}/deeptools/bigwig",   mode: 'copy', overwrite: true, pattern: '*bigwig'

    input:
    tuple(val(id), path(bam), path(bai))

    output:
    path('*bigwig', emit: bigwig)

    script:
    def bigwig=bam.name.replaceFirst("bam","bigwig")
    """
    #!/bin/bash
    bamCoverage -b ${bam} --centerReads -p ${task.cpus} \
                --minMappingQuality 20 \
                --ignoreDuplicates \
                --extendReads 150 \
                --centerReads \
                --binSize 10 \
                -o ${bigwig}
    """
  }

process computeMatrixRefPoint {
  tag { "${type} : ${bigwig}" }

  publishDir "${params.outdir}/deeptools/matrix",  mode: 'copy', overwrite: true, pattern: "*matrix"
  publishDir "${params.outdir}/deeptools/matrix",  mode: 'copy', overwrite: true, pattern: "*tab"
  publishDir "${params.outdir}/deeptools/heatmap", mode: 'copy', overwrite: true, pattern: "*svg"
  publishDir "${params.outdir}/deeptools/heatmap", mode: 'copy', overwrite: true, pattern: "*png"

  input:
  tuple(path(bigwig), val(type), path(bed))

  output:
  path('*.matrix', emit: matrix)
  path('*.tab',    emit: tab)
  path('*.png',    emit: png)
  path('*.svg',    emit: svg)

  script:
  def id = bigwig.name.replaceFirst(".bigwig","")
  """
  #!/bin/bash
  computeMatrix reference-point \
              -S ${bigwig} \
              -R ${bed} \
              -o "${id}.${type}.deeptools.matrix" \
              --outFileNameMatrix "${id}.${type}.deeptools.tab" \
              --referencePoint center \
              -b 2000 \
              -a 2000

  plotHeatmap -m "${id}.${type}.deeptools.matrix" \
              --outFileName "${id}.${type}.deeptools.png" \
              --plotType se \
              --averageTypeSummaryPlot mean \
              --heatmapWidth 10 \
              --xAxisLabel "Distance" \
              --plotFileFormat png

  plotHeatmap -m "${id}.${type}.deeptools.matrix" \
              --outFileName "${id}.${type}.deeptools.svg" \
              --plotType se \
              --averageTypeSummaryPlot mean \
              --heatmapWidth 10 \
              --xAxisLabel "Distance" \
              --plotFileFormat svg
  """
  }

process computeMatrixRefPointMultiple {
  tag { "${bigwig}" }

  publishDir "${params.outdir}/deeptools/matrix",  mode: 'copy', overwrite: true, pattern: "*matrix"
  publishDir "${params.outdir}/deeptools/matrix",  mode: 'copy', overwrite: true, pattern: "*tab"
  publishDir "${params.outdir}/deeptools/heatmap", mode: 'copy', overwrite: true, pattern: "*svg"
  publishDir "${params.outdir}/deeptools/heatmap", mode: 'copy', overwrite: true, pattern: "*png"

  input:
  path(bigwig)
  path(bed)

  output:
  path('*.matrix', emit: matrix)
  path('*.tab',    emit: tab)
  path('*.png',    emit: png)
  path('*.svg',    emit: svg)

  script:
  def id = bigwig.name.replaceFirst(".bigwig","")
  def type = 'all'
  """
  #!/bin/bash
  computeMatrix reference-point \
              -S ${bigwig} \
              -R ${bed} \
              -o "${id}.${type}.deeptools.matrix" \
              --outFileNameMatrix "${id}.${type}.deeptools.tab" \
              --referencePoint center \
              -b 2000 \
              -a 2000

  plotHeatmap -m "${id}.${type}.deeptools.matrix" \
              --outFileName "${id}.${type}.deeptools.png" \
              --plotType se \
              --averageTypeSummaryPlot mean \
              --heatmapWidth 10 \
              --xAxisLabel "Distance" \
              --plotFileFormat png --perGroup

  plotHeatmap -m "${id}.${type}.deeptools.matrix" \
              --outFileName "${id}.${type}.deeptools.svg" \
              --plotType se \
              --averageTypeSummaryPlot mean \
              --heatmapWidth 10 \
              --xAxisLabel "Distance" \
              --plotFileFormat svg --perGroup
  """
  }
