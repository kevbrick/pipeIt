profiles {
  singularity {
    singularity.enabled = true
    singularity.autoMounts = true
    singularity.envWhitelist='https_proxy,http_proxy,ftp_proxy,DISPLAY,SLURM_JOBID,GENOMES,NXF_GENOMES'
    singularity.runOptions = ' -B /data/RDCO -B /data/brickkm -B /lscratch -B /fdb/igenomes/ '
    process.container = "$DSL2DIR/singularity/RDCO.sif"
  }

  local {
    process.maxForks = 1
    process.executor='local'
    env{
      TMPDIR='./'
    }
  }

  slurm {
    process.executor='slurm'
    process.container = "$DSL2DIR/singularity/RDCO.sif"
    process.scratch = '/lscratch/$SLURM_JOBID'
    process.clusterOptions = ' --gres=lscratch:800 --partition=norm'
    env{
      TMPDIR='/lscratch/$SLURM_JOBID'
    }
  }

  none {
    // Add custom configs here
  }
}

params.outdir = './nxfOut'

report {
  enabled = true
  file = "${params.outdir}/nxfReports/${params.name}.report.html"
}

timeline {
  enabled = true
  file = "${params.outdir}/nxfReports/${params.name}.timeline.html"
}

trace {
  enabled = true
  file = "${params.outdir}/nxfReports/${params.name}.trace.txt"
}

manifest {
  description = '2020: Kevin Brick'
}

params.igenomes_base = '/fdb/igenomes/'
includeConfig "igenomes.config.nf"
