params.memory = "3g"
params.cpus = 1
params.output = "."
params.jobs = 5


process BAMSTATS{
  publishDir "${params.output}/06_bamStats", mode: 'copy'
  memory params.memory
  cpus params.cpus
  maxForks params.jobs
        input:
        tuple val(sid), path(bam)
        output:
        path "*"

        script:
        """
        samtools flagstat -@ ${task.cpus} ${bam} > ${sid}.flagstat.tsv
        samtools coverage ${bam} > ${sid}.coverage.txt
      	samtools depth -@ ${task.cpus} -s -d 0 -H ${bam} > ${sid}.depth.tsv
      	mosdepth -t ${task.cpus} -x ${sid} ${bam.toRealPath()}
        """
}

process MASKASSEMBLY{

  script:
  """
  covtobed -x ${params.coverageThreshold} ${bam} > ${sid}_lowcoverage.bed
  subtractBed -a ${sid}_lowcoverage.bed -b ${params.exclude} | mergeBed  > ${sid}.mask.bed
  """
}
