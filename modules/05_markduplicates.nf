params.memory = "3g"
params.cpus = 1
params.output = "."
params.jobs = 5


process MARKDUP{
  publishDir "${params.output}/05_markDuplicates", mode: 'copy'
  memory params.memory
  cpus params.cpus
  maxForks params.jobs
        input:
        tuple val(sid), path(bam)
        output:
        tuple val(sid), path("${sid}.dedup.bam")
        path("${sid}_markdup_metrics.txt"), emit: dedupmtx
        path("${sid}.dedup.bam.bai")
        """
        gatk MarkDuplicatesSpark -I ${bam} -O ${sid}.dedup.bam \
        -M ${sid}_markdup_metrics.txt --create-output-bam-index true --optical-duplicate-pixel-distance ${params.opticalDupPixalDis}  --tmp-dir .
        """
}
