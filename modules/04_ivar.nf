params.memory = "3g"
params.cpus = 1
params.output = "."

process IVAR{
        publishDir "${params.output}/04_primerTrimmed_alignments", mode: 'copy'
        memory params.memory
        cpus params.cpus

        input:
        tuple val(sid), path(bam)
        output:
        tuple val(sid), path("${sid}.ivar_sorted.bam"), emit: ivar_trimmed_bams
        path("${sid}.ivar_sorted.bam.bai")
        script:
        """
        ivar trim -i ${bam} -b ${params.primers} -p ${sid}.ivar.bam
        samtools sort -@ ${task.cpus} -o ${sid}.ivar_sorted.bam ${sid}.ivar.bam
        samtools index -@ ${task.cpus} ${sid}.ivar_sorted.bam
        """
}
