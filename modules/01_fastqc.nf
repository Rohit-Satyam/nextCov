params.memory = "3g"
params.cpus = 1
params.output = "."


process FASTQC{
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/01_rawFastQC", mode: "copy"

    input:
        tuple val(sid), path(reads)

    output:
        path "*"

    """
    fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
    """
}
