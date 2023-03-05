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
        path "*", emit: fastqc

    """
    fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
    """
}

process MULTIQC {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/${name}", mode: "copy"
    input:
    val(name)
    path(filepaths)
    val(filename)
    

    output:
    path "*"

    script:
    """
    multiqc --force --config ${projectDir}/bin/multiqc_config.yaml --filename ${filename} ${filepaths} 
    """
}