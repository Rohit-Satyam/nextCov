
params.memory = "3g"
params.cpus = 1
params.output = "."


process BWAMEM{
        publishDir "${params.output}/03_alignment", mode: 'copy'
        memory params.memory
        cpus params.cpus

        input:
        tuple val(sid), file(reads1), file(reads2)
        path(reference)

        output:
        tuple val(sid), path ("${sid}.sorted.bam"), emit: alignments
        path "${sid}.sorted.bam.bai"

        shell:
        '''
        id=$(zcat !{reads1} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//')
        bwa mem -M -R "$(echo "@RG\\tID:${id}\\tSM:!{sid}\\tPL:ILLUMINA")" -t !{task.cpus} !{reference.toRealPath()} !{reads1} !{reads2} | samtools sort -@ !{task.cpus} -o !{sid}.sorted.bam -
        samtools index -@ !{task.cpus} !{sid}.sorted.bam
        '''
}
