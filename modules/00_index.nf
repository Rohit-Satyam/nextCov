params.memory = "3g"
params.cpus = 1
params.output = "."


process BWAINDEX{
        publishDir "${params.output}/00_index", mode: 'copy'
        memory params.memory
        cpus params.cpus

        input:
        path (fasta)

        output:
        path "$fasta", emit: bwa_idx
        path ("*")

        script:
        """
        bwa index $fasta
        gatk CreateSequenceDictionary -R $fasta
        samtools faidx $fasta
        """
        }
