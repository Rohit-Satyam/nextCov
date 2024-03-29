params.memory = "3g"
params.cpus = 1
params.output = "."
params.mode = "ibex"

process TRIMMOMATIC {
	cpus params.cpus
	memory params.memory
    publishDir "${params.output}/02_adapterTrimming", mode: 'copy'

input:
    tuple val(sid), path(reads)
output:
    tuple val(sid), path(fq_1_paired), path(fq_2_paired), emit:trimmed_reads
		path "*.log"
        script:
    fq_1_paired = sid + '_R1_P.fastq.gz'
    fq_1_unpaired = sid + '_R1_UP.fastq.gz'
    fq_2_paired = sid + '_R2_P.fastq.gz'
    fq_2_unpaired = sid + '_R2_UP.fastq.gz'
if ("${params.mode}" == "local")
        """
        trimmomatic PE \
        -threads ${task.cpus} \
        ${reads[0]} \
        ${reads[1]}\
        $fq_1_paired \
        $fq_1_unpaired \
        $fq_2_paired \
        $fq_2_unpaired \
        ILLUMINACLIP:${params.adapter}:2:30:10 \
        MINLEN:30 2> ${sid}.log
        """
else if ("${params.mode}" == "ibex")
        """
java -jar $TRIMMOMATIC_JAR PE -threads 20 ${reads[0]} ${reads[1]} \
        $fq_1_paired \
        $fq_1_unpaired \
        $fq_2_paired \
        $fq_2_unpaired \
        ILLUMINACLIP:${params.adapter}:2:30:10 \
        MINLEN:30 2> ${sid}.log
        """

}

process FASTP{
	cpus params.cpus
	memory params.memory
	publishDir "${params.output}/02_adapterTrimming", mode: 'copy'

input:
        tuple val(sid), path(reads)

    output:
        tuple val(sid), path(fq_1_paired), path(fq_2_paired)
		path("${sid}.fastp.json"), emit: fastp_logs
		path("${sid}.fastp.html")
		path "*.log"

        script:
    fq_1_paired = sid + '_R1.fastq.gz'
    fq_2_paired = sid + '_R2.fastq.gz'
    """
	fastp \
	--in1 ${reads[0]} \
	--in2 ${reads[1]}\
    --thread ${task.cpus} \
	--out1 $fq_1_paired \
	--out2 $fq_2_paired \
	--json ${sid}.fastp.json \
	--html ${sid}.fastp.html \
    ${params.fastp_ext} \
    2> ${sid}.log

    """
}

process POSTTRIMFASTQC{
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/02_adapterTrimming/postTrimFASTQC", mode: 'copy'

    input:
        tuple val(sid), path(reads1), path(reads2)

    output:
        path "*", emit: postfastqc

    """
    fastqc -t ${task.cpus} ${reads1} ${reads2}
    """
}
