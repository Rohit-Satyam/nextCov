params.memory = "3g"
params.cpus = 1
params.low_frequency_variant_threshold = 0.2
params.subclonal_variant_threshold = 0.8


process VAFATOR {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/08_variantAnnotation/unphased/${caller}", pattern: "${sid}.${caller}.vaf.annot.{vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'

    input:
    tuple val(sid), val(caller), path(vcf), path(bam)

    output:
    tuple val(sid), val(caller), path("${sid}.${caller}.vaf.annot.vcf.gz"), emit: annotated_vcf
    path "${sid}.${caller}.vaf.annot.vcf.gz.tbi"
    path "${sid}.${caller}.vaf.annot.vcf.gz.csi"
    script:
    """
    vafator \
    --input-vcf ${vcf} \
    --output-vcf ${sid}.${caller}.vaf.vcf \
    --bam vafator ${bam.toRealPath()} --mapping-quality 0 --base-call-quality 0

    bgzip -c ${sid}.${caller}.vaf.vcf > ${sid}.${caller}.vaf.vcf.gz
    tabix -p vcf ${sid}.${caller}.vaf.vcf.gz

    # annotates low frequency and subclonal variants
    #bcftools view -Ob ${sid}.${caller}.vaf.vcf.gz | \
    #bcftools filter \
    #--exclude 'INFO/vafator_af < ${params.low_frequency_variant_threshold}' \
    #--soft-filter LOW_FREQUENCY - | \
    #bcftools filter \
    #--exclude 'INFO/vafator_af >= ${params.low_frequency_variant_threshold} && INFO/vafator_af < ${params.subclonal_variant_threshold}' \
    #--soft-filter SUBCLONAL \
    #--output-type v - > ${sid}.${caller}.vaf.annot.vcf

    bcftools view -Ob ${sid}.${caller}.vaf.vcf.gz   | \
    bcftools filter --exclude 'INFO/vafator_af < 0.5 && INFO/vafator_dp < 100 && INFO/vafator_ac < 50' \
    --soft-filter POOR_CALLS --output-type v - | bcftools norm -d both - > ${sid}.${caller}.vaf.annot.vcf

    bgzip -c ${sid}.${caller}.vaf.annot.vcf > ${sid}.${caller}.vaf.annot.vcf.gz
    tabix -p vcf ${sid}.${caller}.vaf.annot.vcf.gz
    bcftools index ${sid}.${caller}.vaf.annot.vcf.gz
    """
}

process PHASING {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/08_variantAnnotation/phased/${caller}", pattern: "${sid}.${caller}.phased.{vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'


    input:
        tuple val(sid), val(caller), path(vcf)
        val(fasta)
        val(gtf)

    output:
    tuple val(sid), val(caller), path("${sid}.${caller}.phased.vcf.gz"), emit: annotated_vcf
    path "${sid}.${caller}.phased.vcf.gz.tbi"
    path "${sid}.${caller}.phased.vcf.gz.csi"

    script:
    """
    python $projectDir/modules/phasing.py \
    --fasta ${fasta} \
    --gtf ${gtf} \
    --input-vcf ${vcf} \
    --output-vcf ${sid}.${caller}.phased.vcf
    bgzip -c ${sid}.${caller}.phased.vcf > ${sid}.${caller}.phased.vcf.gz
    tabix -p vcf ${sid}.${caller}.phased.vcf.gz
    bcftools index ${sid}.${caller}.phased.vcf.gz
    """
}

process FILTERING {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}/08_variantAnnotation/phased/${caller}", pattern: "${sid}.${caller}.phased.{vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'


    input:
        tuple val(sid), val(caller), path(vcf)
        val(fasta)
        val(gtf)

    output:
    tuple val(sid), val(caller), path("${sid}.${caller}.phased.vcf.gz"), emit: annotated_vcf
    path "${sid}.${caller}.phased.vcf.gz.tbi"
    path "${sid}.${caller}.phased.vcf.gz.csi"

    script:
    """
    python $projectDir/modules/phasing.py \
    --fasta ${fasta} \
    --gtf ${gtf} \
    --input-vcf ${vcf} \
    --output-vcf ${sid}.${caller}.phased.vcf
    bgzip -c ${sid}.${caller}.phased.vcf > ${sid}.${caller}.phased.vcf.gz
    tabix -p vcf ${sid}.${caller}.phased.vcf.gz
    bcftools index ${sid}.${caller}.phased.vcf.gz
    """
}
