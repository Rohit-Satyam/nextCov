params.memory = "3g"
params.cpus = 1
params.output = "."
params.jobs = 5
params.min_mapping_quality = 20
params.min_base_quality = 20

process GATK{
	if (params.runbqsr){
	publishDir "${params.output}/07_variantcall/gatk/intermediate", pattern: '*.vcf', mode: 'copy'
	publishDir "${params.output}/07_variantcall/gatk/bqsrfiles",pattern: "*.table", mode: 'copy'
	publishDir "${params.output}/07_variantcall/gatk/bqsrfiles",pattern: "*.pdf", mode: 'copy'
	publishDir "${params.output}/07_variantcall/gatk/raw_vcf",pattern: "*.raw.gatk.{bcf,vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'
	publishDir "${params.output}/07_variantcall/gatk/normalized_vcf", pattern: '*.normalized.gatk.{vcf.gz,vcf.gz.tbi,vcf.gz.csi}', mode: 'copy'
	publishDir "${params.output}/07_variantcall/gatk/log", pattern: '*.log', mode: 'copy'
} else {
	publishDir "${params.output}/07_variantcall/gatk/raw_vcf",pattern: "*.raw.gatk.{bcf,vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'
	publishDir "${params.output}/07_variantcall/gatk/normalized_vcf", pattern: "*.normalized.gatk.{vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'
	publishDir "${params.output}/07_variantcall/gatk/log", pattern: '*.log', mode: 'copy'
}

memory params.memory
cpus params.cpus
maxForks params.jobs

		input:
		tuple val(sid), path(bam)
		path(reference)

		output:
		tuple val(sid), val("gatk"), path("${sid}.normalized.gatk.vcf.gz"), path(bam)
		path "*.bcf"
		path "*.tbi"
		path "*.csi"
		path "*.vcf"
		path "*.pdf"
		path "*.table", emit: gatkRecalTable
		path "*.log"

		script:
if ("${params.runbqsr}")
		"""
		gatk HaplotypeCaller --sample-name ${sid} -ploidy 1 \
		--min-base-quality-score ${params.min_base_quality} --minimum-mapping-quality ${params.min_mapping_quality} \
		--native-pair-hmm-threads 8 -R ${reference.toRealPath()} --annotation AlleleFraction -I ${bam} \
		-O ${sid}.raw.vcf

		gatk VariantFiltration -R ${reference.toRealPath()} -V ${sid}.raw.vcf \
		--filter-expression \"QD<2.0 || FS > 60.0 || SOR > 4.0 || ReadPosRankSum < -8.0\"  -filter-name  \"good_filter\" \
		-O ${sid}.filtered_variants.vcf

		gatk SelectVariants -V ${sid}.filtered_variants.vcf --select-type-to-include SNP   -O ${sid}_known.vcf

		gatk BaseRecalibrator -R ${reference.toRealPath()} -I ${sid}.dedup.bam --known-sites  ${sid}_known.vcf  \
		-O ${sid}_beforerecal_data.table

		gatk ApplyBQSR -R ${reference.toRealPath()} -I ${bam} --bqsr-recal-file ${sid}_beforerecal_data.table \
		-O ${sid}_recal.bam

		gatk BaseRecalibrator -R ${reference.toRealPath()} -I ${sid}_recal.bam --known-sites ${sid}_known.vcf \
		-O ${sid}_afterrecal_data.table

		gatk AnalyzeCovariates -before ${sid}_beforerecal_data.table -after ${sid}_afterrecal_data.table \
		-plots ${sid}_AnalyzeCovariates.pdf

		gatk HaplotypeCaller --sample-name ${sid} -ploidy 1 \
		--min-base-quality-score ${params.min_base_quality} --minimum-mapping-quality ${params.min_mapping_quality} \
		--native-pair-hmm-threads 8 -R ${reference.toRealPath()} --annotation AlleleFraction \
		-I ${sid}_recal.bam -O ${sid}.raw.gatk.vcf 2>${sid}.raw.gatk.vcf.log
		
		bcftools view --output-type b ${sid}.raw.gatk.vcf > ${sid}.raw.gatk.bcf
		bgzip -c ${sid}.raw.gatk.vcf > ${sid}.raw.gatk.vcf.gz
		tabix -p vcf ${sid}.raw.gatk.vcf.gz
		bcftools index ${sid}.raw.gatk.vcf.gz

		bcftools sort ${sid}.raw.gatk.vcf.gz |  bcftools norm \
		--multiallelics -any --check-ref e --fasta-ref ${reference.toRealPath()} --old-rec-tag OLD_CLUMPED --atomize - | \
		bcftools norm --rm-dup exact --output-type z -o ${sid}.normalized.gatk.vcf.gz -

		tabix -p vcf ${sid}.normalized.gatk.vcf.gz
		bcftools index ${sid}.normalized.gatk.vcf.gz


		"""
else
		"""
		touch dummy.table
		gatk HaplotypeCaller --sample-name ${sid} -ploidy 1 \
		--min-base-quality-score ${params.min_base_quality} --minimum-mapping-quality ${params.min_mapping_quality} \
		--native-pair-hmm-threads 8 -R ${reference.toRealPath()} --annotation AlleleFraction -I ${bam} \
		-O ${sid}.raw.gatk.vcf 2>${sid}.raw.gatk.vcf.log
		
		bcftools view --output-type b ${sid}.raw.gatk.vcf > ${sid}.raw.gatk.bcf
		bgzip -c ${sid}.raw.gatk.vcf > ${sid}.raw.gatk.vcf.gz
		tabix -p vcf ${sid}.raw.gatk.vcf.gz
		bcftools index ${sid}.raw.gatk.vcf.gz

		bcftools sort ${sid}.raw.gatk.vcf.gz |  bcftools norm \
		--multiallelics -any --check-ref e --fasta-ref ${reference.toRealPath()} --old-rec-tag OLD_CLUMPED \
		--atomize - | \
		bcftools norm --rm-dup exact --output-type z -o ${sid}.normalized.gatk.vcf.gz -

		tabix -p vcf ${sid}.normalized.gatk.vcf.gz
		bcftools index ${sid}.normalized.gatk.vcf.gz

		"""
}


process BCFTOOLS {
	publishDir "${params.output}/07_variantcall/bcftools/raw_vcf", pattern: "*.raw.bcftools.{bcf,vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'
	publishDir "${params.output}/07_variantcall/bcftools/normalized_vcf", pattern: "*.normalized.bcftools.{vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'
	publishDir "${params.output}/07_variantcall/bcftools/log", pattern: '*.log', mode: 'copy'
	memory params.memory
	cpus params.cpus
	maxForks params.jobs

    input:
        tuple val(sid), path(bam)
        path(reference)

    output:
		tuple val(sid), val("bcftools"), path("${sid}.normalized.bcftools.vcf.gz"), path(bam)
		path "*.bcf"
		path "*.tbi"
		path "*.csi"
		path "*.log"

		script:
if ("${params.experimentalBCF}")
		"""
		bcftools mpileup \
		--redo-BAQ \
		--max-depth 0 \
		--min-BQ ${params.min_base_quality} \
		--min-MQ ${params.min_mapping_quality} \
		--count-orphans \
		--ignore-overlaps \
		--min-ireads 10 \
		--tandem-qual 500\
		--fasta-ref ${reference.toRealPath()} \
		--annotate "AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR" ${bam} | \
		bcftools call \
		--multiallelic-caller \
		--variants-only \
		--ploidy 1 \
		--output-type b - > ${sid}.raw.bcftools.bcf 2>${sid}.raw.bcftools.bcf.log

	bcftools convert -O z -o ${sid}.raw.bcftools.vcf.gz ${sid}.raw.bcftools.bcf
	tabix -p vcf ${sid}.raw.bcftools.vcf.gz
	bcftools index ${sid}.raw.bcftools.vcf.gz

	bcftools sort ${sid}.raw.bcftools.vcf.gz |  bcftools norm \
	--multiallelics -any --check-ref e --fasta-ref ${reference.toRealPath()} --old-rec-tag OLD_CLUMPED --atomize - | \
	bcftools norm --rm-dup exact --output-type z -o ${sid}.normalized.bcftools.vcf.gz -

	tabix -p vcf ${sid}.normalized.bcftools.vcf.gz
	bcftools index ${sid}.normalized.bcftools.vcf.gz
		"""
else
    """
    bcftools mpileup \
    --redo-BAQ \
    --max-depth 0 \
    --min-BQ ${params.min_base_quality} \
    --min-MQ ${params.min_mapping_quality} \
    --count-orphans \
    --fasta-ref ${reference.toRealPath()} \
    --annotate AD ${bam} | \
	bcftools call \
    --multiallelic-caller \
    --variants-only \
    --ploidy 1 \
    --output-type b - > ${sid}.raw.bcftools.bcf 2>${sid}.raw.bcftools.bcf.log

	bcftools convert -O z -o ${sid}.raw.bcftools.vcf.gz ${sid}.raw.bcftools.bcf
	tabix -p vcf ${sid}.raw.bcftools.vcf.gz
	bcftools index ${sid}.raw.bcftools.vcf.gz

	bcftools sort ${sid}.raw.bcftools.vcf.gz |  bcftools norm \
	--multiallelics -any --check-ref e --fasta-ref ${reference.toRealPath()} --old-rec-tag OLD_CLUMPED --atomize - | \
	bcftools norm --rm-dup exact --output-type z -o ${sid}.normalized.bcftools.vcf.gz -

	tabix -p vcf ${sid}.normalized.bcftools.vcf.gz
	bcftools index ${sid}.normalized.bcftools.vcf.gz
    """
}

process LOFREQ {
	publishDir "${params.output}/07_variantcall/lofreq/raw_vcf", pattern: "*.raw.lofreq.{bcf,vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'
	publishDir "${params.output}/07_variantcall/lofreq/normalized_vcf", pattern: "*.normalized.lofreq.{vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'
	publishDir "${params.output}/07_variantcall/lofreq/log", pattern: '*.log', mode: 'copy'

	memory params.memory
	cpus params.cpus
	maxForks params.jobs
    input:
        tuple val(sid), path(bam)
        path(reference)

    output:
		tuple val(sid), val("lofreq"), path("${sid}.normalized.lofreq.vcf.gz"), path(bam)
		path "*.bcf"
		path "*.tbi"
		path "*.csi"
		path "*.log"

		script:
    """
    lofreq call \
    --min-bq ${params.min_base_quality} \
    --min-alt-bq ${params.min_base_quality} \
    --min-mq ${params.min_mapping_quality} \
    --ref ${reference.toRealPath()} \
    --call-indels \
    <( lofreq indelqual --dindel --ref ${reference.toRealPath()} ${bam} ) | bgzip > ${sid}.raw.lofreq.vcf.gz 2>${sid}.raw.lofreq.vcf.log
    
	# NOTE: adding the tabix index is a dirty fix to deal with LoFreq VCF missing the chromosome in the header
    bcftools index ${sid}.raw.lofreq.vcf.gz
	tabix -p vcf ${sid}.raw.lofreq.vcf.gz
    bcftools view --output-type b ${sid}.raw.lofreq.vcf.gz > ${sid}.raw.lofreq.bcf

	bcftools sort ${sid}.raw.lofreq.vcf.gz |  bcftools norm \
	--multiallelics -any --check-ref e --fasta-ref ${reference.toRealPath()} --old-rec-tag OLD_CLUMPED --atomize - | \
	bcftools norm --rm-dup exact --output-type z -o ${sid}.normalized.lofreq.vcf.gz -

	tabix -p vcf ${sid}.normalized.lofreq.vcf.gz
	bcftools index ${sid}.normalized.lofreq.vcf.gz
    """
}


process IVAR_VARIANT_CALL {
	publishDir "${params.output}/07_variantcall/ivar/intermediate", pattern: '*.tsv', mode: 'copy'
	publishDir "${params.output}/07_variantcall/ivar/log", pattern: '*.log', mode: 'copy'
	memory params.memory
	cpus params.cpus
	maxForks params.jobs
    input:
        tuple val(sid), path(bam)
        path(reference)
        path(gff)

    output:
		tuple val(sid), path("${sid}.ivar.tsv"), path(bam)
		path "*.log"

		script:
    """
    samtools mpileup \
    -aa \
    --count-orphans \
    --max-depth 0 \
	--reference ${reference.toRealPath()} \
    --redo-BAQ -x \
    --min-BQ ${params.min_base_quality} \
    --min-MQ ${params.min_mapping_quality} \
    ${bam} | \
	ivar variants \
    -p ${sid}.ivar \
    -q ${params.min_base_quality} \
    -t 0.03 \
    -r ${reference.toRealPath()} \
    -g ${gff} 2>${sid}.raw.ivar.tsv.log
    """
}

process IVAR_TO_VCF{
	publishDir "${params.output}/07_variantcall/ivar/raw_vcf", pattern: "*.raw.ivar.{bcf,vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'
	publishDir "${params.output}/07_variantcall/ivar/normalized_vcf", pattern: "*normalized.ivar.{vcf.gz,vcf.gz.tbi,vcf.gz.csi}", mode: 'copy'

	memory params.memory
	cpus params.cpus
	maxForks params.jobs

	input:
	tuple val(sid), path(tsv),	path(bam)
	path(reference)


	output:
	tuple val(sid), val("ivar"), path("${sid}.normalized.ivar.vcf.gz"), path(bam)
	path "*.bcf"
	path "*.tbi"
	path "*.csi"

	script:
	"""

	python $projectDir/bin/covigatorivar2vcf.py --ivar ${tsv} --fasta ${reference.toRealPath()} \
	--output-vcf ${sid}.raw.ivar.vcf

	bcftools sort -O v --output-file sorted_${sid}.raw.ivar.vcf ${sid}.raw.ivar.vcf
	bcftools view --output-type b sorted_${sid}.raw.ivar.vcf > ${sid}.raw.ivar.bcf
	bgzip -c sorted_${sid}.raw.ivar.vcf > ${sid}.raw.ivar.vcf.gz
	tabix -p vcf ${sid}.raw.ivar.vcf.gz
	bcftools index ${sid}.raw.ivar.vcf.gz


	bcftools norm \
	--multiallelics -any --check-ref e --fasta-ref ${reference.toRealPath()} --old-rec-tag OLD_CLUMPED --atomize ${sid}.raw.ivar.vcf.gz | \
	bcftools norm --rm-dup exact --output-type z -o ${sid}.normalized.ivar.vcf.gz -

	tabix -p vcf ${sid}.normalized.ivar.vcf.gz
	bcftools index ${sid}.normalized.ivar.vcf.gz
	"""
}
