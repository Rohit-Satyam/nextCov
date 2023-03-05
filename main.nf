#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from './modules/01_fastqc'
include {TRIMMOMATIC; FASTP; POSTTRIMFASTQC} from './modules/02_trimmomatic'
include {BWAMEM} from './modules/03_align'
include {IVAR} from './modules/04_ivar'
include {MARKDUP} from './modules/05_markduplicates'
include {BAMSTATS} from './modules/06_bamstats'
include {GATK; BCFTOOLS; LOFREQ; IVAR_VARIANT_CALL; IVAR_TO_VCF } from './modules/07_variantcall'
include {VAFATOR; PHASING} from './modules/08_variantAnnotation'
include {MULTIQC as PRETRIM; MULTIQC as POSTTRIM; MULTIQC as SUMMARISEALL } from './modules/01_fastqc'

params.help= false
params.input = false
params.output= false
params.mode = false
if (params.input != false) {
            Channel.fromFilePairs(params.input, checkIfExists: true )
                .set { input_fastqs }
        }


workflow{

//fastqc
if (input_fastqs) {
		rawfqc_ch=FASTQC(input_fastqs)
    PRETRIM("01_rawFastQC",FASTQC.out.fastqc.collect(),'pre-trimming')
  }

//trimming
if (params.skipTrim) {
  trim_ch=rawfqc_ch
} else {
  FASTP(input_fastqs)
  trim_ch=FASTP.out[0]
  POSTTRIMFASTQC(trim_ch)
  postrim_input=POSTTRIMFASTQC.out.postfastqc.collect()
	POSTTRIM("02_adapterTrimming",postrim_input,'post-trimming')
  }


//alignment,primer trimming and mark duplicates
if (params.skipPrimerTrim){
    BWAMEM(trim_ch, params.ref)
    MARKDUP(BWAMEM.out[0])
} else {
    BWAMEM(trim_ch, params.ref)
    IVAR(BWAMEM.out[0])
    MARKDUP(IVAR.out[0])}
//BAM Stats
    BAMSTATS(MARKDUP.out[0])
//Variant calling and Normalization
    GATK(MARKDUP.out[0],params.ref)
    BCFTOOLS(MARKDUP.out[0],params.ref)
    LOFREQ(MARKDUP.out[0],params.ref)
    IVAR_VARIANT_CALL(MARKDUP.out[0],params.ref,params.gff)
    IVAR_TO_VCF(IVAR_VARIANT_CALL.out[0],params.ref)

//Variant Allele freq 08_variantAnnotation
    normVcf_ch= GATK.out[0].mix(BCFTOOLS.out[0],LOFREQ.out[0],IVAR_TO_VCF.out[0])
    VAFATOR(normVcf_ch)
    PHASING(VAFATOR.out[0],params.ref,params.gff)

if (params.skipTrim) {
	all_combine=BAMSTATS.out[0],mix(GATK.out.gatkRecalTable.collect(),MARKDUP.out.dedupmtx.collect())
	SUMMARISEALL("Summary_Reports",all_combine.collect(),'summary_report')
} else {
	all_combine=BAMSTATS.out[0].mix(GATK.out.gatkRecalTable.collect(), FASTP.out.fastp_logs.collect(),MARKDUP.out.dedupmtx.collect())
	SUMMARISEALL("Summary_Reports",all_combine.collect(),'summary_report')}

}
