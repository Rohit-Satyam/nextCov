#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {BWAINDEX} from './modules/00_index'
include { FASTQC } from './modules/01_fastqc'
include {TRIMMOMATIC; FASTP; POSTTRIMFASTQC} from './modules/02_trimmomatic'
include {BWAMEM} from './modules/03_align'
include {IVAR} from './modules/04_ivar'
include {MARKDUP} from './modules/05_markduplicates'
include {BAMSTATS} from './modules/06_bamstats'
include {GATK; BCFTOOLS; LOFREQ; IVAR_VARIANT_CALL; IVAR_TO_VCF } from './modules/07_variantcall'
include {VAFATOR; PHASING} from './modules/08_variantAnnotation'
params.help= false
params.input = false
params.output= false
params.mode = false
if (params.input != false) {
            Channel.fromFilePairs(params.input, checkIfExists: true )
                .set { input_fastqs }
        }


workflow{
  //indexing
  if (params.index){
    BWAINDEX(params.ref)
  } else {
    reference_ch=params.ref
  }
//fastqc
	if (input_fastqs) {
		rawfqc_ch=FASTQC(input_fastqs)
//trimming
    if (params.trimUse == "trimmomatic"){
		TRIMMOMATIC(input_fastqs)
    trim_ch=TRIMMOMATIC.out[0]
    POSTTRIMFASTQC(trim_ch)
    } else if (params.trimUse == "fastp") {
		FASTP(input_fastqs)
    trim_ch=FASTP.out[0]
    POSTTRIMFASTQC(trim_ch)
    }
//alignment,primer trimming and mark duplicates
    if (params.primerTrim){
    BWAMEM(trim_ch, BWAINDEX.out[0])
    IVAR(BWAMEM.out[0])
    MARKDUP(IVAR.out[0])
    } else {
    BWAMEM(trim_ch, BWAINDEX.out[0])
    MARKDUP(BWAMEM.out[0])
    }
//BAM Stats
    BAMSTATS(MARKDUP.out[0])
//Variant calling and Normalization
    GATK(MARKDUP.out[0],BWAINDEX.out[0])
    BCFTOOLS(MARKDUP.out[0],BWAINDEX.out[0])
    LOFREQ(MARKDUP.out[0],BWAINDEX.out[0])
    IVAR_VARIANT_CALL(MARKDUP.out[0],BWAINDEX.out[0],params.gff)
    IVAR_TO_VCF(IVAR_VARIANT_CALL.out[0],BWAINDEX.out[0])

//Variant Allele freq 08_variantAnnotation
    normVcf_ch= GATK.out[0].mix(BCFTOOLS.out[0],LOFREQ.out[0],IVAR_TO_VCF.out[0])
    normVcf_ch.view()
    VAFATOR(normVcf_ch)
    PHASING(VAFATOR.out[0],BWAINDEX.out[0],params.gff)
	}
}
