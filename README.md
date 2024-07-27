# nextCov

A short read Nextflow pipeline for SARS-Cov2 variant calling and assembly development

## Environment Setup

```
mamba create -n nextcov -c conda-forge -c bioconda \
assembly-stats  \
bcftools \
biopython \
bwa-mem2 \
covtobed  \
fastp \
fastqc  \
fonttools \
gatk4  \
ghostscript \
gtfparse=1.2.1 \
intervene \
ivar  \
lofreq \
matplotlib \
mosdepth  \
multiqc  \
mummer  \
nextclade  \
nextflow  \
picard \
pybedtools \
rtg-tools \
r-upsetr \
r-gplots \
r-gsalib \
samplot \
samtools  \
scikit-learn \
seaborn \
snpeff \
trimmomatic  \
vafator  \
openjdk==8.0.332=h166bdaf_0
```

## Running the pipeline

The one liner command to produce VCF are given below.

```
nextflow run main.nf --input "test/*R{1,2}_001.fastq.gz"
```

## TO do's

1.Impact of using common variants from all variant callers.
2. Find right combination of intersection since IVAR and LOFreq determine low frequency variant and GATK and bcftools capture high frequency variants.
