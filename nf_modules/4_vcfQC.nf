#!/usr/bin/env nextflow

process vcfQC{
    container 'docker://erictdawson/glnexus'

    publishDir '5_vcfQC/', mode: 'copy'

    cpus 8
    memory '64 GB'

    input:
    path input_vcf
    path input_vcf_index


    output:
    path "ExomesQC1.vcf.gz", emit: vcfQC
    path "ExomesQC1.vcf.gz.csi", emit: vcfQC_index


    script:
    """
    ### Intial QC steps at VCF level (remove half calls, genotype quality < 20, missingness >= 0.05, split biallelic sites)

    bcftools view ${input_vcf} -Ov | awk -F'\t' -v OFS='\t' '{ if(\$7==".") \$7="PASS"; print }' - | \
    bcftools +setGT - -Ov -- -t q -n . -i 'FMT/GQ<20' | \
    bcftools view - -m2 -M3 -e 'GT="./0" | GT="0/." | GT="./1" | GT="./2" | GT="1/." | GT="2/."' -Ov | \
    bcftools view - -a --threads 8 -Ov | bcftools norm - -m-both -Oz -o ExomesQC1.vcf.gz

    bcftools index ExomesQC1.vcf.gz


    """

}