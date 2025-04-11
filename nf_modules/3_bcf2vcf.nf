#!/usr/bin/env nextflow

process bcf_2_vcf{
    container 'docker://erictdawson/glnexus'

    publishDir '4_vcf/', mode: 'copy'

    cpus 8
    memory '64 GB'

    input:
    path input_bcf

    output:
    path "Exomes.vcf.gz", emit: vcf_merged
    path "Exomes.vcf.gz.csi", emit: vcf_merged_index

    script:
    """
    bcftools view ${input_bcf}| bgzip -@ 4 -c > Exomes.vcf.gz
    bcftools index Exomes.vcf.gz
    """

}