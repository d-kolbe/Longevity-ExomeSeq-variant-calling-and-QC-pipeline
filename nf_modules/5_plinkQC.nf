#!/usr/bin/env nextflow

process plinkQC{
    container 'docker://zilinli/staarpipeline:0.9.7'

    publishDir '6_plinkQC/', mode: 'copy'

    cpus 4
    memory '32 GB'

    input:
    path input_vcf
    path input_vcf_index
    path plinkQC_rscript

    output:
    path "ExomesQC2.vcf.gz", emit: vcf_plinkQC
    path "ExomesQC2.vcf.gz.csi", emit: vcf_plinkQC_index

    script:
    """
    ### First lets  make the plink files
    plink --vcf ${input_vcf} --const-fid "0" --make-bed --out Exomes
    Rscript ${plinkQC_rscript}


    bcftools view ${input_vcf} -T ^snps2remove.txt -S ^individuals2remove.txt -Oz -o ExomesQC2.vcf.gz
    bcftools index ExomesQC2.vcf.gz



    """

}