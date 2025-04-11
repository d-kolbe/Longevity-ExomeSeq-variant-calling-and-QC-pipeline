#!/usr/bin/env nextflow

process vafQC{
    container 'docker://erictdawson/glnexus'

    publishDir '7_vafQC/', mode: 'copy'

    cpus 8
    memory '64 GB'

    input:
    path input_vcf
    path input_vcf_index
    path vaf_script


    output:
    path "Exomes_chr*.vcf.gz", emit: vcf_single

    script:
    """

    ### Remove variants with lower or higher than 1st and 99th percentile of VAF (i.e., variants that are likely somatic mutations – see CHIP)

    bcftools query -f'[%CHROM\t%POS\t%REF\t%ALT\t%ID\t%SAMPLE\t%GT\t%AD\n]' ExomesQC2.vcf.gz | pigz -p 8 > ExomesQC2_query.tsv
    python ${vaf_script} ExomesQC2_query.tsv ExomesQC2_snps2keep.txt

    bcftools view --threads 8 -R ExomesQC2_snps2keep.txt ExomesQC2.vcf.gz -Oz -o ExomesQC3.vcf.gz
    bcftools index ExomesQC3.vcf.gz



    ### Then split files by single autosomal chromosomes – these are now ready to be converted to .gds file (see STAARpipeline by X. Li)

    chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

    for chromosome in \$(seq 1 22); do
        output_file="Exomes_chr\${chromosome}.vcf.gz"
        bcftools view -f PASS -r "\$chromosome" ExomesQC3.vcf.gz -Oz -o \$output_file
    done
    """

}


    