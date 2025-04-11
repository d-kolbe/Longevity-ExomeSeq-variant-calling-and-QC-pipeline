#!/usr/bin/env nextflow

process cram_2_vcf{
    container 'docker://google/deepvariant'

    publishDir '2_gvcf/', mode: 'copy'

    cpus 8
    memory '64 GB'

    input:
    tuple val(sample),path(crams)
    path ref
    path ref_ind
    path regions

    output:
    path "${sample}.g.vcf.gz"

    script:
    def (cram_in, crai_in) = crams

    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type WES \
    --ref ${ref} \
    --reads ${cram_in} \
    --regions ${regions} \
    --output_gvcf ${sample}.g.vcf.gz \
    --output_vcf ${sample}.vcf.gz \
    --num_shards 8 \

    """
}
