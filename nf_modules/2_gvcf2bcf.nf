#!/usr/bin/env nextflow

process gvcf_2_bcf{
    container 'docker://erictdawson/glnexus'

    publishDir '3_bcf/', mode: 'copy'

    cpus 8
    memory '64 GB'

    input:
    path input_vcfs
    path regions

    output:
    path "Exomes.bcf"

    script:
    """
    glnexus_cli --config DeepVariant --bed ${regions} \
    ${input_vcfs} > Exomes.bcf
    """

}