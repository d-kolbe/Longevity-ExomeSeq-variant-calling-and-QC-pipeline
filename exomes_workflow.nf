//// Define parameters for nextflow
params.ref = "./parfiles/genome.fa" // reference genome
params.ref_ind = "./parfiles/genome.fa.fai" // reference genome index
params.regions = "./parfiles/hg38_exomregions_withoutCHR.bed" // exonic regions
params.plinkQC_rscript = "./parfiles/plinkQC.R" // quality control at plink level (hwe, ibd, het etc.)
params.vaf_script = "./parfiles/Annotate_VAF.py" // python script to annotate and process variant allele fraction (removal of somatic mutations)




// import modules
include { cram_2_vcf } from './nf_modules/1_cram2vcf.nf'
include { gvcf_2_bcf } from './nf_modules/2_gvcf2bcf.nf'
include { bcf_2_vcf } from './nf_modules/3_bcf2vcf.nf'
include { vcfQC } from './nf_modules/4_vcfQC.nf'
include { plinkQC } from './nf_modules/5_plinkQC.nf'
include { vafQC } from './nf_modules/6_vafQC.nf'


// final workflow
workflow {
    def cram = Channel.fromFilePairs("1_cram/*{.cram,.cram.crai}")
    cram_2_vcf(cram, file(params.ref), file(params.ref_ind), file(params.regions))

    def vcfs = cram_2_vcf.out.collect()
    gvcf_2_bcf(vcfs, file(params.regions))

    bcf_2_vcf(gvcf_2_bcf.out)

    vcfQC(bcf_2_vcf.out.vcf_merged, bcf_2_vcf.out.vcf_merged_index)
    plinkQC(vcfQC.out.vcfQC, vcfQC.out.vcfQC_index, file(params.plinkQC_rscript))

    vafQC(plinkQC.out.vcf_plinkQC, )

}

